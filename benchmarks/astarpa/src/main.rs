//! The A*PA side of the comparison: global edit distance on the DNA workloads only, since A*PA
//! supports nothing else.
//!
//! Each aligner runs twice per pair: with traceback, which returns a CIGAR, and without, which
//! returns only the cost, so both of dinara-align's tasks have a like-for-like row.
//!
//! Each measurement is taken warm, in this process, the way `ours.mojo` times dinara-align's
//! bit-parallel rows, so microsecond workloads compare like for like (see `measure`).

use astarpa2::AstarPa2Params;
use pa_heuristic::Prune;
use std::{env, fs, time::Instant};

/// Batches whose average is reported: the fastest of them.
const BATCHES: usize = 20;
/// About how long one batch runs, in seconds.
const BATCH_SECONDS: f64 = 0.01;
/// A run at least this long, in seconds, is timed once: noise is small next to it.
const ONCE_SECONDS: f64 = 0.1;

/// The time of one call of `run`, in seconds, and the value it returned.
///
/// A first call sizes the batches; a short call is then repeated in `BATCHES` batches of about
/// `BATCH_SECONDS` each, and the fastest batch's average is the time, since anything slowing a
/// batch down comes from outside the aligner.
fn measure(mut run: impl FnMut() -> i64) -> (f64, i64) {
    let started = Instant::now();
    let mut value = run();
    let once = started.elapsed().as_secs_f64();
    if once >= ONCE_SECONDS {
        return (once, value);
    }
    let size = ((BATCH_SECONDS / once.max(1e-9)) as usize).max(1);
    let mut best = f64::INFINITY;
    for _ in 0..BATCHES {
        let started = Instant::now();
        for _ in 0..size {
            value = run();
        }
        best = best.min(started.elapsed().as_secs_f64() / size as f64);
    }
    (best, value)
}

/// A*PA2's evaluation datasets: `seq <tool> <budget seconds> <file>...`, every pair aligned once.
///
/// pa-bench's `.seq` files hold pairs as a `>` line and a `<` line. Each pair is aligned with its
/// traceback, once, as pa-bench times every aligner, until the budget is spent; one row per pair,
/// flushed as it finishes, gives its file, its time and its cost, so a run stopped mid-pair still
/// reports the pairs before it. A*PA takes the evaluation's
/// settings as `a*pa r=<r> prune=<start|both>`.
fn seq_mode(arguments: &[String]) {
    let tool = arguments[0].as_str();
    let budget: f64 = arguments[1].parse().unwrap();
    let started = Instant::now();
    while started.elapsed().as_millis() < 50 {
        std::hint::black_box(0);
    }
    let mut align: Box<dyn FnMut(&[u8], &[u8]) -> i64> = match tool {
        "a*pa2-full" | "a*pa2-simple" => {
            let params = if tool == "a*pa2-full" { AstarPa2Params::full() } else { AstarPa2Params::simple() };
            let mut aligner = params.make_aligner(true);
            Box::new(move |a, b| aligner.align_with_stats(a, b).0 as i64)
        }
        _ => {
            // `a*pa r=2 prune=start`, the evaluation's GCSH with diagonal transition and k = 15.
            let r: u8 = tool.split("r=").nth(1).unwrap().split(' ').next().unwrap().parse().unwrap();
            let prune = if tool.ends_with("prune=both") { Prune::Both } else { Prune::Start };
            Box::new(move |a, b| astarpa::astarpa_gcsh(a, b, r, 15, prune).0 as i64)
        }
    };
    let mut spent = 0.0;
    for path in &arguments[2..] {
        if spent >= budget {
            break;
        }
        let text = fs::read_to_string(path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        for pair in lines.chunks(2) {
            if pair.len() < 2 || spent >= budget {
                break;
            }
            let (a, b) = (pair[0][1..].as_bytes(), pair[1][1..].as_bytes());
            let started = Instant::now();
            let cost = align(a, b);
            let seconds = started.elapsed().as_secs_f64();
            spent += seconds;
            // Rust's stdout flushes by line, piped or not.
            println!("{tool}\t{path}\t{seconds}\t{cost}");
        }
    }
}

fn main() {
    let arguments: Vec<String> = env::args().skip(1).collect();
    if arguments.first().map(String::as_str) == Some("seq") {
        return seq_mode(&arguments[1..]);
    }
    if arguments.first().map(String::as_str) == Some("stats") {
        // `stats <file>`: A*PA2-full's own statistics for each pair, for comparing work, not time.
        let text = fs::read_to_string(&arguments[1]).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let mut aligner = AstarPa2Params::full().make_aligner(true);
        for pair in lines.chunks(2) {
            let (a, b) = (pair[0][1..].as_bytes(), pair[1][1..].as_bytes());
            let started = Instant::now();
            let (cost, _, stats) = aligner.align_with_stats(a, b);
            println!(
                "cost {cost} seconds {:.4} tries {} blocks {} incremental {} computed_lanes {} unique_lanes {} t_compute {:.4} t_j_range {:.4} t_fixed {:.4} t_pruning {:.4} t_contours {:.4} t_precomp {:.4}",
                started.elapsed().as_secs_f64(),
                stats.f_max_tries,
                stats.block_stats.num_blocks,
                stats.block_stats.num_incremental_blocks,
                stats.block_stats.computed_lanes,
                stats.block_stats.unique_lanes,
                stats.block_stats.t_compute.as_secs_f64(),
                stats.t_j_range.as_secs_f64(),
                stats.t_fixed_j_range.as_secs_f64(),
                stats.t_pruning.as_secs_f64(),
                stats.t_contours_update.as_secs_f64(),
                stats.t_precomp.as_secs_f64(),
            );
        }
        return;
    }
    let directory = env::args().nth(1).expect("usage: astarpa-runner <data directory>");
    // A short spin first, so the scheduler has moved this process onto a fast core.
    let started = Instant::now();
    while started.elapsed().as_millis() < 200 {
        std::hint::black_box(0);
    }
    for line in fs::read_to_string(format!("{directory}/dna_edit.tsv")).unwrap().lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let (name, a, b) = (fields[0], fields[1].as_bytes(), fields[2].as_bytes());

        for (tool, params) in [
            ("a*pa2-full", AstarPa2Params::full()),
            ("a*pa2-simple", AstarPa2Params::simple()),
            ("a*pa2-nw", AstarPa2Params::nw()),
        ] {
            for (task, trace) in [("score", false), ("alignment", true)] {
                let mut aligner = params.make_aligner(trace);
                let (seconds, cost) = measure(|| aligner.align_with_stats(a, b).0 as i64);
                // The same sum-and-weighted-sum checksum as every runner, for a one-pair workload.
                println!("{tool}\t{name}\t{task}\tcpu\t{seconds}\t{cost}:{cost}");
            }
        }

        // The original A*PA has no cost-only entry point, so it reports its alignment alone.
        let (seconds, cost) = measure(|| astarpa::astarpa(a, b).0 as i64);
        println!("a*pa\t{name}\talignment\tcpu\t{seconds}\t{cost}:{cost}");
    }
}
