//! The A*PA side of the comparison: global edit distance on the DNA workloads only, since A*PA
//! supports nothing else.
//!
//! Each aligner runs twice per pair: with traceback, which returns a CIGAR, and without, which
//! returns only the cost, so both of dinara-align's tasks have a like-for-like row.
//!
//! Each measurement is taken warm, in this process, the way `ours.mojo` times dinara-align's
//! bit-parallel rows, so microsecond workloads compare like for like (see `measure`).

use astarpa2::AstarPa2Params;
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

fn main() {
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
