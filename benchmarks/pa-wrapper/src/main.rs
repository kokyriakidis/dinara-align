//! The exact aligners A*PA2's evaluation compared against, on the edit-distance workloads: Edlib, and
//! WFA2-lib's BiWFA and WFA, each through pa-bench's `pa-wrapper` with the parameters that evaluation
//! used. Every column of the comparison computes the optimum, so WFA-adaptive and Block Aligner, which
//! may not, are left out; so is KSW2, as its kernels are SSE only, and TripleAccel, as its quadratic
//! time makes a 100 kbp pair take minutes.
//!
//! Each aligner runs twice per pair, with traceback and without, as the A*PA runner does, and each
//! measurement is taken warm and in-process the same way every runner takes its own (see `measure`).

use pa_types::CostModel;
use pa_wrapper::{
    wrappers::{
        block_aligner::{BlockAlignerParams, BlockAlignerSize},
        edlib::EdlibParams,
        wfa::WfaParams,
    },
    AlignerParams,
};
use rust_wfa2::aligner::{Heuristic, MemoryModel};
use std::{env, fs, time::Instant};

/// Batches whose average is reported: the fastest of them.
const BATCHES: usize = 20;
/// About how long one batch runs, in seconds.
const BATCH_SECONDS: f64 = 0.01;
/// A run at least this long, in seconds, is timed once: noise is small next to it.
const ONCE_SECONDS: f64 = 0.1;

/// The time of one call of `run`, in seconds, and the value it returned, as every runner measures.
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
/// reports the pairs before it. `biwfa` is the evaluation's WFA2-lib, its lowest-memory mode, and
/// `wfa` its keep-every-front mode, whose memory grows with the square of the distance; the harness
/// caps it (see `pa_bench.py`).
///
/// `wfa-adaptive` and `block-aligner` are the evaluation's two approximate aligners, with its
/// parameters: WFA-adaptive's defaults (10, 50, 10), and Block Aligner's blocks from 0.1 to 1% of the
/// input, at a gap-opening cost of one as it takes only affine costs. Either may return a worse
/// alignment than the optimum at its costs, so the harness holds each against an exact one at the
/// same costs, as the evaluation does: the exact distance, and for Block Aligner `biwfa-affine`,
/// BiWFA at its affine costs.
fn seq_mode(arguments: &[String]) {
    let tool = arguments[0].as_str();
    let budget: f64 = arguments[1].parse().unwrap();
    let started = Instant::now();
    while started.elapsed().as_millis() < 50 {
        std::hint::black_box(0);
    }
    let params = match tool {
        "edlib" => AlignerParams::Edlib(EdlibParams),
        "biwfa" => AlignerParams::Wfa(WfaParams { memory_model: MemoryModel::MemoryUltraLow, heuristic: Heuristic::None }),
        "wfa" => AlignerParams::Wfa(WfaParams { memory_model: MemoryModel::MemoryHigh, heuristic: Heuristic::None }),
        "wfa-adaptive" => AlignerParams::Wfa(WfaParams {
            memory_model: MemoryModel::MemoryUltraLow,
            heuristic: Heuristic::WFadaptive(10, 50, 10),
        }),
        "block-aligner" => {
            AlignerParams::BlockAligner(BlockAlignerParams { size: BlockAlignerSize::Percent(0.001, 0.01) })
        }
        "biwfa-affine" => {
            AlignerParams::Wfa(WfaParams { memory_model: MemoryModel::MemoryUltraLow, heuristic: Heuristic::None })
        }
        other => panic!("unknown tool {other}"),
    };
    // Block Aligner takes only affine costs, a gap's opening one more; its exact reference shares them.
    let costs = match tool {
        "block-aligner" | "biwfa-affine" => CostModel::affine(1, 1, 1),
        _ => CostModel::unit(),
    };
    let mut spent = 0.0;
    for path in &arguments[2..] {
        if spent >= budget {
            break;
        }
        let text = fs::read_to_string(path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let longest = lines.iter().map(|line| line.len()).max().unwrap_or(1);
        let (mut aligner, _exact) = params.build_aligner(costs, true, longest);
        for pair in lines.chunks(2) {
            if pair.len() < 2 || spent >= budget {
                break;
            }
            let (a, b) = (pair[0][1..].as_bytes(), pair[1][1..].as_bytes());
            let before = peak_resident();
            let started = Instant::now();
            let cost = aligner.align(a, b).0;
            let seconds = started.elapsed().as_secs_f64();
            let growth = peak_resident() - before;
            spent += seconds;
            // Rust's stdout flushes by line, piped or not.
            println!("{tool}\t{path}\t{seconds}\t{cost}\t{growth}");
        }
    }
}


/// The process's peak resident memory so far, in bytes: `getrusage` counts kilobytes on Linux and bytes
/// on macOS. Read before and after each alignment, its growth is A*PA2's memory measure.
fn peak_resident() -> i64 {
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage) };
    let peak = usage.ru_maxrss as i64;
    if cfg!(target_os = "macos") { peak } else { peak * 1024 }
}

fn main() {
    let arguments: Vec<String> = env::args().skip(1).collect();
    if arguments.first().map(String::as_str) == Some("seq") {
        return seq_mode(&arguments[1..]);
    }
    let directory = env::args().nth(1).expect("usage: pa-wrapper-runner <data directory>");
    // A short spin first, so the scheduler has moved this process onto a fast core.
    let started = Instant::now();
    while started.elapsed().as_millis() < 200 {
        std::hint::black_box(0);
    }
    let aligners: Vec<(&str, AlignerParams)> = vec![
        ("edlib", AlignerParams::Edlib(EdlibParams)),
        (
            "biwfa",
            AlignerParams::Wfa(WfaParams { memory_model: MemoryModel::MemoryUltraLow, heuristic: Heuristic::None }),
        ),
        (
            "wfa",
            AlignerParams::Wfa(WfaParams { memory_model: MemoryModel::MemoryHigh, heuristic: Heuristic::None }),
        ),
    ];
    for line in fs::read_to_string(format!("{directory}/dna_edit.tsv")).unwrap().lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let (name, a, b) = (fields[0], fields[1].as_bytes(), fields[2].as_bytes());
        let longest = a.len().max(b.len());
        for (tool, params) in &aligners {
            for (task, trace) in [("score", false), ("alignment", true)] {
                // WFA2-lib's score-only mode keeps no memory model but its lowest.
                let params = match params {
                    AlignerParams::Wfa(wfa) if !trace => {
                        AlignerParams::Wfa(WfaParams { memory_model: MemoryModel::MemoryUltraLow, ..*wfa })
                    }
                    other => other.clone(),
                };
                let (mut aligner, _exact) = params.build_aligner(CostModel::unit(), trace, longest);
                let (seconds, cost) = measure(|| aligner.align(a, b).0 as i64);
                println!("{tool}\t{name}\t{task}\tcpu\t{seconds}\t{cost}:{cost}");
            }
        }
    }
}
