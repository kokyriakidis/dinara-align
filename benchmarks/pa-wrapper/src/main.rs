//! The aligners A*PA2's evaluation compared against, on the edit-distance workloads: Edlib, WFA2-lib's
//! BiWFA and WFA, and WFA-adaptive, each through pa-bench's `pa-wrapper` with the parameters that
//! evaluation used. Block Aligner is left out, as it scores affine costs only; KSW2, as its kernels
//! are SSE only; TripleAccel, as its quadratic time makes a 100 kbp pair take minutes.
//!
//! Each aligner runs twice per pair, with traceback and without, as the A*PA runner does, and each
//! measurement is taken warm and in-process the same way every runner takes its own (see `measure`).
//! WFA-adaptive drops lagging diagonals and is not exact; `run.py` leaves its answers out of the
//! agreement check and marks a time whose answer was not optimal.

use pa_types::CostModel;
use pa_wrapper::{
    wrappers::{edlib::EdlibParams, wfa::WfaParams},
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

fn main() {
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
        (
            "wfa-adaptive",
            AlignerParams::Wfa(WfaParams {
                memory_model: MemoryModel::MemoryUltraLow,
                heuristic: Heuristic::WFadaptive(10, 50, 10),
            }),
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
