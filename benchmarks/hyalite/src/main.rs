//! The hyalite side of the comparison: global (`Mode::Nw`) affine-gap alignment of the DNA batches
//! and pairs `run.py` wrote, through `align_pair` (score) and `align` (alignment), on one CPU thread;
//! and with `local <file>`, `local_bench.py`'s local (`Mode::Sw`) or overlap (`Mode::Ov`) alignments.
//!
//! Every workload reads `dna_scoring.txt`, which dinara-align writes from its own default, so both
//! tools score with the same numbers. The edit-distance pairs are left to the bit-parallel aligners:
//! hyalite would answer them with a full quadratic sweep. Prints the same tab-separated rows as every
//! other runner, each measured warm and in-process as they all measure (see `measure`).

use hyalite::{Mode, Scoring, SearchType, align, align_pair};
use std::{env, fs, time::Instant};

/// Traceback working memory hyalite may use before it switches to its checkpoint path.
const BUDGET: usize = 1 << 30;
/// Batches whose average is reported: the fastest of them.
const BATCHES: usize = 20;
/// About how long one batch runs, in seconds.
const BATCH_SECONDS: f64 = 0.01;
/// A run at least this long, in seconds, is timed once: noise is small next to it.
const ONCE_SECONDS: f64 = 0.1;

struct Pairs {
    names: Vec<String>,
    firsts: Vec<Vec<u8>>,
    seconds: Vec<Vec<u8>>,
}

fn read(path: &str, alphabet: &str) -> Pairs {
    let mut table = [255u8; 256];
    for (index, letter) in alphabet.bytes().enumerate() {
        table[letter as usize] = index as u8;
    }
    let encode = |text: &str| text.bytes().map(|b| table[b as usize]).collect::<Vec<u8>>();
    let mut pairs = Pairs { names: Vec::new(), firsts: Vec::new(), seconds: Vec::new() };
    for line in fs::read_to_string(path).unwrap().lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        pairs.names.push(fields[0].to_string());
        pairs.firsts.push(encode(fields[1]));
        pairs.seconds.push(encode(fields[2]));
    }
    pairs
}

/// The sum and a position-weighted sum, the same checksum every runner prints.
fn checksum(values: &[i64]) -> String {
    let total: i64 = values.iter().sum();
    let weighted: i64 = values.iter().enumerate().map(|(i, v)| (i as i64 + 1) * v).sum();
    format!("{total}:{weighted}")
}

/// The time of one call of `run`, in seconds, and the value it returned.
///
/// A first call sizes the batches; a short call is then repeated in `BATCHES` batches of about
/// `BATCH_SECONDS` each, and the fastest batch's average is the time. A long call is timed once.
fn measure<T>(mut run: impl FnMut() -> T) -> (f64, T) {
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

/// Times the chosen pairs as one workload: every score, then every alignment.
fn time(workload: &str, pairs: &Pairs, chosen: &[usize], scoring: &Scoring) {
    let (seconds, scored) = measure(|| {
        chosen
            .iter()
            .map(|&i| align_pair(&pairs.firsts[i], &pairs.seconds[i], scoring, Mode::Nw, SearchType::Score).unwrap().score as i64)
            .collect::<Vec<i64>>()
    });
    println!("hyalite\t{workload}\tscore\tcpu\t{seconds}\t{}", checksum(&scored));
    let (seconds, aligned) = measure(|| {
        chosen
            .iter()
            .map(|&i| align(&pairs.firsts[i], &pairs.seconds[i], scoring, Mode::Nw, BUDGET).unwrap().score as i64)
            .collect::<Vec<i64>>()
    });
    println!("hyalite\t{workload}\talignment\tcpu\t{seconds}\t{}", checksum(&aligned));
}

fn batch(pairs: &Pairs, scoring: &Scoring) {
    let every: Vec<usize> = (0..pairs.names.len()).collect();
    time(&pairs.names[0], pairs, &every, scoring);
}

fn each(pairs: &Pairs, scoring: &Scoring) {
    for index in 0..pairs.names.len() {
        time(&pairs.names[index], pairs, &[index], scoring);
    }
}

/// `local_bench.py`'s workload, local, overlap or infix (`Mode::Hw`) alignment with traceback at a match 2, a mismatch -4
/// and a gap of `k` letters `6 + 2k`, which hyalite charges as `8 + 2 (k - 1)`; the faster of two
/// passes, its mean per pair, and the scores' checksum.
fn local(path: &str) {
    let pairs = read(path, "ACGT");
    let mut matrix = vec![-4; 16];
    for letter in 0..4 {
        matrix[letter * 4 + letter] = 2;
    }
    let scoring = Scoring::new(4, matrix, 8, 2).unwrap();
    let overlap = pairs.names[0].contains("overlap");
    let infix = pairs.names[0].contains("infix");
    let mode = if overlap {
        Mode::Ov
    } else if infix {
        Mode::Hw
    } else {
        Mode::Sw
    };
    let mut best = f64::INFINITY;
    let mut scores = Vec::new();
    for _ in 0..2 {
        let started = Instant::now();
        scores = (0..pairs.names.len())
            .map(|i| align(&pairs.seconds[i], &pairs.firsts[i], &scoring, mode, BUDGET).unwrap().score as i64)
            .collect();
        best = best.min(started.elapsed().as_secs_f64());
    }
    let task = if overlap {
        "overlap"
    } else if infix {
        "infix"
    } else {
        "local"
    };
    println!("hyalite\t{}\t{task}\t{}\t{}", pairs.names[0], best / pairs.names.len() as f64, checksum(&scores));
}

fn main() {
    if env::args().nth(1).as_deref() == Some("local") {
        local(&env::args().nth(2).expect("usage: hyalite-runner local <workload file>"));
        return;
    }
    let directory = env::args().nth(1).expect("usage: hyalite-runner <data directory>");
    let written = fs::read_to_string(format!("{directory}/dna_scoring.txt")).unwrap();
    let mut lines = written.lines();
    let alphabet = lines.next().unwrap().to_string();
    // dinara-align writes its penalties as negative scores; hyalite takes their magnitudes, and
    // charges `open` for a gap's first base and `ext` for each after it, as dinara-align does.
    let open: i32 = -lines.next().unwrap().parse::<i32>().unwrap();
    let ext: i32 = -lines.next().unwrap().parse::<i32>().unwrap();
    let matrix: Vec<i32> = lines.map(|l| l.parse().unwrap()).collect();
    let dna = Scoring::new(alphabet.len(), matrix, open, ext).unwrap();

    // A short spin first, so the scheduler has moved this process onto a fast core.
    let started = Instant::now();
    while started.elapsed().as_millis() < 200 {
        std::hint::black_box(0);
    }
    batch(&read(&format!("{directory}/dna_reads.tsv"), &alphabet), &dna);
    batch(&read(&format!("{directory}/dna_kilobase.tsv"), &alphabet), &dna);
    each(&read(&format!("{directory}/dna_affine.tsv"), &alphabet), &dna);
}
