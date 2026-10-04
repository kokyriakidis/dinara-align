//! The hyalite side of the comparison: global (`Mode::Nw`) alignment of every workload `run.py` wrote,
//! pair by pair through `align_pair` (score) and `align` (alignment), on one CPU thread.
//!
//! The affine workloads read `dna_scoring.txt`, which dinara-align writes from its own default, so
//! both tools score with the same numbers. Prints the same tab-separated rows as every other runner.

use hyalite::{Mode, Scoring, SearchType, align, align_pair};
use std::{env, fs, time::Instant};

/// Traceback working memory hyalite may use before it switches to its checkpoint path.
const BUDGET: usize = 1 << 30;

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

fn emit(workload: &str, task: &str, started: Instant, values: &[i64]) {
    let seconds = started.elapsed().as_secs_f64();
    println!("hyalite\t{workload}\t{task}\tcpu\t{seconds}\t{}", checksum(values));
}

/// Times the chosen pairs as one workload: every score first, then every alignment.
fn time(workload: &str, pairs: &Pairs, chosen: &[usize], scoring: &Scoring, sign: i64) {
    let started = Instant::now();
    let scored: Vec<i64> = chosen
        .iter()
        .map(|&i| {
            let hit = align_pair(&pairs.firsts[i], &pairs.seconds[i], scoring, Mode::Nw, SearchType::Score);
            sign * hit.unwrap().score as i64
        })
        .collect();
    emit(workload, "score", started, &scored);
    let started = Instant::now();
    let aligned: Vec<i64> = chosen
        .iter()
        .map(|&i| sign * align(&pairs.firsts[i], &pairs.seconds[i], scoring, Mode::Nw, BUDGET).unwrap().score as i64)
        .collect();
    emit(workload, "alignment", started, &aligned);
}

fn batch(pairs: &Pairs, scoring: &Scoring) {
    let every: Vec<usize> = (0..pairs.names.len()).collect();
    time(&pairs.names[0], pairs, &every, scoring, 1);
}

fn each(pairs: &Pairs, scoring: &Scoring, sign: i64) {
    for index in 0..pairs.names.len() {
        time(&pairs.names[index], pairs, &[index], scoring, sign);
    }
}

fn main() {
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

    batch(&read(&format!("{directory}/dna_reads.tsv"), &alphabet), &dna);
    batch(&read(&format!("{directory}/dna_kilobase.tsv"), &alphabet), &dna);
    each(&read(&format!("{directory}/dna_affine.tsv"), &alphabet), &dna, 1);

    let mut unit = vec![-1; 16];
    for i in 0..4 {
        unit[i * 4 + i] = 0;
    }
    let unit = Scoring::new(4, unit, 1, 1).unwrap();
    each(&read(&format!("{directory}/dna_edit.tsv"), "ACGT"), &unit, -1);
}
