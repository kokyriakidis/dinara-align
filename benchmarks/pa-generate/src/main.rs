//! Writes A*PA2's synthetic datasets as pa-bench generates them, under the names pa-bench gives them.
//!
//!     pa-generate-runner <directory> <seed> <total size> <error rate>... -- <length>...
//!
//! Each `(error rate, length)` is one file of uniform-error pairs whose lengths add up to the total
//! size, as `GeneratedDataset::to_generator` sets the generator up in pa-bench.

use pa_generate::{DatasetGenerator, ErrorModel, SeqPairGenerator};
use std::{env, path::Path};

fn main() {
    let arguments: Vec<String> = env::args().skip(1).collect();
    let directory = Path::new(&arguments[0]);
    let seed: u64 = arguments[1].parse().unwrap();
    let total_size: usize = arguments[2].parse().unwrap();
    let split = arguments.iter().position(|a| a == "--").expect("missing -- before the lengths");
    let rates: Vec<f32> = arguments[3..split].iter().map(|a| a.parse().unwrap()).collect();
    let lengths: Vec<usize> = arguments[split + 1..].iter().map(|a| a.parse().unwrap()).collect();
    std::fs::create_dir_all(directory).unwrap();
    for &error_rate in &rates {
        for &length in &lengths {
            let name = format!("{:?}-t{}-n{}-e{}.seq", ErrorModel::Uniform, total_size, length, error_rate);
            let path = directory.join(&name);
            if path.exists() && path.metadata().unwrap().len() > 0 {
                continue;
            }
            DatasetGenerator {
                settings: SeqPairGenerator { length, error_rate, error_model: ErrorModel::Uniform, pattern_length: None },
                seed: Some(seed),
                cnt: None,
                size: Some(total_size),
            }
            .generate_file(&path);
            println!("{name}");
        }
    }
}
