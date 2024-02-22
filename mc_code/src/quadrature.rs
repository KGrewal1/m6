use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
use std::collections::VecDeque;

pub fn trapezoidal<F: Fn(f64) -> f64 + Sync>(mut range: VecDeque<f64>, function: F) -> f64 {
    let n_segments = range.len() - 1;
    let n_0 = range.pop_front().unwrap_or(0.);
    let f_0 = function(n_0) / 2.;
    let n_fin = range.pop_back().unwrap_or(0.);
    let f_fin = function(n_fin) / 2.;
    let sum_centre = range.par_iter().map(|x| function(*x)).sum::<f64>();

    return (f_0 + f_fin + sum_centre) / (n_segments as f64);
}

pub fn uniform_sample<F: Fn(f64) -> f64 + Sync>(
    range_bottom: f64,
    range_top: f64,
    n_samples: usize,
    function: F,
    seed: u64,
) -> f64 {
    let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
    let sum = (0..n_samples)
        .map(|_| {
            let x = rng.gen_range(range_bottom..range_top);
            function(x)
        })
        .sum::<f64>();
    sum / (n_samples as f64)
}

pub fn importance_sample<F: Fn(f64) -> f64 + Sync, G: Fn(f64) -> f64 + Sync>(
    range_bottom: f64,
    range_top: f64,
    mut n_samples: usize,
    function: F,
    weight_function: G,
    seed: u64,
) -> f64 {
    let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
    let mut sample_vals = Vec::with_capacity(n_samples);
    while sample_vals.len() != n_samples {
        let x = rng.gen_range(range_bottom..range_top);
        let y = rng.gen_range(0.0..1.0);
        let wf = weight_function(x);
        if y < wf {
            sample_vals.push(function(x) / weight_function(x));
        }
    }
    sample_vals.par_iter().sum::<f64>() / (n_samples as f64)
}
