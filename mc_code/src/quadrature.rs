use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Uniform};
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

    (f_0 + f_fin + sum_centre) / (n_segments as f64)
}

pub fn uniform_sample<F: Fn(f64) -> f64 + Sync>(
    range_bottom: f64,
    range_top: f64,
    n_samples: usize,
    function: F,
    seed: u64,
) -> f64 {
    let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
    let range = Uniform::from(range_bottom..range_top);
    let sum = (0..n_samples)
        .map(|_| {
            let x = range.sample(&mut rng);
            function(x)
        })
        .sum::<f64>();

    sum / (n_samples as f64)
}

/// importance sampling, stochastically sampling the weight function,
/// and then averaging the function values weighted by the inverse of the weight function
pub fn importance_sample<F: Fn(f64) -> f64 + Sync, G: Fn(f64) -> f64 + Sync>(
    range_bottom: f64,
    range_top: f64,
    n_samples: usize,
    function: F,
    pdf: G,
    seed: u64,
) -> f64 {
    let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
    let range = Uniform::from(range_bottom..range_top);
    let mut sample_vals = Vec::with_capacity(n_samples);
    while sample_vals.len() != n_samples {
        let x = range.sample(&mut rng);
        let y: f64 = rng.gen();
        let wf = pdf(x);
        if y < wf {
            sample_vals.push(function(x) / wf);
        }
    }

    sample_vals.par_iter().sum::<f64>() / (n_samples as f64)
}

/// importance sampling transforming a uniform distribution over [0,1] to the weight function
/// by using its inverse cdf, and then averaging the function values weighted by the inverse of the weight function
pub fn importance_sample_alt<
    F: Fn(f64) -> f64 + Sync,
    G: Fn(f64) -> f64 + Sync,
    H: Fn(f64) -> f64 + Sync,
>(
    n_samples: usize,
    function: F,
    pdf: G,
    inverse_cdf: H,
    seed: u64,
) -> f64 {
    let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
    let sum = (0..n_samples)
        .map(|_| {
            let x = rng.gen::<f64>();
            let x = inverse_cdf(x);
            function(x) / pdf(x)
        })
        .sum::<f64>();

    sum / (n_samples as f64)
}
