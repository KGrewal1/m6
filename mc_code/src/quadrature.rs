use rand_distr::Uniform;
use rayon::prelude::*;
use std::collections::VecDeque;

use crate::monad_rng::{MonadicRng, UniformMonad};

pub fn trapezoidal<F: Fn(f64) -> f64 + Sync>(mut range: VecDeque<f64>, function: F) -> f64 {
    let n_segments = range.len() - 1;
    let f_0 = range.pop_front().map(|x| function(x) / 2.).unwrap_or(0.);

    let f_fin = range.pop_back().map(|x| function(x) / 2.).unwrap_or(0.);

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
    let rng = MonadicRng::new(seed);
    let range = UniformMonad::new(Uniform::from(range_bottom..range_top));
    let (sum, _rng) = (0..n_samples).fold((0., rng), |(sum, rng), _| {
        let (x, rng) = range.sample(rng);
        (sum + function(x), rng)
    });

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
    let rng = MonadicRng::new(seed);
    let range = UniformMonad::new(Uniform::from(range_bottom..range_top));

    let (init_pos, rng) = range.sample(rng);
    let init_wf = pdf(init_pos);
    let (sum, _final_pos, _final_wf, _rng) =
        (0..n_samples).fold((0_f64, init_pos, init_wf, rng), |(sum, pos, wf, rng), _| {
            let val = function(pos) / wf;
            let (new_pos, rng) = range.sample(rng);
            let new_wf = pdf(new_pos);
            let (rand, rng) = rng.gen_val::<f64>();
            if rand < new_wf / wf {
                (sum + val, new_pos, new_wf, rng)
            } else {
                (sum + val, pos, wf, rng)
            }
        });
    sum / (n_samples as f64)
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
    let rng = MonadicRng::new(seed);
    let (sum, _rng) = (0..n_samples).fold((0_f64, rng), |(sum, rng), _| {
        let (x, rng) = rng.gen_val::<f64>();
        let x = inverse_cdf(x);
        let val = function(x) / pdf(x);
        (sum + val, rng)
    });

    sum / (n_samples as f64)
}
