use std::{collections::VecDeque, fs, path::Path, time::Instant};

use quadrature::trapezoidal;
use rayon::iter::{IntoParallelRefIterator, ParallelBridge, ParallelIterator};

use crate::quadrature::{importance_sample, importance_sample_alt, uniform_sample};

mod photon_gas;
mod quadrature;

fn main() {
    let plots_dir = Path::new("plots");
    if !plots_dir.exists() {
        fs::create_dir(plots_dir).expect("Unable to create directory");
    }
    //-----------------
    // MC q 1
    //----------------
    {
        let now = Instant::now();
        photon_gas::photon_gas(plots_dir);
        println!("MC Ex 1 took {} ms\n", now.elapsed().as_millis());
    }

    //-----------------
    // MC q 2
    //----------------
    {
        let now = Instant::now();
        println!("Exact integral 1.0");

        {
            let now = Instant::now();
            let range: VecDeque<f64> = (0..1_000_000).map(|i| 1e-6 * f64::from(i)).collect();

            let integral = trapezoidal(range, three_x_sq);

            println!("Trapezoidal integral: {}", integral);
            println!("Trapezoidal took {} ms\n", now.elapsed().as_millis());
        }

        {
            let now = Instant::now();
            let integrals: Vec<f64> = (1..=1000)
                .par_bridge()
                .map(|i| uniform_sample(0., 1., 1_000_000, three_x_sq, i))
                .collect();

            let expected_integral = integrals.par_iter().sum::<f64>() / 1000.;

            let variance = integrals
                .par_iter()
                .map(|x| (x - expected_integral).powi(2))
                .sum::<f64>()
                / 1000.;

            println!(
                "Uniform Sampling integral: {}, variance {}",
                expected_integral, variance
            );
            println!("Uniform took {} ms\n", now.elapsed().as_millis());
        }

        {
            let now = Instant::now();
            let integrals: Vec<f64> = (1..=1000)
                .par_bridge()
                .map(|i| importance_sample(0., 1., 1_000_000, three_x_sq, two_x, i))
                .collect();

            let expected_integral = integrals.par_iter().sum::<f64>() / 1000.;

            let variance = integrals
                .par_iter()
                .map(|x| (x - expected_integral).powi(2))
                .sum::<f64>()
                / 1000.;

            println!(
                "2x Sampling integral: {}, variance {}",
                expected_integral, variance
            );
            println!("2x took {} ms\n", now.elapsed().as_millis());
        }

        {
            let now = Instant::now();
            let integrals: Vec<f64> = (1..=1000)
                .par_bridge()
                .map(|i| importance_sample(0., 1., 1_000_000, three_x_sq, four_x_cubed, i))
                .collect();

            let expected_integral = integrals.par_iter().sum::<f64>() / 1000.;

            let variance = integrals
                .par_iter()
                .map(|x| (x - expected_integral).powi(2))
                .sum::<f64>()
                / 1000.;

            println!(
                "4x^3 Sampling integral: {}, variance {}",
                expected_integral, variance
            );
            println!("4x^3 took {} ms\n", now.elapsed().as_millis());
        }

        {
            let now = Instant::now();
            let integrals: Vec<f64> = (1..=1000)
                .par_bridge()
                .map(|i| importance_sample(0., 1., 1_000_000, three_x_sq, three_x_sq, i))
                .collect();

            let expected_integral = integrals.par_iter().sum::<f64>() / 1000.;

            let variance = integrals
                .par_iter()
                .map(|x| (x - expected_integral).powi(2))
                .sum::<f64>()
                / 1000.;

            println!(
                "3x^2 Sampling integral: {}, variance {}",
                expected_integral, variance
            );
            println!("3x^2 took {} ms\n", now.elapsed().as_millis());
        }

        // using inversion method

        {
            let now = Instant::now();
            let integrals: Vec<f64> = (1..=1000)
                .par_bridge()
                .map(|i| {
                    importance_sample_alt(0., 1., 1_000_000, three_x_sq, two_x, two_x_inv_cdf, i)
                })
                .collect();

            let expected_integral = integrals.par_iter().sum::<f64>() / 1000.;

            let variance = integrals
                .par_iter()
                .map(|x| (x - expected_integral).powi(2))
                .sum::<f64>()
                / 1000.;

            println!(
                "2x Sampling integral: {}, variance {}",
                expected_integral, variance
            );
            println!("2x took {} ms\n", now.elapsed().as_millis());
        }

        {
            let now = Instant::now();
            let integrals: Vec<f64> = (1..=1000)
                .par_bridge()
                .map(|i| {
                    importance_sample_alt(
                        0.,
                        1.,
                        1_000_000,
                        three_x_sq,
                        four_x_cubed,
                        four_x_cubed_inv_cdf,
                        i,
                    )
                })
                .collect();

            let expected_integral = integrals.par_iter().sum::<f64>() / 1000.;

            let variance = integrals
                .par_iter()
                .map(|x| (x - expected_integral).powi(2))
                .sum::<f64>()
                / 1000.;

            println!(
                "4x^3 Sampling integral: {}, variance {}",
                expected_integral, variance
            );
            println!("4x^3 took {} ms\n", now.elapsed().as_millis());
        }

        {
            let now = Instant::now();
            let integrals: Vec<f64> = (1..=1000)
                .par_bridge()
                .map(|i| {
                    importance_sample_alt(
                        0.,
                        1.,
                        1_000_000,
                        three_x_sq,
                        three_x_sq,
                        three_x_sq_inv_cdf,
                        i,
                    )
                })
                .collect();

            let expected_integral = integrals.par_iter().sum::<f64>() / 1000.;

            let variance = integrals
                .par_iter()
                .map(|x| (x - expected_integral).powi(2))
                .sum::<f64>()
                / 1000.;

            println!(
                "3x^2 Sampling integral: {}, variance {}",
                expected_integral, variance
            );
            println!("3x^2 took {} ms\n", now.elapsed().as_millis());
        }

        println!("MC Ex 2 took {} ms\n", now.elapsed().as_millis());
    }
}

// 3x2 pdf
fn three_x_sq(x: f64) -> f64 {
    3. * x.powi(2)
}

/// 3x2 inverse cdf

fn three_x_sq_inv_cdf(x: f64) -> f64 {
    x.cbrt()
}

/// 2x pdf
fn two_x(x: f64) -> f64 {
    2. * x
}

/// two x inverse cdf
fn two_x_inv_cdf(x: f64) -> f64 {
    x.sqrt()
}

/// 4x^3 pdf
fn four_x_cubed(x: f64) -> f64 {
    4. * x.powi(3)
}

/// 4x^3 inverse cdf
fn four_x_cubed_inv_cdf(x: f64) -> f64 {
    x.powf(1. / 4.)
}
