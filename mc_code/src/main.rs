use std::{collections::VecDeque, time::Instant};

use quadrature::trapezoidal;

use crate::quadrature::{importance_sample, uniform_sample};

mod photon_gas;
mod quadrature;

fn main() {
    //-----------------
    // MC q 1
    //----------------
    {
        let now = Instant::now();
        photon_gas::photon_gas();
        println!("MC Ex 1 took {} ms\n", now.elapsed().as_millis());
    }

    {
        let now = Instant::now();
        println!("Exact integral 1.0");

        let range: VecDeque<f64> = (0..1_000_000).map(|i| 1e-6 * i as f64).collect();

        let integral = trapezoidal(range, three_x_sq);

        println!("Trapezoidal integral: {}", integral);

        let integral = uniform_sample(0., 1., 1_000_000, three_x_sq);

        println!("Uniform Sampling integral: {}", integral);

        let integral = importance_sample(0., 1., 1_000_000, three_x_sq, two_x);

        println!("2x Sampling integral: {}", integral);

        let integral = importance_sample(0., 1., 1_000_000, three_x_sq, four_x_cubed);

        println!("4x^3 Sampling integral: {}", integral);

        println!("MC Ex 2 took {} ms\n", now.elapsed().as_millis());
    }
}

fn three_x_sq(x: f64) -> f64 {
    3. * x.powi(2)
}
fn two_x(x: f64) -> f64 {
    2. * x
}

fn four_x_cubed(x: f64) -> f64 {
    4. * x.powi(3)
}
