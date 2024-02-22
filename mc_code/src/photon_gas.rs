#[cfg(feature = "png")]
use plotly::ImageFormat;
use plotly::{common::Title, layout::Axis, Layout, Plot, Scatter};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
use std::path::Path;

const MC_STEPS: usize = 25_000_000;
#[allow(non_upper_case_globals)]
const e_j: f64 = 1.;

pub fn photon_gas(path: &Path) {
    let vals = (1..=20)
        .par_bridge()
        .map(|i| 0.1 * f64::from(i))
        .map(|beta| (beta, mc_calc::<MC_STEPS>(beta), 1. / (beta * e_j).exp_m1()))
        .collect::<Vec<(f64, f64, f64)>>();
    let (beta, mc_j, exact_j) = transpose(vals);
    let delta = mc_j
        .iter()
        .zip(exact_j.iter())
        .map(|(a, b)| (a - b))
        .collect::<Vec<f64>>();

    let mut plot = Plot::new();
    let trace_1 = Scatter::new(beta.clone(), mc_j).name("MC value");
    let trace_2 = Scatter::new(beta.clone(), exact_j).name("Exact value");
    let layout = Layout::new()
        .title(Title::new("Photon gas"))
        .x_axis(Axis::new().title(Title::new("Beta")))
        .y_axis(Axis::new().title(Title::new("<nj>")));
    plot.set_layout(layout);
    plot.add_trace(trace_1);
    plot.add_trace(trace_2);
    #[cfg(feature = "png")]
    plot.write_image(path.join("m6_mc_1.png"), ImageFormat::PNG, 800, 600, 1.0);
    plot.write_html(path.join("m6_mc_1.html"));

    let mut plot = Plot::new();
    let trace_1 = Scatter::new(beta.clone(), delta).name("MC delta");
    let layout = Layout::new()
        .title(Title::new("Photon gas"))
        .x_axis(Axis::new().title(Title::new("Beta")))
        .y_axis(Axis::new().title(Title::new("<nj> difference")));
    plot.set_layout(layout);
    plot.add_trace(trace_1);
    #[cfg(feature = "png")]
    plot.write_image(
        path.join("m6_mc_1_delta.png"),
        ImageFormat::PNG,
        800,
        600,
        1.0,
    );
    plot.write_html(path.join("m6_mc_1_delta.html"));
}

// calculate the occupation number for the photon gas in a specific state
fn mc_calc<const NSTEPS: usize>(beta: f64) -> f64 {
    let init_state: u16 = 10;
    let mut rng = Xoshiro256StarStar::seed_from_u64(42);
    let (sum, _final_state) = (0..NSTEPS).fold((0_f64, init_state), move |(sum, state), _| {
        let direction: bool = rng.gen();
        let new_state = if direction {
            state.saturating_add(1)
        } else {
            state.saturating_sub(1) // saturating sub prevents overflow (unsigned cannot be negative)
        };
        let u_init = e_j * f64::from(state);
        let u_final = e_j * f64::from(new_state);
        let p = (-beta * (u_final - u_init)).exp();
        let final_state = if rng.gen::<f64>() < p {
            new_state
        } else {
            state
        };

        (sum + (f64::from(final_state) / NSTEPS as f64), final_state)
    });
    sum
}

// transpose from vec of (f64, f64, f64) to three vecs
fn transpose(mut v: Vec<(f64, f64, f64)>) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let n_elements = v.len();
    let mut a_vec: Vec<f64> = Vec::with_capacity(n_elements);
    let mut b_vec: Vec<f64> = Vec::with_capacity(n_elements);
    let mut c_vec: Vec<f64> = Vec::with_capacity(n_elements);
    for (a, b, c) in v {
        a_vec.push(a);
        b_vec.push(b);
        c_vec.push(c);
    }
    (a_vec, b_vec, c_vec)
}
