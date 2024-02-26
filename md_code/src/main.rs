#![allow(non_upper_case_globals)]
#[cfg(feature = "png")]
use plotly::ImageFormat;
use plotly::{
    common::Title,
    layout::{Axis, GridPattern, LayoutGrid, RowOrder},
    Layout, Plot, Scatter,
};
use rayon::prelude::*;
use std::{fs, path::Path, time::Instant};
use xyz_tools::{Atom, MolSystem, MolSystems};

const k_b: f64 = 1.380_649; // Joules per kelvin

// conversion factor from kcal to joules
const kcal_to_j: f64 = 4184.;

// for diffusion question
const sigma_diff: f64 = 3.405; // Angstroms

#[allow(clippy::too_many_lines)]
fn main() {
    let m6_dir = Path::new("M6_files");
    let m6_files = m6_dir.join("enviar");
    let now_global = Instant::now();

    let plots_dir = Path::new("plots");
    if !plots_dir.exists() {
        fs::create_dir(plots_dir).expect("Unable to create directory");
    }
    // --------------------------------
    // 2.1 RDF of an ideal gas
    // --------------------------------

    {
        let now = Instant::now();
        let systems: MolSystems = fs::read_to_string(m6_files.join("ideal.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();

        // we know all atoms are the same so default to not selecting which atom for slightly faster path
        let (r, g_r) = systems.rdf(None, 1000, 20.);

        let mut plot = Plot::new();
        let trace = Scatter::new(r, g_r).name("Autocorrelation Function");
        let layout = Layout::new()
            .title(Title::new("Autocorrelation of ideal.xyz"))
            .x_axis(Axis::new().title(Title::new("Distance / Angstrom")))
            .y_axis(Axis::new().title(Title::new("g(r)")));
        plot.set_layout(layout);
        plot.add_trace(trace);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex2_1.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex2_1.html"));
        println!("Ex 2.1 took {} ms\n", now.elapsed().as_millis());
    }

    // --------------------------------
    // 2.2 RDF of bulk water
    // --------------------------------

    {
        let now = Instant::now();
        let systems: MolSystems = fs::read_to_string(m6_files.join("water.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();

        let (r, g_r) = systems.rdf(Some((1, 1)), 1000, 15.);

        let mut plot = Plot::new();
        let trace = Scatter::new(r, g_r).name("Autocorrelation Function");
        let layout = Layout::new()
            .title(Title::new("O-O Autocorrelation of water.xyz"))
            .x_axis(Axis::new().title(Title::new("Distance (Angstroms)")))
            .y_axis(Axis::new().title(Title::new("g(r)")));
        plot.set_layout(layout);
        plot.add_trace(trace);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex2_2_a.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex2_2_a.html"));

        let (r, g_r) = systems.rdf(Some((1, 2)), 1000, 15.);

        let mut plot = Plot::new();
        let trace = Scatter::new(r, g_r).name("Autocorrelation Function");
        let layout = Layout::new()
            .title(Title::new("O-H Autocorrelation of water.xyz"))
            .x_axis(Axis::new().title(Title::new("Distance (Angstroms)")))
            .y_axis(Axis::new().title(Title::new("g(r)")).range(vec![0., 2.00]));
        plot.set_layout(layout);
        plot.add_trace(trace);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex2_2_b.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex2_2_b.html"));
        println!("Ex 2.2 took {} ms\n", now.elapsed().as_millis());
    }

    // --------------------------------
    // 2.3 Characterising different crystal and liquid phases through their RDFs
    // --------------------------------

    {
        let now = Instant::now();
        let systems: MolSystems = fs::read_to_string(m6_files.join("set1.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();

        let (r, g_r) = systems.rdf(None, 1000, 20.);

        // let mut plot = Plot::new();
        let trace_1 = Scatter::new(r, g_r)
            .name("Autocorrelation Function")
            .name("set1.xyz");

        let systems: MolSystems = fs::read_to_string(m6_files.join("set2.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();

        let (r, g_r) = systems.rdf(None, 1000, 20.);

        // let mut plot = Plot::new();
        let trace_2 = Scatter::new(r, g_r)
            .name("Autocorrelation Function")
            .x_axis("x2")
            .y_axis("y2")
            .name("set2.xyz");

        let systems: MolSystems = fs::read_to_string(m6_files.join("set3.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();

        let (r, g_r) = systems.rdf(None, 1000, 20.);

        // let mut plot = Plot::new();
        let trace_3 = Scatter::new(r, g_r)
            .name("Autocorrelation Function")
            .x_axis("x3")
            .y_axis("y3")
            .name("set3.xyz");

        let mut plot = Plot::new();
        plot.add_trace(trace_1);
        plot.add_trace(trace_2);
        plot.add_trace(trace_3);
        let layout = Layout::new()
            .grid(
                LayoutGrid::new()
                    .rows(3)
                    .columns(1)
                    .pattern(GridPattern::Independent)
                    .row_order(RowOrder::TopToBottom),
            )
            .height(2000)
            .width(2500)
            .title(Title::new("Autocorrelation of Sets"))
            .x_axis3(Axis::new().title(Title::new("Distance / Angstroms")))
            .y_axis(Axis::new().title(Title::new("g(r)")))
            .y_axis2(Axis::new().title(Title::new("g(r)")))
            .y_axis3(Axis::new().title(Title::new("g(r)")));
        plot.set_layout(layout);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex2_3_comb.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex2_3_comb.html"));
        println!("Ex 2.3 took {} ms\n", now.elapsed().as_millis());
    }

    // --------------------------------
    // 3.1 Potential Energy
    // --------------------------------

    {
        let now = Instant::now();
        let system: MolSystem = fs::read_to_string(m6_files.join("conf.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();
        let energy = system.potential_enegry(lj_potential);
        println!("LJ energy (J): {:.2}", energy);

        let energy = system.potential_enegry(phs_potential);
        println!("Pseudo Hard Sphere energy (J): {:.2}", energy);

        let energy = system.potential_enegry(yd_potential);
        println!("Yukawa Debye energy (J): {:.2}", energy);
        println!("Ex 3.1 took {} ms\n", now.elapsed().as_millis());
    }

    // --------------------------------
    // 3.2 Kinetic Energy
    // --------------------------------

    {
        let now = Instant::now();
        let mut system: MolSystem = fs::read_to_string(m6_files.join("conf.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();
        let temp = 179.81;
        let m = 6.63e-26;
        system.random_v(k_b, temp, m, 3_624_360);
        let ke_1 = system.mean_ke(m);

        system.random_v(k_b, temp, m, 42);
        let ke_2 = system.mean_ke(m);

        system.random_v(k_b, temp, m, 0);
        let ke_3 = system.mean_ke(m);

        let exact_ke = 0.5 * (3. * k_b * temp);

        println!("Mean KE 1 (J): {:.2}", ke_1);
        println!("Mean KE 2 (J): {:.2}", ke_2);
        println!("Mean KE 3 (J): {:.2}", ke_3);
        println!("Exact KE (J): {:.2}", exact_ke);
        println!("Ex 3.2 took {} ms\n", now.elapsed().as_millis());
    } //1.126e-17J expected KE

    // --------------------------------
    // 4.1 Pressure through the virial expression
    // --------------------------------
    // angstroms per emtosecond

    {
        let now = Instant::now();
        let systems: MolSystems = fs::read_to_string(m6_files.join("pres.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();
        let (times, pressure) = systems.virial_pressure(rdlj_potential);

        let mut plot = Plot::new();
        let trace = Scatter::new(times, pressure).name("LJ Pressure");
        plot.add_trace(trace);
        let layout = Layout::new()
            .title(Title::new("LJ Pressure"))
            .x_axis(Axis::new().title(Title::new("Time")))
            .y_axis(Axis::new().title(Title::new("Pressure / Joules per cubic Angstrom")));
        plot.set_layout(layout);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex4_1_a.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex4_1_a.html"));

        let (times, pressure) = systems.virial_pressure(rdphs_potential);

        let mut plot = Plot::new();
        let trace = Scatter::new(times, pressure).name("PHS Pressure");
        plot.add_trace(trace);

        let layout = Layout::new()
            .title(Title::new("PHS Pressure"))
            .x_axis(Axis::new().title(Title::new("Time")))
            .y_axis(Axis::new().title(Title::new("Pressure / Joules per cubic Angstrom")));
        plot.set_layout(layout);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex4_1_b.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex4_1_b.html"));

        println!("Ex 4.1 took {} ms\n", now.elapsed().as_millis());
    }

    // --------------------------------
    // 4.2 Time evolution of enthalpy
    // --------------------------------

    {
        let now = Instant::now();
        let systems: MolSystems = fs::read_to_string(m6_files.join("pres.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();
        let (times, pressure) = systems.enthalpy(lj_potential, rdlj_potential);

        let mut plot = Plot::new();
        let trace = Scatter::new(times, pressure).name("LJ Enthalpy");
        plot.add_trace(trace);
        let layout = Layout::new()
            .title(Title::new("LJ Enthalpy"))
            .x_axis(Axis::new().title(Title::new("Time")))
            .y_axis(Axis::new().title(Title::new("H / Joules")));
        plot.set_layout(layout);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex4_2_a.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex4_2_a.html"));

        let (times, pressure) = systems.enthalpy(phs_potential, rdphs_potential);

        let mut plot = Plot::new();
        let trace = Scatter::new(times, pressure).name("PHS Enthalpy");
        plot.add_trace(trace);
        let layout = Layout::new()
            .title(Title::new("PHS Enthalpy"))
            .x_axis(Axis::new().title(Title::new("Time")))
            .y_axis(Axis::new().title(Title::new("H / Joules")));
        plot.set_layout(layout);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex4_2_b.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex4_2_b.html"));

        println!("Ex 4.2 took {} ms\n", now.elapsed().as_millis());
    }

    // --------------------------------
    // 5 Dynamic properties: Diffusion coefficient
    // --------------------------------

    {
        let now = Instant::now();
        // for diffusion question: fp maths not callable from const context
        let tau_diff: f64 = (6.63e-26 * sigma_diff.powi(2) / (0.24 * kcal_to_j)).sqrt() * 1e12; // femtoseconds
        let mut systems: MolSystems = fs::read_to_string(m6_files.join("diffusionA.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();

        let times = (0..systems.len()).map(|i| 2 * i).collect::<Vec<_>>();
        let msds = systems.diffusion();
        let tt = times.par_iter().map(|t| (t * t) as f64).sum::<f64>();
        let mt = times
            .par_iter()
            .zip(msds.par_iter())
            .map(|(t, msd)| *t as f64 * msd)
            .sum::<f64>();
        println!("slope = {:.5}", mt / tt);
        println!("D = {:.5}", (mt / sigma_diff) / (tt * 2. / tau_diff));
        let mut plot = Plot::new();
        let trace = Scatter::new(times, msds).name("Diffuson of A");
        let layout = Layout::new()
            .title(Title::new("MSD of diffusionA.xyz"))
            .x_axis(Axis::new().title(Title::new("Time / fs")))
            .y_axis(Axis::new().title(Title::new("MSD / A")));
        plot.set_layout(layout);
        plot.add_trace(trace);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex5_a.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex5_a.html"));
        println!("Ex 5a took {} ms\n", now.elapsed().as_millis());
    }

    {
        let now = Instant::now();
        // for diffusion question: fp maths not callable from const context
        let tau_diff: f64 = (6.63e-26 * sigma_diff.powi(2) / (0.24 * kcal_to_j)).sqrt() * 1e12; // femtoseconds
        let mut systems: MolSystems = fs::read_to_string(m6_files.join("diffusionB.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();

        let times = (0..systems.len()).map(|i| 2 * i).collect::<Vec<_>>();
        let msds = systems.diffusion();
        let tt = times.par_iter().map(|t| (t * t) as f64).sum::<f64>();
        let mt = times
            .par_iter()
            .zip(msds.par_iter())
            .map(|(t, msd)| *t as f64 * msd)
            .sum::<f64>();
        println!("slope = {:.5}", mt / tt);
        println!("D = {:.5}", (mt / sigma_diff) / (tt * 2. / tau_diff));
        println!("This however may not be the most valid straight line fit");
        let mut plot = Plot::new();
        let trace = Scatter::new(times, msds).name("Diffuson of B");
        let layout = Layout::new()
            .title(Title::new("MSD of diffusionB.xyz"))
            .x_axis(Axis::new().title(Title::new("Time / fs")))
            .y_axis(Axis::new().title(Title::new("MSD / A")));
        plot.set_layout(layout);
        plot.add_trace(trace);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex5_b.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex5_b.html"));
        println!("Ex 5b took {} ms\n", now.elapsed().as_millis());
    }

    // --------------------------------
    // 6.1 Number of contacts through a local order parameter
    // --------------------------------

    {
        let now = Instant::now();
        let t264k: MolSystems = fs::read_to_string(m6_files.join("264K.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();
        // let n = systems.len();/
        let nc264 = t264k.mean_contact(4.1);
        println!("Mean contact at 264K: {:.2}", nc264);

        let t276k: MolSystems = fs::read_to_string(m6_files.join("276K.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();
        // let n = systems.len();/
        let nc276 = t276k.mean_contact(4.1);
        println!("Mean contact at 278K: {:.2}", nc276);

        let t288k: MolSystems = fs::read_to_string(m6_files.join("288K.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();
        // let n = systems.len();/
        let nc288 = t288k.mean_contact(4.1);
        println!("Mean contact at 288K: {:.2}", nc288);

        let t300k: MolSystems = fs::read_to_string(m6_files.join("300K.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();
        // let n = systems.len();/
        let nc300 = t300k.mean_contact(4.1);
        println!("Mean contact at 300K: {:.2}", nc300);

        let temps = vec![264, 276, 288, 300];
        let ncs = vec![nc264, nc276, nc288, nc300];
        let mut plot = Plot::new();
        let trace = Scatter::new(temps, ncs).name("Diffuson of B");
        let layout = Layout::new()
            .title(Title::new("Mean contact number of polymer"))
            .x_axis(Axis::new().title(Title::new("Temp / K")))
            .y_axis(Axis::new().title(Title::new("Mean Contact Number")));
        plot.set_layout(layout);
        plot.add_trace(trace);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex6_1.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex6_1.html"));
        println!("Ex 6.1 took {} ms\n", now.elapsed().as_millis());
    }

    // --------------------------------
    // 6.2 Computing the phase diagram
    // --------------------------------

    {
        // The filters for the two means are based solely on graphs
        let mut p_t: Vec<(f64, f64)> = Vec::with_capacity(8);
        let now = Instant::now();
        let t264k: MolSystems = fs::read_to_string(m6_files.join("264K.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();
        let (boxes, count) = t264k.z_dist(100, 3.4);
        let highs: Vec<f64> = count.iter().filter(|x| **x > 0.45).copied().collect();
        let lows: Vec<f64> = count.iter().filter(|x| **x < 0.05).copied().collect();

        let high_avg = highs.iter().sum::<f64>() / highs.len() as f64;
        let low_avg = lows.iter().sum::<f64>() / lows.len() as f64;
        println!("Mean high density 266K: {}", high_avg);
        println!("Mean low density 266K: {}", low_avg);
        p_t.push((264., high_avg));
        p_t.push((264., low_avg));

        let mut plot = Plot::new();
        let trace = Scatter::new(boxes, count).name("Dispersion at 264 K");
        let layout = Layout::new()
            .title(Title::new("Dispersion at 264 K"))
            .x_axis(Axis::new().title(Title::new("Box")))
            .y_axis(Axis::new().title(Title::new("Number Density")));
        plot.set_layout(layout);
        plot.add_trace(trace);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex6_2_a.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex6_2_a.html"));

        let t276k: MolSystems = fs::read_to_string(m6_files.join("276K.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();
        let (boxes, count) = t276k.z_dist(100, 3.4);
        let highs: Vec<f64> = count.iter().filter(|x| **x > 0.45).copied().collect();
        let lows: Vec<f64> = count.iter().filter(|x| **x < 0.05).copied().collect();

        let high_avg = highs.iter().sum::<f64>() / highs.len() as f64;
        let low_avg = lows.iter().sum::<f64>() / lows.len() as f64;
        println!("Mean high density 276K: {}", high_avg);
        println!("Mean low density 276K: {}", low_avg);
        p_t.push((276., high_avg));
        p_t.push((276., low_avg));

        let mut plot = Plot::new();
        let trace = Scatter::new(boxes, count).name("Dispersion at 276 K");
        let layout = Layout::new()
            .title(Title::new("Dispersion at 276 K"))
            .x_axis(Axis::new().title(Title::new("Box")))
            .y_axis(Axis::new().title(Title::new("Number Density")));
        plot.set_layout(layout);
        plot.add_trace(trace);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex6_2_b.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex6_2_b.html"));

        let t288k: MolSystems = fs::read_to_string(m6_files.join("288K.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();
        let (boxes, count) = t288k.z_dist(100, 3.4);

        let highs: Vec<f64> = count.iter().filter(|x| **x > 0.35).copied().collect();
        let lows: Vec<f64> = count.iter().filter(|x| **x < 0.05).copied().collect();
        let high_avg = highs.iter().sum::<f64>() / highs.len() as f64;
        let low_avg = lows.iter().sum::<f64>() / lows.len() as f64;
        println!("Mean high density 288K: {}", high_avg);
        println!("Mean low density 288K: {}", low_avg);
        p_t.push((288., high_avg));
        p_t.push((288., low_avg));

        let mut plot = Plot::new();
        let trace = Scatter::new(boxes, count).name("Dispersion at 288 K");
        let layout = Layout::new()
            .title(Title::new("Dispersion at 288 K"))
            .x_axis(Axis::new().title(Title::new("Box")))
            .y_axis(Axis::new().title(Title::new("Number Density")));
        plot.set_layout(layout);
        plot.add_trace(trace);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex6_2_c.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex6_2_c.html"));

        let t300k: MolSystems = fs::read_to_string(m6_files.join("300K.xyz"))
            .expect("Unable to read the file")
            .parse()
            .unwrap();
        let (boxes, count) = t300k.z_dist(100, 3.4);

        let highs: Vec<f64> = count.iter().filter(|x| **x > 0.215).copied().collect();
        let lows: Vec<f64> = count.iter().filter(|x| **x < 0.18).copied().collect();

        let high_avg = highs.iter().sum::<f64>() / highs.len() as f64;
        let low_avg = lows.iter().sum::<f64>() / lows.len() as f64;
        println!("Mean high density 300K: {}", high_avg);
        println!("Mean low density 300K: {}", low_avg);
        p_t.push((300., high_avg));
        p_t.push((300., low_avg));

        let mut plot = Plot::new();
        let trace = Scatter::new(boxes, count).name("Dispersion at 300 K");
        let layout = Layout::new()
            .title(Title::new("Dispersion at 300 K"))
            .x_axis(Axis::new().title(Title::new("Box")))
            .y_axis(Axis::new().title(Title::new("Number Density")));
        plot.set_layout(layout);
        plot.add_trace(trace);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex6_2_d.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex6_2_d.html"));

        // println!("{:?}", p_t);
        p_t.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        // println!("{:?}", p_t);
        let mut pressures: Vec<f64> = Vec::with_capacity(8);
        let mut temperatures: Vec<f64> = Vec::with_capacity(8);
        for (t, p) in p_t {
            temperatures.push(t);
            pressures.push(p);
        }
        let mut plot = Plot::new();
        let trace = Scatter::new(pressures, temperatures).name("PT");
        let layout = Layout::new()
            .title(Title::new("PT phase diagram"))
            .x_axis(Axis::new().title(Title::new("Number Density")))
            .y_axis(Axis::new().title(Title::new("Temperature / K")));
        plot.set_layout(layout);
        plot.add_trace(trace);
        #[cfg(feature = "png")]
        plot.write_image(
            plots_dir.join("m6_ex6_2_e.png"),
            ImageFormat::PNG,
            800,
            600,
            1.0,
        );
        plot.write_html(plots_dir.join("m6_ex6_2_e.html"));
        println!("Ex 6.2 took {} ms\n", now.elapsed().as_millis());
    }

    println!("MD exercises took {} ms", now_global.elapsed().as_millis());
}

// Potential energy functions
fn lj_potential(system: &MolSystem, atom_1: &Atom, atom_2: &Atom) -> f64 {
    let r = system.distance(atom_1, atom_2);
    let sigma = 3.405;
    let epsilon = 119.87 * k_b; // in K
    if r < 3. * sigma {
        let s6 = (sigma / r).powi(6);
        4. * epsilon * (s6.powi(2) - s6)
    } else {
        0.
    }
}

fn rdlj_potential(system: &MolSystem, atom_1: &Atom, atom_2: &Atom) -> f64 {
    let r = system.distance(atom_1, atom_2);
    let sigma = 3.405;
    let epsilon = 119.87 * k_b; // in K
    if r < 3. * sigma {
        let s6r7 = sigma.powi(6) / r.powi(7);
        r * 4. * epsilon * (6. * s6r7 - 12. * r * s6r7.powi(2))
    } else {
        0.
    }
}

fn phs_potential(system: &MolSystem, atom_1: &Atom, atom_2: &Atom) -> f64 {
    let r = system.distance(atom_1, atom_2);
    let lambda_a = 49_f64;
    let lambda_r = 50_f64;
    let sigma = 3.405;
    // in K
    let epsilon_r = 119.87 * k_b;
    // no need for another cutoff as 0 before regardless
    if r < (lambda_r / lambda_a) * sigma {
        lambda_r
            * (lambda_a / lambda_r).powf(lambda_a)
            * epsilon_r
            * ((sigma / r).powf(lambda_r) - (sigma / r).powf(lambda_a))
            + epsilon_r
    } else {
        0.
    }
}

fn rdphs_potential(system: &MolSystem, atom_1: &Atom, atom_2: &Atom) -> f64 {
    let r = system.distance(atom_1, atom_2);
    let lambda_a = 49_f64;
    let lambda_r = 50_f64;
    let sigma = 3.405;
    // in K
    let epsilon_r = 119.87 * k_b;
    // no need for another cutoff as 0 before regardless
    if r < (lambda_r / lambda_a) * sigma {
        r * (lambda_r
            * (lambda_a / lambda_r).powf(lambda_a)
            * epsilon_r
            * (lambda_a * sigma.powf(lambda_a) / r.powf(lambda_a + 1.)
                - lambda_r * sigma.powf(lambda_r) / r.powf(lambda_r + 1.))
            + epsilon_r)
    } else {
        0.
    }
}

fn yd_potential(system: &MolSystem, atom_1: &Atom, atom_2: &Atom) -> f64 {
    let r = system.distance(atom_1, atom_2);
    let sigma = 3.05;
    // in K
    let epsilon = 119.87 * k_b;
    let kappa = 5_f64;
    match r {
        r if r < sigma => 0.,
        r if ((sigma <= r) & (r <= 2.5 * sigma)) => {
            let sign = if atom_1.at_type == atom_2.at_type {
                1.
            } else {
                -1.
            };
            sign * epsilon * sigma * (-kappa * (r - sigma)).exp() / r
        }
        _ => 0_f64,
    }
}
