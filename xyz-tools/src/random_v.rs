use crate::MolSystem;
use rand::prelude::*;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
impl MolSystem {
    pub fn random_v(&mut self, k: f64, temp: f64, mass: f64, seed: u64) {
        let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
        let normal = rand_distr::Normal::new(0.0, (k * temp / mass).sqrt()).unwrap();
        self.atoms.iter_mut().for_each(|atom| {
            atom.vx = normal.sample(&mut rng);
            atom.vy = normal.sample(&mut rng);
            atom.vz = normal.sample(&mut rng);
        });
    }

    pub fn mean_ke(&self, mass: f64) -> f64 {
        self.atoms
            .par_iter()
            .map(|atom| 0.5 * mass * (atom.vx.powi(2) + atom.vy.powi(2) + atom.vz.powi(2)))
            .sum::<f64>()
            / (self.n_atoms as f64)
    }
}
