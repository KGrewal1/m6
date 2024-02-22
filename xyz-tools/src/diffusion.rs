//! tools to calculate the rate of diffusion of
//! a system evolving in time
use crate::{MolSystem, MolSystems};
use rayon::prelude::*;

impl MolSystems {
    #[must_use]
    /// method to calculate the diffusion of a system
    /// over time by calculating the mean square displacement
    /// from the original configuration
    pub fn diffusion(&mut self) -> Vec<f64> {
        self.0.par_iter_mut().for_each(|system| {
            system.centre_atoms();
            system.sort_atoms();
        });
        let initial_state = &self.0[0];
        self.0
            .par_iter()
            .map(|system| mean_sq_disp(initial_state, system))
            .collect()
    }
}

fn mean_sq_disp(system1: &MolSystem, system2: &MolSystem) -> f64 {
    system1
        .atoms
        .iter()
        .zip(system2.atoms.iter())
        .map(|(a, b)| system1.square_distance(a, b))
        .sum::<f64>()
        / (system1.n_atoms as f64)
}
