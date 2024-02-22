use crate::{Atom, MolSystem, MolSystems};
use rayon::prelude::*;
impl MolSystem {
    /// Calculate the virial pressure of a system given the derivative of the
    /// pairwise potential
    pub fn virial_pressure<F: Fn(&MolSystem, &Atom, &Atom) -> f64 + Sync>(
        &self,
        rdpotential: &F,
    ) -> f64 {
        self.pv(rdpotential) / (self.bb_volume())
    }

    // returns PV
    fn pv<F: Fn(&MolSystem, &Atom, &Atom) -> f64 + Sync>(&self, rdpotential: &F) -> f64 {
        self.atoms
            .iter()
            .enumerate()
            .map(|(i, atom_1)| -> f64 {
                tke(&atom_1)
                    + self
                        .atoms
                        .iter()
                        .skip(i + 1)
                        .map(|atom_2| -> f64 { rdpotential(self, atom_1, atom_2) })
                        .sum::<f64>()
            })
            .sum::<f64>()
            / (3.)
    }

    fn total_energy_int<F: Fn(&MolSystem, &Atom, &Atom) -> f64 + Sync>(&self, potential: F) -> f64 {
        self.atoms
            .par_iter()
            .enumerate()
            .map(|(i, atom_1)| {
                tke(&atom_1) / 2.
                    + self
                        .atoms
                        .iter()
                        .skip(i + 1)
                        .map(|atom_2| -> f64 { potential(self, atom_1, atom_2) })
                        .sum::<f64>()
            })
            .sum()
    }
}

impl MolSystems {
    /// Calculate the virial pressure of all systems given the derivative of the
    /// pairwise potential: see [`MolSystem::virial_pressure`]
    pub fn virial_pressure<F: Fn(&MolSystem, &Atom, &Atom) -> f64 + Sync + Send>(
        &self,
        rdpotential: F,
    ) -> (Vec<u32>, Vec<f64>) {
        (
            self.iter().map(|system| system.time).collect(),
            self.par_iter()
                .map(|system| system.virial_pressure(&rdpotential))
                .collect(),
        )
    }

    /// Calculate the enthalpy of all systems given the pairwise potential and the derivative of the
    /// pairwise potential: see [`MolSystem::virial_pressure`] and [`MolSystem::total_energy_int`]
    pub fn enthalpy<
        F: Fn(&MolSystem, &Atom, &Atom) -> f64 + Sync + Send,
        G: Fn(&MolSystem, &Atom, &Atom) -> f64 + Sync + Send,
    >(
        &self,
        potential: G,
        rdpotential: F,
    ) -> (Vec<u32>, Vec<f64>) {
        (
            self.iter().map(|system| system.time).collect(),
            self.par_iter()
                .map(|system| system.pv(&rdpotential) + system.total_energy_int(&potential))
                .collect(),
        )
    }
}

fn tke(atom: &Atom) -> f64 {
    atom.vx.powi(2) + atom.vy.powi(2) + atom.vz.powi(2)
}
