use crate::{Atom, MolSystem, MolSystems, Periodicity, Scaling};
use rayon::prelude::*;

impl MolSystem {
    fn n_molecules(&self) -> usize {
        // works on the assumption that counting starts at 1 and finishes at n molecules with no gaps
        self.atoms.par_iter().map(|atom| atom.mol).max().unwrap() as usize
    }

    /// get the mean number of contacts for atoms in the system, excluding bonded atoms
    pub fn mean_contact(&self, cutoff: f64) -> f64 {
        self.atoms
            .par_iter()
            .enumerate()
            .flat_map_iter(|(i, atom_1)| {
                self.atoms.iter().skip(i + 1).map(|atom_2| -> f64 {
                    if (self.distance(atom_1, atom_2) < cutoff) && !is_bonded(atom_1, atom_2) {
                        1.
                    } else {
                        0.
                    }
                })
            })
            .sum::<f64>()
            / (self.n_molecules() as f64)
    }

    /// get the z axis distribution of the system
    pub fn z_dist(&self, resolution: usize, lunit: f64) -> (Vec<f64>, Vec<f64>) {
        let boxes = (0..resolution).map(|i| i as f64).collect();
        let dr = 1.0 / ((resolution - 1) as f64);
        (boxes, self.z_dist_int(resolution, dr, lunit.powi(3)))
    }

    fn z_dist_int(&self, resolution: usize, dr: f64, vunit: f64) -> Vec<f64> {
        assert_eq!(
            self.periodicity,
            Periodicity::Periodic,
            "Non-periodic systems are not supported"
        );

        assert_eq!(
            self.scaling,
            Scaling::Scaled,
            "Unscaled systems are not supported"
        );

        let volume = self.bb_volume() * dr / vunit;

        let g_r = self.atoms.iter().map(|atom_1| atom_1.z).fold(
            vec![0.0; resolution],
            |mut g_r: Vec<f64>, z| {
                let i = (z / dr) as usize;
                if i < resolution {
                    g_r[i] += 1.0;
                };
                g_r
            },
        );

        // #[allow(clippy::cast_precision_loss)]
        // g_r.iter_mut()
        //     .zip(volumes.iter())
        //     .for_each(|(g, v)| *g /= (self.n_atoms.pow(2)) as f64 * v / self.bb_volume());
        g_r.into_iter().map(|n_atoms| n_atoms / volume).collect()
    }
}

impl MolSystems {
    /// get the mean contact number for all atoms in the systems (see [`MolSystem::mean_contact`] for more details)
    pub fn mean_contact(&self, cutoff: f64) -> f64 {
        self.0
            .par_iter()
            .map(|system| system.mean_contact(cutoff))
            .sum::<f64>()
            / (self.0.len() as f64)
    }

    /// get the z axis distribution for all systems (see [`MolSystem::z_dist`] for more details)
    pub fn z_dist(&self, resolution: usize, lunit: f64) -> (Vec<f64>, Vec<f64>) {
        let boxes = (0..resolution).map(|i| i as f64).collect();
        let dr = 1.0 / ((resolution) as f64);
        let g_r = self
            .0
            .par_iter()
            .map(|system| system.z_dist_int(resolution, dr, lunit.powi(3)))
            .reduce(
                || vec![0.0; resolution],
                |mut g_r_1: Vec<f64>, g_r_2: Vec<f64>| {
                    g_r_1
                        .iter_mut()
                        .zip(g_r_2.iter())
                        .for_each(|(g_1, g_2)| *g_1 += g_2);
                    g_r_1
                },
            );

        (
            boxes,
            g_r.into_iter()
                .map(|n_atoms| n_atoms / (self.0.len() as f64))
                .collect(),
        )
    }
}

fn is_bonded(atom_1: &Atom, atom_2: &Atom) -> bool {
    // atoms are bonded if their id is adjacent and they are in the same molecule
    (atom_1.id.abs_diff(atom_2.id) == 1) && (atom_1.mol == atom_2.mol)
}
