use crate::{Atom, MolSystem, MolSystems, Periodicity, Scaling};
use rayon::prelude::*;
use std::{f64::consts::PI, vec};

impl MolSystem {
    /// radial distribuion function of a single system
    #[must_use]
    pub fn rdf(
        &self,
        atom_types: Option<(u8, u8)>,
        resolution: usize,
        cutoff: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        #[allow(clippy::cast_precision_loss)]
        let dr = cutoff / ((resolution - 1) as f64);
        #[allow(clippy::cast_precision_loss)]
        let radii: Vec<f64> = (0..resolution).map(|i| dr * (i as f64)).collect();
        #[allow(clippy::cast_precision_loss)]
        let volumes: Vec<f64> = (0..resolution)
            .map(|i| 4. / 3. * PI * ((dr * ((i + 1) as f64)).powi(3) - (dr * (i as f64)).powi(3)))
            .collect();

        let g_r = self.rdf_internal(atom_types, resolution, &volumes, dr);
        (radii, g_r)
    }

    fn bucket(&self, atom_1s: &[&Atom], atom_2s: &[&Atom], resolution: usize, dr: f64) -> Vec<f64> {
        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss
        )]
        atom_1s
            .iter()
            .enumerate()
            .flat_map(|(i, atom_1)| {
                atom_2s
                    .iter()
                    .skip(i + 1)
                    .map(|atom_2| -> f64 { self.distance(atom_1, atom_2) })
            })
            .fold(vec![0.0; resolution], |mut g_r: Vec<f64>, dist| {
                let i = (dist / dr) as usize; // should never lose sign as both strictly positive
                if i < resolution {
                    g_r[i] += 2.0;
                };
                g_r
            })
    }

    fn rdf_internal(
        &self,
        atom_types: Option<(u8, u8)>,
        resolution: usize,
        volumes: &[f64],
        dr: f64,
    ) -> Vec<f64> {
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

        if let Some((atom_type_1, atom_type_2)) = atom_types {
            let mut atom_1s: Vec<&Atom> = vec![];
            let mut atom_2s: Vec<&Atom> = vec![];
            for atom in &self.atoms {
                if atom.at_type == atom_type_1 {
                    atom_1s.push(atom);
                } else if atom.at_type == atom_type_2 {
                    atom_2s.push(atom);
                }
            }
            let n_atoms_1 = atom_1s.len();
            let n_atoms_2;

            let mut g_r = if atom_type_1 == atom_type_2 {
                n_atoms_2 = n_atoms_1;
                // println!("threads: {}", current_num_threads());
                self.bucket(&atom_1s, &atom_1s, resolution, dr)
            } else {
                n_atoms_2 = atom_2s.len();
                self.bucket(&atom_1s, &atom_2s, resolution, dr)
            };

            #[allow(clippy::cast_precision_loss)]
            g_r.iter_mut()
                .zip(volumes.iter())
                .for_each(|(g, v)| *g /= (n_atoms_1 * n_atoms_2) as f64 * v / self.bb_volume());

            g_r
        } else {
            #[allow(
                clippy::cast_precision_loss,
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss
            )]
            let mut g_r = self
                .atoms
                .iter()
                .enumerate()
                .flat_map(|(i, atom_1)| {
                    self.atoms
                        .iter()
                        .skip(i + 1)
                        .map(|atom_2| -> f64 { self.distance(atom_1, atom_2) })
                })
                .fold(vec![0.0; resolution], |mut g_r: Vec<f64>, dist| {
                    let i = (dist / dr) as usize;
                    if i < resolution {
                        g_r[i] += 2.0;
                    };
                    g_r
                });
            #[allow(clippy::cast_precision_loss)]
            g_r.iter_mut()
                .zip(volumes.iter())
                .for_each(|(g, v)| *g /= (self.n_atoms.pow(2)) as f64 * v / self.bb_volume());

            g_r
        }
    }
}

impl MolSystems {
    /// avergaes the RDF across systems
    /// makes asssumption that they're all the same system at different times
    #[must_use]
    pub fn rdf(
        &self,
        atom_types: Option<(u8, u8)>,
        resolution: usize,
        cutoff: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        #[allow(clippy::cast_precision_loss)]
        let n_systems = self.0.len() as f64;
        // let mut sum_g_r = vec![0.0; resolution];

        // let cutoff = 0.5;
        #[allow(clippy::cast_precision_loss)]
        let dr = cutoff / ((resolution) as f64); //- 1
        #[allow(clippy::cast_precision_loss)]
        let radii: Vec<f64> = (0..resolution).map(|i| dr * (i as f64)).collect();
        #[allow(clippy::cast_precision_loss)]
        let volumes: Vec<f64> = (0..resolution)
            .map(|i| 4. / 3. * PI * ((dr * ((i + 1) as f64)).powi(3) - (dr * (i as f64)).powi(3)))
            .collect();

        let sum_g_r = self
            .0
            .par_iter()
            .map(|system| system.rdf_internal(atom_types, resolution, &volumes, dr))
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

        // self.0.iter().for_each(|system| {
        //     let g_r = system.rdf_internal(atom_types, resolution, &volumes, dr);
        //     sum_g_r
        //         .iter_mut()
        //         .zip(g_r.iter())
        //         .for_each(|(sum, g)| *sum += g);
        // });

        (radii, sum_g_r.iter().map(|g| g / n_systems).collect())
    }
}
