mod deserializer;
pub mod diffusion;
pub mod polymer;
pub mod potential_energy;
pub mod pressure;
pub mod radial_distribution;
pub mod random_v;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Periodicity {
    Periodic,
    NonPeriodic,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Scaling {
    Scaled,
    Unscaled,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Atom {
    pub id: u32,
    pub mol: u32,
    pub at_type: u8,
    pub charge: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MolSystem {
    pub time: u32,
    pub n_atoms: usize,
    pub periodicity: Periodicity,
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
    pub z0: f64,
    pub z1: f64,
    pub scaling: Scaling,
    pub atoms: Vec<Atom>,
}

impl MolSystem {
    /// Returns the volume of the bounding box in cubic angstroms
    #[must_use]
    pub fn bb_volume(&self) -> f64 {
        (self.x1 - self.x0) * (self.y1 - self.y0) * (self.z1 - self.z0)
    }

    #[must_use]
    pub fn distance(&self, atom_a: &Atom, atom_b: &Atom) -> f64 {
        self.square_distance(atom_a, atom_b).sqrt()
    }

    #[must_use]
    fn square_distance(&self, atom_a: &Atom, atom_b: &Atom) -> f64 {
        // TODO: Adjust by scale factor
        let dx = (atom_a.x - atom_b.x).abs();
        let dy = (atom_a.y - atom_b.y).abs();
        let dz = (atom_a.z - atom_b.z).abs();
        // println!("dx: {}, dy: {}, dz: {}", dx, dy, dz);
        if self.scaling == Scaling::Scaled {
            let dist_x = self.x1 - self.x0;
            let dist_y = self.y1 - self.y0;
            let dist_z = self.z1 - self.z0;

            let (dx, dy, dz) = {
                (
                    dist_x * dx.min((1. - dx).abs()),
                    dist_y * dy.min((1. - dy).abs()),
                    dist_z * dz.min((1. - dz).abs()),
                )
            };
            dx * dx + dy * dy + dz * dz
        } else {
            dx * dx + dy * dy + dz * dz
        }
    }

    fn sort_atoms(&mut self) {
        self.atoms.sort_by(|a, b| a.id.cmp(&b.id));
    }

    fn com(&self) -> (f64, f64, f64) {
        let (x, y, z) = self.atoms.iter().fold((0., 0., 0.), |acc, atom| {
            (acc.0 + atom.x, acc.1 + atom.y, acc.2 + atom.z)
        });
        (
            x / self.n_atoms as f64,
            y / self.n_atoms as f64,
            z / self.n_atoms as f64,
        )
    }

    fn centre_atoms(&mut self) {
        let (x, y, z) = self.com();
        self.atoms.iter_mut().for_each(|atom| {
            atom.x -= x;
            atom.y -= y;
            atom.z -= z;
        });
    }
}
pub struct MolSystems(Vec<MolSystem>);

impl std::ops::Deref for MolSystems {
    type Target = Vec<MolSystem>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn com_test() {
        // test on LJ3
        let input = "ITEM: NUMBER OF ATOMS
        4 #number of atoms
        ITEM: BOX BOUNDS pp pp pp
        0e+00 1e+00 #xmin xmax
        0e+00 1e+00 #ymin ymax
        0e+00 1e+00 #zmin zmax
        ITEM: ATOMS id mol type q xu yu zu #q=charge, coordenates are scaled by the box side in each dimension (i.e. xmax-xmin)
        4 4 2 0 0 0 0
        2 2 2 0 0 2 4
        5 5 2 0 2 0 0
        8 8 2 0 2 2 6";

        let mut res: MolSystem = input.parse().unwrap();
        let (x, y, z) = res.com();
        assert_approx_eq!(x, 1.);
        assert_approx_eq!(y, 1.);
        assert_approx_eq!(z, 2.5);
        res.sort_atoms();
        res.centre_atoms();

        let atom2 = Atom {
            id: 2,
            mol: 2,
            at_type: 2,
            charge: 0.,
            x: -1.,
            y: 1.,
            z: 1.5,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };

        let atom4 = Atom {
            id: 4,
            mol: 4,
            at_type: 2,
            charge: 0.,
            x: -1.,
            y: -1.,
            z: -2.5,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };

        let atom5 = Atom {
            id: 5,
            mol: 5,
            at_type: 2,
            charge: 0.,
            x: 1.,
            y: -1.,
            z: -2.5,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };

        let atom8 = Atom {
            id: 8,
            mol: 8,
            at_type: 2,
            charge: 0.,
            x: 1.,
            y: 1.,
            z: 3.5,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        assert_eq!(res.atoms[0], atom2);
        assert_eq!(res.atoms[1], atom4);
        assert_eq!(res.atoms[2], atom5);
        assert_eq!(res.atoms[3], atom8);
    }
}
