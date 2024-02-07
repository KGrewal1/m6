use crate::{Atom, MolSystem};
use rayon::prelude::*;
impl MolSystem {
    pub fn potential_enegry<F: Fn(&MolSystem, &Atom, &Atom) -> f64 + Sync>(
        &self,
        potential: F,
    ) -> f64 {
        self.atoms
            .par_iter()
            .enumerate()
            .flat_map_iter(|(i, atom_1)| {
                self.atoms
                    .iter()
                    .skip(i + 1)
                    .map(|atom_2| -> f64 { potential(self, atom_1, atom_2) })
            })
            .sum()
    }

    // pub fn random_vel(&mut self, temperature: f64) {
    //     let mut rng = rand::thread_rng();
    //     let boltzmann = 8.617333262145e-5; // eV/K
    //     let mass = 1.0; // assume mass is 1
    //     let sigma = (boltzmann * temperature / mass).sqrt();
    //     self.atoms.iter_mut().for_each(|atom| {
    //         atom.vx = rng.;
    //         atom.vy = rand::distributions::Normal::new(0.0, sigma).sample(&mut rng);
    //         atom.vz = rand::distributions::Normal::new(0.0, sigma).sample(&mut rng);
    //     });
    // }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    fn lj_potential(system: &MolSystem, atom_1: &Atom, atom_2: &Atom) -> f64 {
        let r = system.distance(atom_1, atom_2);
        // println!("r: {}", r);
        4. * (1. / r.powi(12) - 1. / r.powi(6))
    }

    #[test]
    fn lj_clusters_test() {
        // test on LJ3
        let input = "ITEM: NUMBER OF ATOMS
        3 #number of atoms
        ITEM: BOX BOUNDS pp pp pp
        0e+00 1e+00 #xmin xmax
        0e+00 1e+00 #ymin ymax
        0e+00 1e+00 #zmin zmax
        ITEM: ATOMS id mol type q xu yu zu #q=charge, coordenates are scaled by the box side in each dimension (i.e. xmax-xmin)
        1 1 2 0  0.4391356726        0.1106588251       -0.4635601962
        2 2 2 0 -0.5185079933        0.3850176090        0.0537084789
        3 3 2 0  0.0793723207       -0.4956764341        0.4098517173";

        let res: MolSystem = input.parse().unwrap();
        let energy = res.potential_enegry(lj_potential);
        assert_approx_eq!(energy, -3.0, 1e-6);

        // test on LJ38
        let input = "ITEM: NUMBER OF ATOMS
        38 #number of atoms
        ITEM: BOX BOUNDS pp pp pp
        0e+00 1e+00 #xmin xmax
        0e+00 1e+00 #ymin ymax
        0e+00 1e+00 #zmin zmax
        ITEM: ATOMS id mol type q xu yu zu #q=charge, coordenates are scaled by the box side in each dimension (i.e. xmax-xmin)
        1  1  2 0  0.1947679907        0.3306365642        1.7069272101
        2  2  2 0  1.1592174250       -1.1514615100       -0.6254746298
        3  3  2 0  1.4851406793       -0.0676273830        0.9223060046
        4  4  2 0 -0.1498046416        1.4425168343       -0.9785553065
        5  5  2 0  1.4277261305        0.3530265376       -0.9475378022
        6  6  2 0 -0.6881246261       -1.5737014419       -0.3328844168
        7  7  2 0 -1.4277352637       -0.3530034531        0.9475270683
        8  8  2 0  0.6881257085        1.5736904826        0.3329032458
        9  9  2 0 -1.1592204530        1.1514535263        0.6254777879
        10 10 2 0  0.1498035273       -1.4424985165        0.9785685322
        11 11 2 0 -1.4851196066        0.0676193562       -0.9223231092
        12 12 2 0 -0.7057028384        0.6207073550       -1.4756523155
        13 13 2 0 -0.8745359533        0.4648140463        1.4422103492
        14 14 2 0 -0.9742077067       -0.8837261792       -1.1536019836
        15 15 2 0 -0.1947765396       -0.3306358487       -1.7069179299
        16 16 2 0  0.3759933035       -1.7072373106       -0.0694439840
        17 17 2 0 -1.7124296000        0.3336352522        0.1307959669
        18 18 2 0  0.9143159284        1.3089975397       -0.7151210582
        19 19 2 0 -0.3759920260        1.7072300336        0.0694634263
        20 20 2 0  1.7124281219       -0.3336312342       -0.1308207313
        21 21 2 0 -0.9143187026       -1.3089785474        0.7151290509
        22 22 2 0  0.9742085109        0.8837023041        1.1536069633
        23 23 2 0  0.7057104439       -0.6206907639        1.4756502961
        24 24 2 0  0.8745319670       -0.4648127187       -1.4422106957
        25 25 2 0 -1.1954804901       -0.6171923123       -0.1021449363
        26 26 2 0  0.0917363053       -1.0144887859       -0.8848410405
        27 27 2 0  0.9276243144       -0.8836123311        0.4234140820
        28 28 2 0  1.1954744473        0.6171883800        0.1021399054
        29 29 2 0 -0.9276176774        0.8836123556       -0.4234173533
        30 30 2 0 -0.3595942315       -0.4863167551        1.2061133825
        31 31 2 0  0.3595891589        0.4863295901       -1.2061152849
        32 32 2 0 -0.0917352078        1.0144694592        0.8848400639
        33 33 2 0  0.6410702480       -0.1978633363       -0.3898095439
        34 34 2 0 -0.4162942817       -0.0651798741       -0.6515502084
        35 35 2 0  0.1334019604        0.7474406294       -0.1600033264
        36 36 2 0 -0.6410732823        0.1978593218        0.3898012337
        37 37 2 0  0.4162968444        0.0651733322        0.6515490914
        38 38 2 0 -0.1333998872       -0.7474445984        0.1600019961
        ";

        let res: MolSystem = input.parse().unwrap();
        let energy = res.potential_enegry(lj_potential);
        assert_approx_eq!(energy, -173.928427, 1e-6);
        // panic!("Total energy: {}", energy);
    }
}
