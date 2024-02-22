use std::str::FromStr;
use winnow::ascii::{digit1, float, multispace0, space0, till_line_ending};
use winnow::combinator::{opt, repeat};
use winnow::error::{AddContext, ErrMode, ErrorKind, ParserError, StrContext};
use winnow::stream::Stream;
use winnow::token::take_until;
use winnow::{PResult, Parser};

use crate::{Atom, MolSystem, MolSystems, Periodicity, Scaling};

fn number_of_atoms(input: &mut &str) -> PResult<u32> {
    "ITEM:".parse_next(input)?;
    space0.parse_next(input)?;
    "NUMBER OF ATOMS".parse_next(input)?;
    till_line_ending.parse_next(input)?;
    multispace0.parse_next(input)?;
    digit1.parse_to().parse_next(input)
}

fn parse_periodicity(input: &mut &str) -> PResult<Periodicity> {
    let expected = "pp pp pp";
    if input.len() < expected.len() {
        return Err(ErrMode::from_error_kind(input, ErrorKind::Slice));
    }
    let actual = input.next_slice(expected.len());
    if actual != expected {
        return Err(ErrMode::from_error_kind(input, ErrorKind::Verify));
    }
    Ok(Periodicity::Periodic)
}

fn parse_bounds_values(input: &mut &str) -> PResult<(Periodicity, f64, f64, f64, f64, f64, f64)> {
    space0.parse_next(input)?;
    let mut period = till_line_ending.parse_next(input)?;
    let periodicity = parse_periodicity(&mut period)?;
    till_line_ending.parse_next(input)?;
    multispace0.parse_next(input)?;

    let x_0 = float_parser(&mut take_until(0.., " ").parse_next(input)?)?;
    space0.parse_next(input)?;
    let x_1 = float_parser(&mut till_line_ending.parse_next(input)?)?;
    till_line_ending.parse_next(input)?;
    multispace0.parse_next(input)?;

    let y_0 = float_parser(&mut take_until(0.., " ").parse_next(input)?)?;
    space0.parse_next(input)?;
    let y_1 = float_parser(&mut till_line_ending.parse_next(input)?)?;
    till_line_ending.parse_next(input)?;
    multispace0.parse_next(input)?;

    let z_0 = float_parser(&mut take_until(0.., " ").parse_next(input)?)?;
    space0.parse_next(input)?;
    let z_1 = float_parser(&mut till_line_ending.parse_next(input)?)?;
    till_line_ending.parse_next(input)?;
    multispace0.parse_next(input)?;

    Ok((periodicity, x_0, x_1, y_0, y_1, z_0, z_1))
}

fn parse_full_bounds(input: &mut &str) -> PResult<(Periodicity, f64, f64, f64, f64, f64, f64)> {
    "ITEM:".parse_next(input)?;
    space0.parse_next(input)?;
    "BOX BOUNDS".parse_next(input)?;
    parse_bounds_values.parse_next(input)
}

fn parse_timestep_value(input: &mut &str) -> PResult<u32> {
    if opt("ITEM: TIMESTEP").parse_next(input)?.is_some() {
        till_line_ending.parse_next(input)?;
        multispace0.parse_next(input)?;
        let n = digit1.parse_to().parse_next(input)?;
        till_line_ending.parse_next(input)?;
        Ok(n)
    } else {
        Ok(0)
    }
}

fn float_parser(s: &mut &str) -> PResult<f64> {
    float(s)
}

fn parse_id_mol_type(input: &mut &str) -> PResult<(u32, u32, u8)> {
    let start = input.checkpoint();
    // space0.parse_next(input)?;
    let id: u32 = digit1.parse_to().parse_next(input).map_err(|err| {
        err.add_context(
            input,
            &start,
            StrContext::Expected(winnow::error::StrContextValue::Description("Valid Atom Id")),
        )
    })?;

    space0.parse_next(input)?;

    let mol: u32 = digit1.parse_to().parse_next(input).map_err(|err| {
        err.add_context(
            input,
            &start,
            StrContext::Expected(winnow::error::StrContextValue::Description("Valid Mol ID")),
        )
    })?;
    space0.parse_next(input)?;
    let at_type: u8 = digit1.parse_to().parse_next(input).map_err(|err| {
        err.add_context(
            input,
            &start,
            StrContext::Expected(winnow::error::StrContextValue::Description(
                "Valid Atom Type",
            )),
        )
    })?;
    Ok((id, mol, at_type))
}

fn parse_xyz(input: &mut &str) -> PResult<(f64, f64, f64)> {
    let start = input.checkpoint();
    let x = float_parser(&mut take_until(0.., " ").parse_next(input)?).map_err(|err| {
        err.add_context(
            input,
            &start,
            StrContext::Expected(winnow::error::StrContextValue::Description("Valid x pos")),
        )
    })?;

    space0.parse_next(input)?;
    let y = float_parser(&mut take_until(0.., " ").parse_next(input)?).map_err(|err| {
        err.add_context(
            input,
            &start,
            StrContext::Expected(winnow::error::StrContextValue::Description("Valid y pos")),
        )
    })?;

    space0.parse_next(input)?;
    let z = float_parser(&mut till_line_ending.parse_next(input)?).map_err(|err| {
        err.add_context(
            input,
            &start,
            StrContext::Expected(winnow::error::StrContextValue::Description("Valid z pos")),
        )
    })?;

    Ok((x, y, z))
}

fn parse_pos_vs(input: &mut &str) -> PResult<(f64, f64, f64, f64, f64, f64)> {
    let start = input.checkpoint();
    let x = float_parser(&mut take_until(0.., " ").parse_next(input)?).map_err(|err| {
        err.add_context(
            input,
            &start,
            StrContext::Expected(winnow::error::StrContextValue::Description("Valid x pos")),
        )
    })?;

    space0.parse_next(input)?;
    let y = float_parser(&mut take_until(0.., " ").parse_next(input)?).map_err(|err| {
        err.add_context(
            input,
            &start,
            StrContext::Expected(winnow::error::StrContextValue::Description("Valid y pos")),
        )
    })?;

    space0.parse_next(input)?;
    let z = float_parser(&mut take_until(0.., " ").parse_next(input)?).map_err(|err| {
        err.add_context(
            input,
            &start,
            StrContext::Expected(winnow::error::StrContextValue::Description("Valid z pos")),
        )
    })?;

    space0.parse_next(input)?;
    let vx = float_parser(&mut take_until(0.., " ").parse_next(input)?).map_err(|err| {
        err.add_context(
            input,
            &start,
            StrContext::Expected(winnow::error::StrContextValue::Description("Valid x vel")),
        )
    })?;

    space0.parse_next(input)?;
    let vy = float_parser(&mut take_until(0.., " ").parse_next(input)?).map_err(|err| {
        err.add_context(
            input,
            &start,
            StrContext::Expected(winnow::error::StrContextValue::Description("Valid y vel")),
        )
    })?;

    space0.parse_next(input)?;
    let vz = float_parser(&mut till_line_ending.parse_next(input)?).map_err(|err| {
        err.add_context(
            input,
            &start,
            StrContext::Expected(winnow::error::StrContextValue::Description("Valid z vel")),
        )
    })?;

    Ok((x, y, z, vx, vy, vz))
}

fn parse_atom(input: &mut &str) -> PResult<Atom> {
    let start = input.checkpoint();
    multispace0.parse_next(input)?;

    // space0.parse_next(input)?;
    let (id, mol, at_type) = parse_id_mol_type(input)?;

    space0.parse_next(input)?;
    let charge: f64 =
        float_parser(&mut take_until(0.., " ").parse_next(input)?).map_err(|err| {
            err.add_context(
                input,
                &start,
                StrContext::Expected(winnow::error::StrContextValue::Description("Valid Charge")),
            )
        })?;
    // .map_err(|err| err.append(&"Failed charge", ErrorKind::Fail))?;
    space0.parse_next(input)?;

    let (x, y, z) = parse_xyz(input)?;

    Ok(Atom {
        id,
        mol,
        at_type,
        charge,
        x,
        y,
        z,
        vx: 0.,
        vy: 0.,
        vz: 0.,
    })
}

fn parse_atom_uncharged(input: &mut &str) -> PResult<Atom> {
    multispace0.parse_next(input)?;

    let (id, mol, at_type) = parse_id_mol_type(input)?;

    space0.parse_next(input)?;
    let charge = 0.;

    let (x, y, z) = parse_xyz(input)?;

    Ok(Atom {
        id,
        mol,
        at_type,
        charge,
        x,
        y,
        z,
        vx: 0.,
        vy: 0.,
        vz: 0.,
    })
}

fn parse_atom_vel(input: &mut &str) -> PResult<Atom> {
    let start = input.checkpoint();
    multispace0.parse_next(input)?;

    // space0.parse_next(input)?;
    let (id, mol, at_type) = parse_id_mol_type(input)?;

    space0.parse_next(input)?;
    let charge: f64 =
        float_parser(&mut take_until(0.., " ").parse_next(input)?).map_err(|err| {
            err.add_context(
                input,
                &start,
                StrContext::Expected(winnow::error::StrContextValue::Description("Valid Charge")),
            )
        })?;
    // .map_err(|err| err.append(&"Failed charge", ErrorKind::Fail))?;
    space0.parse_next(input)?;

    let (x, y, z, vx, vy, vz) = parse_pos_vs(input)?;

    Ok(Atom {
        id,
        mol,
        at_type,
        charge,
        x,
        y,
        z,
        vx,
        vy,
        vz,
    })
}

fn parse_atom_unchaged_vel(input: &mut &str) -> PResult<Atom> {
    multispace0.parse_next(input)?;

    // space0.parse_next(input)?;
    let (id, mol, at_type) = parse_id_mol_type(input)?;

    space0.parse_next(input)?;
    let charge = 0.;

    let (x, y, z, vx, vy, vz) = parse_pos_vs(input)?;

    Ok(Atom {
        id,
        mol,
        at_type,
        charge,
        x,
        y,
        z,
        vx,
        vy,
        vz,
    })
}

fn parse_atom_list(input: &mut &str, n: usize, charge: bool, velocity: bool) -> PResult<Vec<Atom>> {
    if charge {
        if velocity {
            repeat(n..=n, parse_atom_vel)
                .context(winnow::error::StrContext::Label("Parsing atom"))
                .parse_next(input)
        } else {
            repeat(n..=n, parse_atom)
                .context(winnow::error::StrContext::Label("Parsing atom"))
                .parse_next(input)
        }
    } else {
        if velocity {
            repeat(n..=n, parse_atom_unchaged_vel)
                .context(winnow::error::StrContext::Label("Parsing atom"))
                .parse_next(input)
        } else {
            repeat(n..=n, parse_atom_uncharged)
                .context(winnow::error::StrContext::Label("Parsing atom"))
                .parse_next(input)
        }
    }
}

fn parse_scaling(input: &mut &str) -> PResult<Scaling> {
    let expected_scaled = "xs ys zs";
    let expected_unscaled = "xu yu zu";

    if input.len() < expected_scaled.len() {
        return Err(ErrMode::from_error_kind(input, ErrorKind::Slice));
    }

    let actual = input.next_slice(expected_scaled.len());

    if actual == expected_scaled {
        Ok(Scaling::Scaled)
    } else if actual == expected_unscaled {
        Ok(Scaling::Unscaled)
    } else {
        Err(ErrMode::from_error_kind(input, ErrorKind::Verify))
    }
}

fn vs_parse<'s>(input: &mut &'s str) -> PResult<&'s str> {
    space0.parse_next(input)?;
    "vx".parse_next(input)?;
    space0.parse_next(input)?;
    "vy".parse_next(input)?;
    space0.parse_next(input)?;
    "vz".parse_next(input)
}

fn parse_atoms_with_scaling(input: &mut &str, n: usize) -> PResult<(Scaling, Vec<Atom>)> {
    "ITEM:"
        .context(winnow::error::StrContext::Label("Parsing ITEM:"))
        .parse_next(input)?;
    space0.parse_next(input)?;
    "ATOMS"
        .context(winnow::error::StrContext::Label("Parsing ATOMS"))
        .parse_next(input)?;
    space0.parse_next(input)?;
    "id".context(winnow::error::StrContext::Label("Parsing id"))
        .parse_next(input)?;
    space0.parse_next(input)?;
    "mol"
        .context(winnow::error::StrContext::Label("Parsing mol"))
        .parse_next(input)?;
    space0.parse_next(input)?;
    "type"
        .context(winnow::error::StrContext::Label("Parsing type"))
        .parse_next(input)?;
    space0.parse_next(input)?;
    let charge = opt("q").parse_next(input)?.is_some();
    space0.parse_next(input)?;
    let scale = parse_scaling
        .context(winnow::error::StrContext::Label("Parsing scale"))
        .parse_next(input)?;
    let velocity = opt(vs_parse).parse_next(input)?.is_some();
    till_line_ending.parse_next(input)?;
    let atoms = parse_atom_list(input, n, charge, velocity)?;
    Ok((scale, atoms))
}

fn parse_full(input: &mut &str) -> PResult<MolSystem> {
    multispace0.parse_next(input)?;
    let time = parse_timestep_value
        .context(winnow::error::StrContext::Label("Parsing Time"))
        .parse_next(input)?;

    multispace0.parse_next(input)?;
    let n_atoms = number_of_atoms
        .context(winnow::error::StrContext::Label("Number of Atoms"))
        .parse_next(input)?;
    till_line_ending.parse_next(input)?;
    multispace0.parse_next(input)?;
    let (periodicity, x0, x1, y0, y1, z0, z1) = parse_full_bounds
        .context(winnow::error::StrContext::Label("Bounds"))
        .parse_next(input)?;
    multispace0.parse_next(input)?;

    let (scale, atoms) = parse_atoms_with_scaling(input, n_atoms as usize)?;

    Ok(MolSystem {
        time,
        n_atoms: n_atoms as usize,
        periodicity,
        x0,
        x1,
        y0,
        y1,
        z0,
        z1,
        scaling: scale,
        atoms,
    })
}

impl FromStr for MolSystem {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut input = s;
        parse_full(&mut input).map_err(|e| e.to_string())
    }
}

impl FromStr for MolSystems {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut input = s;
        Ok(MolSystems(
            repeat(0.., parse_full)
                .context(winnow::error::StrContext::Label("Parsing systems"))
                .parse_next(&mut input)
                .map_err(|e| e.to_string())?,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn timestep_test() {
        let mut input = "ITEM: TIMESTEP \n 1000";
        let time = parse_timestep_value(&mut input).unwrap();
        assert_eq!(time, 1000);
    }

    #[test]
    fn atom_test() {
        let mut input =
            "222 222 2 0 0.0463876 0.149607 0.143497 \n202 202 2 0 0.0455102 0.0531442 0.138979";
        let atom_res = parse_atom(&mut input).unwrap();
        let atom = Atom {
            id: 222,
            mol: 222,
            at_type: 2,
            charge: 0.,
            x: 0.046_387_6,
            y: 0.149_607,
            z: 0.143_497,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };

        assert_eq!(atom_res, atom);
        println!("{}", input);
        let atom_res = parse_atom(&mut input).unwrap();
        let atom = Atom {
            id: 202,
            mol: 202,
            at_type: 2,
            charge: 0.,
            x: 0.045_510_2,
            y: 0.053_144_2,
            z: 0.138_979,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        println!("{}", input);
        assert_eq!(atom_res, atom);
    }

    #[test]
    fn atom_list_test() {
        let mut input =
            "222 222 2 0 0.0463876 0.149607 0.143497\n202 202 2 0 0.0455102 0.0531442 0.138979\n";
        let atom_res = parse_atom_list(&mut input, 2, true, false).unwrap();
        let atom1 = Atom {
            id: 222,
            mol: 222,
            at_type: 2,
            charge: 0.,
            x: 0.046_387_6,
            y: 0.149_607,
            z: 0.143_497,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        let atom2 = Atom {
            id: 202,
            mol: 202,
            at_type: 2,
            charge: 0.,
            x: 0.045_510_2,
            y: 0.053_144_2,
            z: 0.138_979,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        assert_eq!(atom_res[0], atom1);
        assert_eq!(atom_res[1], atom2);
    }

    #[test]
    fn atom_test2() {
        let mut input = "222 222 2q 0 0.0463876 0.149607 0.143497\n";
        let atom_res = parse_atom(&mut input);
        assert!(atom_res.is_err());
        // panic!("{:?}", atom_res);
    }

    #[test]
    fn bounds_test() {
        let mut input = "ITEM: BOX BOUNDS pp pp pp\n5.3704472032230051e+00 4.9629552796755235e+01 #xmin xmax\n5.3704472032230051e+00 4.9629552796755235e+01 #ymin ymax\n5.3704472032230051e+00 4.9629552796755235e+01 #zmin zmax";
        let (periodicity, x0, x1, y0, y1, z0, z1) = parse_full_bounds(&mut input).unwrap();
        assert_eq!(periodicity, Periodicity::Periodic);

        assert_eq!(x0, 5.370_447_203_223_005);
        assert_eq!(x1, 4.962_955_279_675_523_5e1);
        assert_eq!(y0, 5.370_447_203_223_005);
        assert_eq!(y1, 4.962_955_279_675_523_5e1);
        assert_eq!(z0, 5.370_447_203_223_005);
        assert_eq!(z1, 4.962_955_279_675_523_5e1);
        // assert!(atom_res.is_err());
        // panic!("{:?}", atom_res);
    }

    #[test]
    fn n_atoms_test() {
        let mut input = "ITEM: NUMBER OF ATOMS\n 2000";
        let n = number_of_atoms(&mut input).unwrap();
        assert_eq!(n, 2000);
    }

    #[test]
    fn scaling_test() {
        let mut input = "xs ys zs";
        let scale = parse_scaling(&mut input).unwrap();
        assert_eq!(scale, Scaling::Scaled);

        let mut input = "xu yu zu";
        let scale = parse_scaling(&mut input).unwrap();
        assert_eq!(scale, Scaling::Unscaled);
    }

    #[test]
    fn atoms_parse_test() {
        let mut input = "ITEM: ATOMS id mol type xs ys zs
        222 222 2 0.0463876 0.149607 0.143497
        202 202 2 0.0455102 0.0531442 0.138979\n ";
        let (scale, atom_res) = parse_atoms_with_scaling(&mut input, 2).unwrap();
        // let scale = parse_scaling(&mut input).unwrap();
        assert_eq!(scale, Scaling::Scaled);

        let atom1 = Atom {
            id: 222,
            mol: 222,
            at_type: 2,
            charge: 0.,
            x: 0.046_387_6,
            y: 0.149_607,
            z: 0.143_497,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        let atom2 = Atom {
            id: 202,
            mol: 202,
            at_type: 2,
            charge: 0.,
            x: 0.045_510_2,
            y: 0.053_144_2,
            z: 0.138_979,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        assert_eq!(atom_res[0], atom1);
        assert_eq!(atom_res[1], atom2);

        // let mut input = "xu yu zu";
        // let scale = parse_scaling(&mut input).unwrap();
        // assert_eq!(scale, Scaling::Unscaled);
    }

    #[test]
    fn system_parse_test() {
        let mut input = "ITEM: TIMESTEP
        19520000
        ITEM: NUMBER OF ATOMS
        2
        ITEM: BOX BOUNDS pp pp pp
        3.4973054145025765e+02 4.0026945854974235e+02
        3.4973054145025765e+02 4.0026945854974235e+02
        1.1377176072643133e+03 1.3622823927356853e+03
        ITEM: ATOMS id mol type xs ys zs
        1 1 1 1.854626 0.593319 0.278698
        2 1 1 0.890698 0.645283 0.283455 ";
        let res = parse_full(&mut input).unwrap();

        assert_eq!(res.time, 19_520_000);
        assert_eq!(res.n_atoms, 2);
        assert_eq!(res.periodicity, Periodicity::Periodic);
        assert_eq!(res.x0, 3.497_305_414_502_576_5e2);
        assert_eq!(res.x1, 4.002_694_585_497_423_5e2);
        assert_eq!(res.y0, 3.497_305_414_502_576_5e2);
        assert_eq!(res.y1, 4.002_694_585_497_423_5e2);
        assert_eq!(res.z0, 1.137_717_607_264_313_3e3);
        assert_eq!(res.z1, 1.362_282_392_735_685_3e3);
        assert_eq!(res.scaling, Scaling::Scaled);
        let atom1 = Atom {
            id: 1,
            mol: 1,
            at_type: 1,
            charge: 0.,
            x: 1.854_626,
            y: 0.593_319,
            z: 0.278_698,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        let atom2 = Atom {
            id: 2,
            mol: 1,
            at_type: 1,
            charge: 0.,
            x: 0.890_698,
            y: 0.645_283,
            z: 0.283_455,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        assert_eq!(res.atoms[0], atom1);
        assert_eq!(res.atoms[1], atom2);
    }

    #[test]
    fn system_parse_fromstr_test() {
        let input = "ITEM: TIMESTEP
        0
        ITEM: NUMBER OF ATOMS
        2
        ITEM: BOX BOUNDS pp pp pp
        5.3704472032230051e+00 4.9629552796755235e+01 #xmin xmax
        5.3704472032230051e+00 4.9629552796755235e+01 #ymin ymax
        5.3704472032230051e+00 4.9629552796755235e+01 #zmin zmax
        ITEM: ATOMS id mol type q xs ys zs
        1655 1655 1 0 0.155146 0.00866324 0.0272779
1435 1435 1 0 0.0420433 0.117914 0.0418938 ";

        let res: MolSystem = input.parse().unwrap();

        assert_eq!(res.time, 0);
        assert_eq!(res.n_atoms, 2);
        assert_eq!(res.periodicity, Periodicity::Periodic);
        assert_eq!(res.x0, 5.370_447_203_223_005);
        assert_eq!(res.x1, 4.962_955_279_675_523_5e1);
        assert_eq!(res.y0, 5.370_447_203_223_005);
        assert_eq!(res.y1, 4.962_955_279_675_523_5e1);
        assert_eq!(res.z0, 5.370_447_203_223_005);
        assert_eq!(res.z1, 4.962_955_279_675_523_5e1);
        assert_eq!(res.scaling, Scaling::Scaled);
        let atom1 = Atom {
            id: 1655,
            mol: 1655,
            at_type: 1,
            charge: 0.,
            x: 0.155_146,
            y: 0.008_663_24,
            z: 0.027_277_9,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        let atom2 = Atom {
            id: 1435,
            mol: 1435,
            at_type: 1,
            charge: 0.,
            x: 0.042_043_3,
            y: 0.117_914,
            z: 0.041_893_8,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        assert_eq!(res.atoms[0], atom1);
        assert_eq!(res.atoms[1], atom2);
    }

    #[test]
    fn system_parse_vec_fromstr_test() {
        let input = "ITEM: TIMESTEP
        0
        ITEM: NUMBER OF ATOMS
        2
        ITEM: BOX BOUNDS pp pp pp
        5.3704472032230051e+00 4.9629552796755235e+01 #xmin xmax
        5.3704472032230051e+00 4.9629552796755235e+01 #ymin ymax
        5.3704472032230051e+00 4.9629552796755235e+01 #zmin zmax
        ITEM: ATOMS id mol type q xs ys zs
        1655 1655 1 0 0.155146 0.00866324 0.0272779
1435 1435 1 0 0.0420433 0.117914 0.0418938
ITEM: TIMESTEP
        0
        ITEM: NUMBER OF ATOMS
        2
        ITEM: BOX BOUNDS pp pp pp
        5.3704472032230051e+00 4.9629552796755235e+01 #xmin xmax
        5.3704472032230051e+00 4.9629552796755235e+01 #ymin ymax
        5.3704472032230051e+00 4.9629552796755235e+01 #zmin zmax
        ITEM: ATOMS id mol type q xs ys zs
        1655 1655 1 0 0.155146 0.00866324 0.0272779
1435 1435 1 0 0.0420433 0.117914 0.0418938  ";

        let res: MolSystems = input.parse().unwrap();
        let res = res.0;
        assert_eq!(res.len(), 2);

        assert_eq!(res[0].time, 0);
        assert_eq!(res[0].n_atoms, 2);
        assert_eq!(res[0].periodicity, Periodicity::Periodic);
        assert_eq!(res[0].x0, 5.370_447_203_223_005);
        assert_eq!(res[0].x1, 4.962_955_279_675_523_5e1);
        assert_eq!(res[0].y0, 5.370_447_203_223_005);
        assert_eq!(res[0].y1, 4.962_955_279_675_523_5e1);
        assert_eq!(res[0].z0, 5.370_447_203_223_005);
        assert_eq!(res[0].z1, 4.962_955_279_675_523_5e1);
        assert_eq!(res[0].scaling, Scaling::Scaled);
        let atom1 = Atom {
            id: 1655,
            mol: 1655,
            at_type: 1,
            charge: 0.,
            x: 0.155_146,
            y: 0.008_663_24,
            z: 0.027_277_9,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        let atom2 = Atom {
            id: 1435,
            mol: 1435,
            at_type: 1,
            charge: 0.,
            x: 0.042_043_3,
            y: 0.117_914,
            z: 0.041_893_8,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        assert_eq!(res[0].atoms[0], atom1);
        assert_eq!(res[0].atoms[1], atom2);
    }

    #[test]
    fn system_parse_fromstr_water_test() {
        let input = "ITEM: NUMBER OF ATOMS
        2 #number of atoms
        ITEM: BOX BOUNDS pp pp pp
        7.3598534393802213e+00 4.4401174149901216e+01 #xmin xmax
        7.3598534393802213e+00 4.4401174149901216e+01 #ymin ymax
        7.3598534393802213e+00 4.4401174149901216e+01 #zmin zmax
        ITEM: ATOMS id mol type q xs ys zs #q=charge, coordenates are scaled by the box side in each dimension (i.e. xmax-xmin)
        3985 1329 1 -1.1128 0.104458 0.13216 0.061756
        3986 1329 2 0.5564 0.0966852 0.13005 0.0372016";

        let res: MolSystem = input.parse().unwrap();

        assert_eq!(res.time, 0);
        assert_eq!(res.n_atoms, 2);
        assert_eq!(res.periodicity, Periodicity::Periodic);
        assert_eq!(res.x0, 7.3598534393802213e+00);
        assert_eq!(res.x1, 4.4401174149901216e+01);
        assert_eq!(res.y0, 7.3598534393802213e+00);
        assert_eq!(res.y1, 4.4401174149901216e+01);
        assert_eq!(res.z0, 7.3598534393802213e+00);
        assert_eq!(res.z1, 4.4401174149901216e+01);
        assert_eq!(res.scaling, Scaling::Scaled);
        let atom1 = Atom {
            id: 3985,
            mol: 1329,
            at_type: 1,
            charge: -1.1128,
            x: 0.104458,
            y: 0.13216,
            z: 0.061756,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        let atom2 = Atom {
            id: 3986,
            mol: 1329,
            at_type: 2,
            charge: 0.5564,
            x: 0.0966852,
            y: 0.13005,
            z: 0.0372016,
            vx: 0.,
            vy: 0.,
            vz: 0.,
        };
        assert_eq!(res.atoms[0], atom1);
        assert_eq!(res.atoms[1], atom2);
    }

    #[test]
    fn system_parse_fromstr_vel_test() {
        let input = "ITEM: TIMESTEP
        0
        ITEM: NUMBER OF ATOMS
        2
        ITEM: BOX BOUNDS pp pp pp
        7.7743652754171228e+02 8.1721825312885755e+02
        7.7743652754171228e+02 8.1721825312885755e+02
        7.7743652754171228e+02 8.1721825312885755e+02
        ITEM: ATOMS id mol type q xs ys zs vx vy vz
        1 229 1 0 0.214763 0.209893 0.0131017 -0.000827139 0.000588265 -0.000445509
        2 55 1 0 0.0626675 0.829244 0.0980052 -0.000348235 -0.00170769 0.000164523 ";

        let res: MolSystem = input.parse().unwrap();

        assert_eq!(res.time, 0);
        assert_eq!(res.n_atoms, 2);
        assert_eq!(res.periodicity, Periodicity::Periodic);
        assert_eq!(res.x0, 7.7743652754171228e+02);
        assert_eq!(res.x1, 8.1721825312885755e+02);
        assert_eq!(res.y0, 7.7743652754171228e+02);
        assert_eq!(res.y1, 8.1721825312885755e+02);
        assert_eq!(res.z0, 7.7743652754171228e+02);
        assert_eq!(res.z1, 8.1721825312885755e+02);
        assert_eq!(res.scaling, Scaling::Scaled);
        let atom1 = Atom {
            id: 1,
            mol: 229,
            at_type: 1,
            charge: 0.,
            x: 0.214763,
            y: 0.209893,
            z: 0.0131017,
            vx: -0.000827139,
            vy: 0.000588265,
            vz: -0.000445509,
        };
        let atom2 = Atom {
            id: 2,
            mol: 55,
            at_type: 1,
            charge: 0.,
            x: 0.0626675,
            y: 0.829244,
            z: 0.0980052,
            vx: -0.000348235,
            vy: -0.00170769,
            vz: 0.000164523,
        };
        assert_eq!(res.atoms[0], atom1);
        assert_eq!(res.atoms[1], atom2);
    }
}
