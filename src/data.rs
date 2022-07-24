use lazy_static::lazy_static;
use nalgebra::{Matrix3, Matrix4};
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::str::FromStr;
use std::string::String;

pub static ZERO_TOL: f32 = 1e-6; // Zero value tolerance

#[derive(PartialEq, Debug, Eq, Hash)]
pub enum Centering {
    P,
    C,
    I,
    F,
    R,
}

impl Display for Centering {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        match &self {
            Centering::P => write!(f, "P"),
            Centering::C => write!(f, "C"),
            Centering::I => write!(f, "I"),
            Centering::F => write!(f, "F"),
            Centering::R => write!(f, "R"),
        }
    }
}

impl FromStr for Centering {
    type Err = ();

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        match input {
            "P" => Ok(Centering::P),
            "C" => Ok(Centering::C),
            "I" => Ok(Centering::I),
            "F" => Ok(Centering::F),
            "R" => Ok(Centering::R),
            _ => Err(()),
        }
    }
}

pub enum LatticeSystem {
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Rhombohedral,
    Hexagonal,
    Cubic,
}

impl Display for LatticeSystem {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        match &self {
            LatticeSystem::Triclinic => write!(f, "Triclinic"),
            LatticeSystem::Monoclinic => write!(f, "Monoclinic"),
            LatticeSystem::Orthorhombic => write!(f, "Orthorhombic"),
            LatticeSystem::Tetragonal => write!(f, "Tetragonal"),
            LatticeSystem::Rhombohedral => write!(f, "Rhombohedral"),
            LatticeSystem::Hexagonal => write!(f, "Hexagonal"),
            LatticeSystem::Cubic => write!(f, "Cubic"),
        }
    }
}

impl FromStr for LatticeSystem {
    type Err = ();

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        match input {
            "Triclinic" => Ok(LatticeSystem::Triclinic),
            "Monoclinic" => Ok(LatticeSystem::Monoclinic),
            "Orthorhombic" => Ok(LatticeSystem::Orthorhombic),
            "Tetragonal" => Ok(LatticeSystem::Tetragonal),
            "Rhombohedral" => Ok(LatticeSystem::Rhombohedral),
            "Hexagonal" => Ok(LatticeSystem::Hexagonal),
            "Cubic" => Ok(LatticeSystem::Cubic),
            "a" => Ok(LatticeSystem::Triclinic),
            "m" => Ok(LatticeSystem::Monoclinic),
            "o" => Ok(LatticeSystem::Orthorhombic),
            "t" => Ok(LatticeSystem::Tetragonal),
            "r" => Ok(LatticeSystem::Rhombohedral),
            "h" => Ok(LatticeSystem::Hexagonal),
            "c" => Ok(LatticeSystem::Cubic),
            _ => Err(()),
        }
    }
}

#[derive(Debug, Clone, Copy)]
#[allow(non_camel_case_types)]
pub enum BravaisType {
    aP,
    mP,
    mC,
    mI,
    oP,
    oC,
    oI,
    oF,
    tP,
    tI,
    hR,
    hP,
    cP,
    cI,
    cF,
}

impl BravaisType {
    pub fn centering(&self) -> Centering {
        let centering: String = self.to_string().chars().last().unwrap().to_string();

        return Centering::from_str(&centering).unwrap();
    }

    pub fn lattice_system(&self) -> LatticeSystem {
        let system: String = self.to_string().chars().nth(0).unwrap().to_string();

        return LatticeSystem::from_str(&system).unwrap();
    }
}

impl Display for BravaisType {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        match &self {
            BravaisType::aP => write!(f, "aP"),
            BravaisType::mP => write!(f, "mP"),
            BravaisType::mC => write!(f, "mC"),
            BravaisType::mI => write!(f, "mI"),
            BravaisType::oP => write!(f, "oP"),
            BravaisType::oC => write!(f, "oC"),
            BravaisType::oI => write!(f, "oI"),
            BravaisType::oF => write!(f, "oF"),
            BravaisType::tP => write!(f, "tP"),
            BravaisType::tI => write!(f, "tI"),
            BravaisType::hR => write!(f, "hR"),
            BravaisType::hP => write!(f, "hP"),
            BravaisType::cP => write!(f, "cP"),
            BravaisType::cI => write!(f, "cI"),
            BravaisType::cF => write!(f, "cF"),
        }
    }
}

lazy_static! {
    // Bravais symbol lookup from lattice character number
    pub static ref LATTICE_CHAR_TO_BRAVAIS: HashMap<u8, BravaisType> = {
        HashMap::from([
            // First set, A=B=C
            (1, BravaisType::cF),
            (2, BravaisType::hR),
            (3, BravaisType::cP),
            (5, BravaisType::cI),
            (4, BravaisType::hR),
            (6, BravaisType::tI),
            (7, BravaisType::tI),
            (8, BravaisType::oI),
            // Second set, A=B
            (9, BravaisType::hR),
            (10, BravaisType::mC),
            (11, BravaisType::tP),
            (12, BravaisType::hP),
            (13, BravaisType::oC),
            (15, BravaisType::tI),
            (16, BravaisType::oF),
            (14, BravaisType::mC),
            (17, BravaisType::mC),
            // Third set, B=C
            (18, BravaisType::tI),
            (19, BravaisType::oI),
            (20, BravaisType::mC),
            (21, BravaisType::tP),
            (22, BravaisType::hP),
            (23, BravaisType::oC),
            (24, BravaisType::hR),
            (25, BravaisType::mC),
            // Fourth set
            (26, BravaisType::oF),
            (27, BravaisType::mC),
            (28, BravaisType::mC),
            (29, BravaisType::mC),
            (30, BravaisType::mC),
            (31, BravaisType::aP),
            (32, BravaisType::oP),
            (40, BravaisType::oC),
            (35, BravaisType::mP),
            (36, BravaisType::oC),
            (33, BravaisType::mP),
            (38, BravaisType::oC),
            (34, BravaisType::mP),
            (42, BravaisType::oI),
            (41, BravaisType::mC),
            (37, BravaisType::mC),
            (39, BravaisType::mC),
            (43, BravaisType::mI),
            (44, BravaisType::aP),
        ])
    };

    // Conventional cell transformation matrix from lattice character number
    pub static ref LATTICE_CHAR_TO_CONV_TRANS: HashMap<u8, Matrix3<i8>> = {
        // Note the ITC matrices (T) appear to be such that a transformation of a
        // lattice matrix (L) is given by:
        // L' = L * transpose(T)
        HashMap::from([
            // First set, A=B=C
            (1, Matrix3::new(1,-1,1,1,1,-1,-1,1,1)),
            (2, Matrix3::new(1,-1,0,-1,0,1,-1,-1,-1)),
            (3, Matrix3::identity()),
            (5, Matrix3::new(1,0,1,1,1,0,0,1,1)),
            (4, Matrix3::new(1,-1,0,-1,0,1,-1,-1,-1)),
            (6, Matrix3::new(0,1,1,1,0,1,1,1,0)),
            (7, Matrix3::new(1,0,1,1,1,0,0,1,1)),
            (8, Matrix3::new(-1,-1,0,-1,0,-1,0,-1,-1)),
            // Second set, A=B
            (9, Matrix3::new(1,0,0,-1,1,0,-1,-1,3)),
            (10, Matrix3::new(1,1,0,1,-1,0,0,0,-1)),
            (11, Matrix3::identity()),
            (12, Matrix3::identity()),
            (13, Matrix3::new(1,1,0,-1,1,0,0,0,1)),
            (15, Matrix3::new(1,0,0,0,1,0,1,1,2)),
            (16, Matrix3::new(-1,-1,0,1,-1,0,1,1,2)),
            (14, Matrix3::new(1,1,0,-1,1,0,0,0,1)),
            (17, Matrix3::new(1,-1,0,1,1,0,-1,0,-1)),
            // Third set, B=C
            (18, Matrix3::new(0,-1,1,1,-1,-1,1,0,0)),
            (19, Matrix3::new(-1,0,0,0,-1,1,-1,1,1)),
            (20, Matrix3::new(0,1,1,0,1,-1,-1,0,0)),
            (21, Matrix3::new(0,1,0,0,0,1,1,0,0)),
            (22, Matrix3::new(0,1,0,0,0,1,1,0,0)),
            (23, Matrix3::new(0,1,1,0,-1,1,1,0,0)),
            (24, Matrix3::new(1,2,1,0,-1,1,1,0,0)),
            (25, Matrix3::new(0,1,1,0,-1,1,1,0,0)),
            // Fourth set
            (26, Matrix3::new(1,0,0,-1,2,0,-1,0,2)),
            (27, Matrix3::new(-1,2,0,-1,0,0,0,-1,1)),
            (28, Matrix3::new(-1,0,0,-1,0,2,0,1,0)),
            (29, Matrix3::new(1,0,0,1,-2,0,0,0,-1)),
            (30, Matrix3::new(0,1,0,0,1,-2,-1,0,0)),
            (31, Matrix3::identity()),
            (32, Matrix3::identity()),
            (40, Matrix3::new(0,-1,0,0,1,2,-1,0,0)),
            (35, Matrix3::new(0,-1,0,-1,0,0,0,0,-1)),
            (36, Matrix3::new(1,0,0,-1,0,-2,0,1,0)),
            (33, Matrix3::identity()),
            (38, Matrix3::new(-1,0,0,1,2,0,0,0,-1)),
            (34, Matrix3::new(-1,0,0,0,0,-1,0,-1,0)),
            (42, Matrix3::new(-1,0,0,0,-1,0,1,1,2)),
            (41, Matrix3::new(0,-1,-2,0,-1,0,-1,0,0)),
            (37, Matrix3::new(1,0,2,1,0,0,0,1,0)),
            (39, Matrix3::new(-1,-2,0,-1,0,0,0,0,-1)),
            (43, Matrix3::new(-1,0,0,-1,-1,-2,0,-1,0)),
            (44, Matrix3::identity()),
        ])
    };

    // Conventional to primitive cell transformation lookup from centering
    pub static ref CENTERING_TO_PRIM_TRANS: HashMap<Centering, Matrix4<isize>> =
    HashMap::from([
        (
            Centering::P,
            Matrix4::identity()
        ),
        (
            Centering::C,
            Matrix4::new(1, 1, 0, 0, -1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2)
        ),
        (
            Centering::I,
            Matrix4::new(-1, 1, 1, 0, 1, -1, 1, 0, 1, 1, -1, 0, 0, 0, 0, 2)
        ),
        (
            Centering::F,
            Matrix4::new(0,1,1,0,1,0,1,0,1,1,0,0,0,0,0,2)
        ),
        (
            Centering::R,
            Matrix4::new(-1,2,-1,0,-2,1,1,0,1,1,1,0,0,0,0,3)
        ),
        ]);
}
