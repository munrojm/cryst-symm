use lazy_static::lazy_static;
use nalgebra::{Matrix3, Vector6};
use std::collections::HashMap;
use std::string::String;

lazy_static! {
    pub static ref SELLING_TO_BRAVAIS: HashMap<Vector6<usize>, (String, String, Matrix3<isize>)> = {
        HashMap::from([
            (
                Vector6::new(1, 1, 1, 1, 1, 1),
                (
                    String::from("Cubic"),
                    String::from("cI"),
                    Matrix3::new(0, 1, 1, 1, 0, 1, 1, 1, 0),
                ),
            ),
            (
                Vector6::new(0, 1, 1, 1, 1, 0),
                (
                    String::from("Cubic"),
                    String::from("cF"),
                    Matrix3::new(1, -1, 1, 1, 1, 1, 0, 0, 2),
                ),
            ),
            (
                Vector6::new(0, 0, 1, 1, 1, 0),
                (
                    String::from("Cubic"),
                    String::from("cP"),
                    Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),
                ),
            ),
            (
                Vector6::new(0, 0, 1, 0, 1, 1),
                (
                    String::from("Cubic"),
                    String::from("cP"),
                    Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),
                ),
            ),
            (
                Vector6::new(1, 0, 1, 0, 1, 2),
                (
                    String::from("Hexagonal"),
                    String::from("hP"),
                    Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),
                ),
            ),
            (
                Vector6::new(1, 1, 2, 1, 2, 2),
                (
                    String::from("Hexagonal"),
                    String::from("hR"),
                    Matrix3::new(1, 0, 1, -1, 1, 1, 0, -1, 1),
                ),
            ),
            (
                Vector6::new(0, 1, 1, 1, 2, 0),
                (
                    String::from("Hexagonal"),
                    String::from("hR"),
                    Matrix3::new(1, 0, 1, 0, 0, 3, 0, 1, 2),
                ),
            ),
            (
                Vector6::new(1, 2, 2, 2, 2, 1),
                (
                    String::from("Tetragonal"),
                    String::from("tI"),
                    Matrix3::new(0, 1, 1, 1, 0, 1, 1, 1, 0),
                ),
            ),
            (
                Vector6::new(0, 1, 1, 1, 1, 2),
                (
                    String::from("Tetragonal"),
                    String::from("tI"),
                    Matrix3::new(1, 0, 1, 0, 1, 1, 0, 0, 2),
                ),
            ),
            (
                Vector6::new(0, 0, 1, 0, 1, 2),
                (
                    String::from("Tetragonal"),
                    String::from("tP"),
                    Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),
                ),
            ),
            (
                Vector6::new(0, 0, 1, 1, 2, 0),
                (
                    String::from("Tetragonal"),
                    String::from("tP"),
                    Matrix3::new(1, 0, 0, 0, 0, 1, 0, 1, 1),
                ),
            ),
            (
                Vector6::new(0, 0, 1, 2, 0, 2),
                (
                    String::from("Tetragonal"),
                    String::from("tP"),
                    Matrix3::new(0, 0, 1, 1, 1, 0, 0, 1, 0),
                ),
            ),
            (
                Vector6::new(1, 2, 2, 2, 2, 3),
                (
                    String::from("Orthorhombic"),
                    String::from("oF"),
                    Matrix3::new(1, -1, 1, 1, 1, 1, 0, 0, 2),
                ),
            ),
            (
                Vector6::new(1, 2, 3, 3, 2, 1),
                (
                    String::from("Orthorhombic"),
                    String::from("oI"),
                    Matrix3::new(0, 1, 1, 1, 0, 1, 1, 1, 0),
                ),
            ),
            (
                Vector6::new(0, 1, 1, 2, 2, 3),
                (
                    String::from("Orthorhombic"),
                    String::from("oI"),
                    Matrix3::new(1, 0, 1, 0, 1, 1, 0, 0, 2),
                ),
            ),
            (
                Vector6::new(0, 1, 2, 2, 1, 0),
                (
                    String::from("Orthorhombic"),
                    String::from("oI"),
                    Matrix3::new(0, 1, 1, 1, 0, 1, 1, 1, 0),
                ),
            ),
            (
                Vector6::new(0, 1, 1, 2, 2, 0),
                (
                    String::from("Orthorhombic"),
                    String::from("oI"),
                    Matrix3::new(1, 0, 1, 0, 1, 1, 0, 0, 2),
                ),
            ),
            (
                Vector6::new(1, 0, 2, 0, 1, 3),
                (
                    String::from("Orthorhombic"),
                    String::from("o(AB)C"),
                    Matrix3::new(2, 0, 0, 1, 1, 0, 0, 0, 1),
                ),
            ),
            (
                Vector6::new(1, 0, 2, 0, 2, 3),
                (
                    String::from("Orthorhombic"),
                    String::from("o(AB)C"),
                    Matrix3::new(1, 1, 0, -1, 1, 0, 0, 0, 1),
                ),
            ),
            (
                Vector6::new(0, 0, 1, 0, 2, 3),
                (
                    String::from("Orthorhombic"),
                    String::from("oP"),
                    Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),
                ),
            ),
            (
                Vector6::new(0, 0, 1, 2, 3, 0),
                (
                    String::from("Orthorhombic"),
                    String::from("oP"),
                    Matrix3::new(1, 0, 0, 0, 0, 1, 0, 1, 1),
                ),
            ),
            (
                Vector6::new(1, 2, 3, 2, 3, 4),
                (
                    String::from("Monoclinic"),
                    String::from("m(AC)I"),
                    Matrix3::new(-1, 1, 0, -1, -1, 0, -1, 0, 1),
                ),
            ),
            (
                Vector6::new(1, 2, 3, 3, 2, 4),
                (
                    String::from("Monoclinic"),
                    String::from("m(AC)I"),
                    Matrix3::new(0, 1, -1, 1, 1, 0, 1, 0, -1),
                ),
            ),
            (
                Vector6::new(0, 1, 2, 3, 3, 4),
                (
                    String::from("Monoclinic"),
                    String::from("m(AC)I"),
                    Matrix3::new(-1, 0, 1, -1, 1, 0, -2, 0, 0),
                ),
            ),
            (
                Vector6::new(0, 1, 2, 2, 1, 3),
                (
                    String::from("Monoclinic"),
                    String::from("m(AC)I"),
                    Matrix3::new(0, 1, -1, 1, 1, 0, 1, 0, -1),
                ),
            ),
            (
                Vector6::new(0, 1, 2, 1, 2, 3),
                (
                    String::from("Monoclinic"),
                    String::from("m(AC)I"),
                    Matrix3::new(-1, 1, 0, -1, -1, 0, -1, 0, 1),
                ),
            ),
            (
                Vector6::new(0, 1, 2, 3, 3, 0),
                (
                    String::from("Monoclinic"),
                    String::from("m(AC)I"),
                    Matrix3::new(-1, 0, 1, -1, 1, 0, -2, 0, 0),
                ),
            ),
            (
                Vector6::new(0, 1, 2, 3, 1, 0),
                (
                    String::from("Monoclinic"),
                    String::from("m(AC)I"),
                    Matrix3::new(1, 0, -1, 1, -1, 0, 0, -1, -1),
                ),
            ),
            (
                Vector6::new(0, 1, 2, 0, 3, 4),
                (
                    String::from("Monoclinic"),
                    String::from("mP"),
                    Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),
                ),
            ),
            (
                Vector6::new(1, 2, 3, 4, 5, 6),
                (
                    String::from("Triclinic"),
                    String::from("aP"),
                    Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),
                ),
            ),
            (
                Vector6::new(0, 1, 2, 3, 4, 5),
                (
                    String::from("Triclinic"),
                    String::from("aP"),
                    Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),
                ),
            ),
            (
                Vector6::new(0, 1, 2, 3, 4, 0),
                (
                    String::from("Triclinic"),
                    String::from("aP"),
                    Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),
                ),
            ),
        ])
    };
}
