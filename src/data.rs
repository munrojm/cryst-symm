use lazy_static::lazy_static;
use nalgebra::{Matrix3, Matrix4, Vector6};
use std::collections::HashMap;
use std::string::String;

pub static ZERO_TOL: f32 = 1e-6; // Zero value tolerance

lazy_static! {
    // Bravais symbol lookup from integer Selling vector
    pub static ref SELLING_TO_BRAVAIS: HashMap<Vector6<usize>, String> = {
        HashMap::from([
            (Vector6::new(1, 1, 1, 1, 1, 1), String::from("cI")),
            (Vector6::new(0, 1, 1, 1, 1, 0), String::from("cF")),
            (Vector6::new(0, 0, 1, 1, 1, 0), String::from("cP")),
            (Vector6::new(0, 0, 1, 0, 1, 1), String::from("cP")),
            (Vector6::new(1, 0, 1, 0, 1, 2), String::from("hP")),
            (Vector6::new(1, 1, 2, 1, 2, 2), String::from("hR")),
            (Vector6::new(0, 1, 1, 1, 2, 0), String::from("hR")),
            (Vector6::new(1, 2, 2, 2, 2, 1), String::from("tI")),
            (Vector6::new(0, 1, 1, 1, 1, 2), String::from("tI")),
            (Vector6::new(0, 0, 1, 0, 1, 2), String::from("tP")),
            (Vector6::new(0, 0, 1, 1, 2, 0), String::from("tP")),
            (Vector6::new(0, 0, 1, 2, 0, 2), String::from("tP")),
            (Vector6::new(1, 2, 2, 2, 2, 3), String::from("oF")),
            (Vector6::new(1, 2, 3, 3, 2, 1), String::from("oI")),
            (Vector6::new(0, 1, 1, 2, 2, 3), String::from("oI")),
            (Vector6::new(0, 1, 2, 2, 1, 0), String::from("oI")),
            (Vector6::new(0, 1, 1, 2, 2, 0), String::from("oI")),
            (Vector6::new(1, 0, 2, 0, 1, 3), String::from("o(AB)C")),
            (Vector6::new(1, 0, 2, 0, 2, 3), String::from("o(AB)C")),
            (Vector6::new(0, 0, 1, 0, 2, 3), String::from("oP")),
            (Vector6::new(0, 0, 1, 2, 3, 0), String::from("oP")),
            (Vector6::new(1, 2, 3, 2, 3, 4), String::from("m(AC)I")),
            (Vector6::new(1, 2, 3, 3, 2, 4), String::from("m(AC)I")),
            (Vector6::new(0, 1, 2, 3, 3, 4), String::from("m(AC)I")),
            (Vector6::new(0, 1, 2, 2, 1, 3), String::from("m(AC)I")),
            (Vector6::new(0, 1, 2, 1, 2, 3), String::from("m(AC)I")),
            (Vector6::new(0, 1, 2, 3, 3, 0), String::from("m(AC)I")),
            (Vector6::new(0, 1, 2, 3, 1, 0), String::from("m(AC)I")),
            (Vector6::new(0, 1, 2, 0, 3, 4), String::from("mP")),
            (Vector6::new(1, 2, 3, 4, 5, 6), String::from("aP")),
            (Vector6::new(0, 1, 2, 3, 4, 5), String::from("aP")),
            (Vector6::new(0, 1, 2, 3, 4, 0), String::from("aP")),
        ])
    };
    // Conventional cell transformation lookup from integer Selling vector
    pub static ref SELLING_TO_CONV_TRANS: HashMap<Vector6<usize>, Matrix3<isize>> =
        HashMap::from([
            (
                Vector6::new(1, 1, 1, 1, 1, 1),
                Matrix3::new(0, 1, 1, 1, 0, 1, 1, 1, 0)
            ),
            (
                Vector6::new(0, 1, 1, 1, 1, 0),
                Matrix3::new(1, -1, 1, 1, 1, 1, 0, 0, 2)
            ),
            (
                Vector6::new(0, 0, 1, 1, 1, 0),
                Matrix3::new(1, 0, 0, 0, 0, 1, 0, 1, 1)
            ),
            (
                Vector6::new(0, 0, 1, 0, 1, 1),
                Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1)
            ),
            (
                Vector6::new(1, 0, 1, 0, 1, 2),
                Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1)
            ),
            (
                Vector6::new(1, 1, 2, 1, 2, 2),
                Matrix3::new(1, 0, 1, -1, 1, 1, 0, -1, 1)
            ),
            (
                Vector6::new(0, 1, 1, 1, 2, 0),
                Matrix3::new(1, 0, 1, 0, 0, 3, 0, 1, 2)
            ),
            (
                Vector6::new(1, 2, 2, 2, 2, 1),
                Matrix3::new(0, 1, 1, 1, 0, 1, 1, 1, 0)
            ),
            (
                Vector6::new(0, 1, 1, 1, 1, 2),
                Matrix3::new(1, 0, 1, 0, 1, 1, 0, 0, 2)
            ),
            (
                Vector6::new(0, 0, 1, 0, 1, 2),
                Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1)
            ),
            (
                Vector6::new(0, 0, 1, 1, 2, 0),
                Matrix3::new(1, 0, 0, 0, 0, 1, 0, 1, 1)
            ),
            (
                Vector6::new(0, 0, 1, 2, 0, 2),
                Matrix3::new(0, 0, 1, 1, 1, 0, 0, 1, 0)
            ),
            (
                Vector6::new(1, 2, 2, 2, 2, 3),
                Matrix3::new(1, -1, 1, 1, 1, 1, 0, 0, 2)
            ),
            (
                Vector6::new(1, 2, 3, 3, 2, 1),
                Matrix3::new(0, 1, 1, 1, 0, 1, 1, 1, 0)
            ),
            (
                Vector6::new(0, 1, 1, 2, 2, 3),
                Matrix3::new(1, 0, 1, 0, 1, 1, 0, 0, 2)
            ),
            (
                Vector6::new(0, 1, 2, 2, 1, 0),
                Matrix3::new(0, 1, 1, 1, 0, 1, 1, 1, 0)
            ),
            (
                Vector6::new(0, 1, 1, 2, 2, 0),
                Matrix3::new(1, 0, 1, 0, 1, 1, 0, 0, 2)
            ),
            (
                Vector6::new(1, 0, 2, 0, 1, 3),
                Matrix3::new(2, 0, 0, 1, 1, 0, 0, 0, 1)
            ),
            (
                Vector6::new(1, 0, 2, 0, 2, 3),
                Matrix3::new(1, 1, 0, -1, 1, 0, 0, 0, 1)
            ),
            (
                Vector6::new(0, 0, 1, 0, 2, 3),
                Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1)
            ),
            (
                Vector6::new(0, 0, 1, 2, 3, 0),
                Matrix3::new(1, 0, 0, 0, 0, 1, 0, 1, 1)
            ),
            (
                Vector6::new(1, 2, 3, 2, 3, 4),
                Matrix3::new(-1, 1, 0, -1, -1, 0, -1, 0, 1)
            ),
            (
                Vector6::new(1, 2, 3, 3, 2, 4),
                Matrix3::new(0, 1, -1, 1, 1, 0, 1, 0, -1)
            ),
            (
                Vector6::new(0, 1, 2, 3, 3, 4),
                Matrix3::new(-1, 0, 1, -1, 1, 0, -2, 0, 0)
            ),
            (
                Vector6::new(0, 1, 2, 2, 1, 3),
                Matrix3::new(0, 1, -1, 1, 1, 0, 1, 0, -1)
            ),
            (
                Vector6::new(0, 1, 2, 1, 2, 3),
                Matrix3::new(-1, 1, 0, -1, -1, 0, -1, 0, 1)
            ),
            (
                Vector6::new(0, 1, 2, 3, 3, 0),
                Matrix3::new(-1, 0, 1, -1, 1, 0, -2, 0, 0)
            ),
            (
                Vector6::new(0, 1, 2, 3, 1, 0),
                Matrix3::new(1, 0, -1, 1, -1, 0, 0, -1, -1)
            ),
            (
                Vector6::new(0, 1, 2, 0, 3, 4),
                Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1)
            ),
            (
                Vector6::new(1, 2, 3, 4, 5, 6),
                Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1)
            ),
            (
                Vector6::new(0, 1, 2, 3, 4, 5),
                Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1)
            ),
            (
                Vector6::new(0, 1, 2, 3, 4, 0),
                Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1)
            ),
        ]);

    // Conventional to primitive cell transformation lookup from centering
    pub static ref CENTERING_TO_PRIM_TRANS: HashMap<String, Matrix4<isize>> =
    HashMap::from([
        (
            String::from("P"),
            Matrix4::identity()
        ),
        (
            String::from("C"),
            Matrix4::new(1, 1, 0, 0, -1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2)
        ),
        (
            String::from("S"),
            Matrix4::new(1, 1, 0, 0, -1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2)
        ),
        (
            String::from("I"),
            Matrix4::new(-1, -1, 1, 0, 1, -1, 1, 0, 1, 1, -1, 0, 0, 0, 0, 2)
        ),
        (
            String::from("F"),
            Matrix4::new(0,1,1,0,1,0,1,0,1,1,0,0,0,0,0,2)
        ),
        (
            String::from("R"),
            Matrix4::new(-1,2,-1,0,-2,1,1,0,1,1,1,0,0,0,0,3)
        ),
        ]);
}
