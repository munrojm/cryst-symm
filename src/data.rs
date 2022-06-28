use lazy_static::lazy_static;
use nalgebra::{Matrix3, Vector6};
use std::collections::HashMap;
use std::string::String;

lazy_static! {
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
    pub static ref BRAVAIS_TO_TRANS: HashMap<String, Matrix3<isize>> = HashMap::from([
        (String::from("cI"), Matrix3::new(0, 1, 1, 1, 0, 1, 1, 1, 0),),
        (String::from("cF"), Matrix3::new(1, -1, 1, 1, 1, 1, 0, 0, 2),),
        (String::from("cP"), Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),),
        (String::from("cP"), Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),),
        (String::from("hP"), Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),),
        (
            String::from("hR"),
            Matrix3::new(1, 0, 1, -1, 1, 1, 0, -1, 1),
        ),
        (String::from("hR"), Matrix3::new(1, 0, 1, 0, 0, 3, 0, 1, 2),),
        (String::from("tI"), Matrix3::new(0, 1, 1, 1, 0, 1, 1, 1, 0),),
        (String::from("tI"), Matrix3::new(1, 0, 1, 0, 1, 1, 0, 0, 2),),
        (String::from("tP"), Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),),
        (String::from("tP"), Matrix3::new(1, 0, 0, 0, 0, 1, 0, 1, 1),),
        (String::from("tP"), Matrix3::new(0, 0, 1, 1, 1, 0, 0, 1, 0),),
        (String::from("oF"), Matrix3::new(1, -1, 1, 1, 1, 1, 0, 0, 2),),
        (String::from("oI"), Matrix3::new(0, 1, 1, 1, 0, 1, 1, 1, 0),),
        (String::from("oI"), Matrix3::new(1, 0, 1, 0, 1, 1, 0, 0, 2),),
        (String::from("oI"), Matrix3::new(0, 1, 1, 1, 0, 1, 1, 1, 0),),
        (String::from("oI"), Matrix3::new(1, 0, 1, 0, 1, 1, 0, 0, 2),),
        (
            String::from("o(AB)C"),
            Matrix3::new(2, 0, 0, 1, 1, 0, 0, 0, 1),
        ),
        (
            String::from("o(AB)C"),
            Matrix3::new(1, 1, 0, -1, 1, 0, 0, 0, 1),
        ),
        (String::from("oP"), Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),),
        (String::from("oP"), Matrix3::new(1, 0, 0, 0, 0, 1, 0, 1, 1),),
        (
            String::from("m(AC)I"),
            Matrix3::new(-1, 1, 0, -1, -1, 0, -1, 0, 1),
        ),
        (
            String::from("m(AC)I"),
            Matrix3::new(0, 1, -1, 1, 1, 0, 1, 0, -1),
        ),
        (
            String::from("m(AC)I"),
            Matrix3::new(-1, 0, 1, -1, 1, 0, -2, 0, 0),
        ),
        (
            String::from("m(AC)I"),
            Matrix3::new(0, 1, -1, 1, 1, 0, 1, 0, -1),
        ),
        (
            String::from("m(AC)I"),
            Matrix3::new(-1, 1, 0, -1, -1, 0, -1, 0, 1),
        ),
        (
            String::from("m(AC)I"),
            Matrix3::new(-1, 0, 1, -1, 1, 0, -2, 0, 0),
        ),
        (
            String::from("m(AC)I"),
            Matrix3::new(1, 0, -1, 1, -1, 0, 0, -1, -1),
        ),
        (String::from("mP"), Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),),
        (String::from("aP"), Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),),
        (String::from("aP"), Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),),
        (String::from("aP"), Matrix3::new(1, 0, 0, 0, 1, 0, 0, 0, 1),),
    ]);
}
