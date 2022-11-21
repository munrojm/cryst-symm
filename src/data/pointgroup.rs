use crate::data::core::BravaisType;
use lazy_static::lazy_static;
use std::collections::HashMap;

lazy_static! {
    // International point group symbol from number
    pub static ref PG_NUM_TO_SYMBOL: HashMap<u8, &'static str> = {
        HashMap::from([
            (1, "1"),
            (2, "-1"),
            (3, "-3"),
            (4, "-3m"),
            (5, "-4"),
            (6, "-42m"),
            (7, "-43m"),
            (8, "-6"),
            (9, "-6m2"),
            (10, "2"),
            (11, "2/m"),
            (12, "222"),
            (13, "23"),
            (14, "3"),
            (15, "32"),
            (16, "3m"),
            (17, "4"),
            (18, "4/m"),
            (19, "4/mmm"),
            (20, "422"),
            (21, "432"),
            (22, "4mm"),
            (23, "6"),
            (24, "6/m"),
            (25, "6/mmm"),
            (26, "622"),
            (27, "6mm"),
            (28, "m"),
            (29, "m-3"),
            (30, "m-3m"),
            (31, "mm2"),
            (32, "mmm"),
        ])
    };

    // Holohedry point group number lookup from Bravais symbol
    pub static ref BRAVAIS_TO_HOLOHEDRY_NUM: HashMap<BravaisType, u8> = {
        HashMap::from([
            (BravaisType::aP, 2),
            (BravaisType::mP, 11),
            (BravaisType::mC, 11),
            (BravaisType::mI, 11),
            (BravaisType::oP, 32),
            (BravaisType::oC, 32),
            (BravaisType::oI, 32),
            (BravaisType::oF, 32),
            (BravaisType::tP, 19),
            (BravaisType::tI, 19),
            (BravaisType::hR, 4),
            (BravaisType::hP, 25),
            (BravaisType::cP, 30),
            (BravaisType::cI, 30),
            (BravaisType::cF, 30)
        ])
    };

    // Encoded conventional setting point group generators lookup from number
    pub static ref PG_NUM_TO_GENERATOR_MATRICES: HashMap<u8, Vec<u32>> = {
        HashMap::from([
            (1, Vec::from([16484])),
            (2, Vec::from([16484, 3198])),
            (3, Vec::from([16484, 7817, 3198])),
            (4, Vec::from([16484, 7817, 12270, 3198])),
            (5, Vec::from([16484, 11784])),
            (6, Vec::from([16484, 3360, 11784])),
            (7, Vec::from([16484, 10816, 11784, 12272])),
            (8, Vec::from([16484, 16482, 7817])),
            (9, Vec::from([16484, 16482, 7412, 7817])),
            (10, Vec::from([16484, 3360])),
            (11, Vec::from([16484, 3360, 3198])),
            (12, Vec::from([16484, 3200, 3360])),
            (13, Vec::from([16484, 3360, 10816])),
            (14, Vec::from([16484, 7817])),
            (15, Vec::from([16484, 7817, 12270])),
            (16, Vec::from([16484, 7817, 7412])),
            (17, Vec::from([16484, 7898])),
            (18, Vec::from([16484, 7898, 3198])),
            (19, Vec::from([16484, 3360, 7898, 3198])),
            (20, Vec::from([16484, 3360, 7898])),
            (21, Vec::from([16484, 10816, 7898])),
            (22, Vec::from([16484, 7898, 16322])),
            (23, Vec::from([16484, 3200, 7817])),
            (24, Vec::from([16484, 3200, 3198, 7817])),
            (25, Vec::from([16484, 3200, 12270, 7817, 3198])),
            (26, Vec::from([16484, 3200, 12270, 7817])),
            (27, Vec::from([16484, 3200, 7412, 7817])),
            (28, Vec::from([16484, 16322])),
            (29, Vec::from([16484, 3360, 10816, 3198])),
            (30, Vec::from([16484, 10816, 7898, 3198, 12272])),
            (31, Vec::from([16484, 3200, 16322])),
            (32, Vec::from([16484, 3200, 3360, 3198])),
        ])
    };
}
