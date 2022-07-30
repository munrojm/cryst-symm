use crate::data::core::BravaisType;
use lazy_static::lazy_static;
use std::collections::HashMap;
use std::string::String;

lazy_static! {
    // International point group symbol from number
    pub static ref PG_NUM_TO_SYMBOL: HashMap<u8, String> = {
        HashMap::from([
            (1, "1".to_string()),
            (2, "-1".to_string()),
            (3, "-3".to_string()),
            (4, "-3m".to_string()),
            (5, "-4".to_string()),
            (6, "-42m".to_string()),
            (7, "-43m".to_string()),
            (8, "-6".to_string()),
            (9, "-6m2".to_string()),
            (10, "2".to_string()),
            (11, "2/m".to_string()),
            (12, "222".to_string()),
            (13, "23".to_string()),
            (14, "3".to_string()),
            (15, "32".to_string()),
            (16, "3m".to_string()),
            (17, "4".to_string()),
            (18, "4/m".to_string()),
            (19, "4/mmm".to_string()),
            (20, "422".to_string()),
            (21, "432".to_string()),
            (22, "4mm".to_string()),
            (23, "6".to_string()),
            (24, "6/m".to_string()),
            (25, "6/mmm".to_string()),
            (26, "622".to_string()),
            (27, "6mm".to_string()),
            (28, "m".to_string()),
            (29, "m-3".to_string()),
            (30, "m-3m".to_string()),
            (31, "mm2".to_string()),
            (32, "mmm".to_string()),
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

    // Encoded point group generator lookup from number
    pub static ref PG_NUM_TO_HOLOHEDRY_GEN: HashMap<u8, Vec<u16>> = {
        HashMap::from([
            (1, Vec::from([16484])),
            (2, Vec::from([16484, 3198])),
            (3, Vec::from([16484, 3198, 7817])),
            (4, Vec::from([16484, 7410, 3198, 7817])),
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
            (15, Vec::from([16484, 12270, 7817])),
            (16, Vec::from([16484, 7412, 7817])),
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
