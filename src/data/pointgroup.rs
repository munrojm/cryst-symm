use crate::data::core::BravaisType;
use lazy_static::lazy_static;
use std::collections::HashMap;
use std::string::String;

// {'-1': [3198],
// '-3': [3198, 7817],
// '-3m': [7410, 3198, 7817],
// '-4': [11784],
// '-42m': [3360, 11784],
// '-43m': [10816, 11784, 12272],
// '-6': [16482, 7817],
// '-6m2': [16482, 7412, 7817],
// '1': [16484],
// '2': [3360],
// '2/m': [3360, 3198],
// '222': [3200, 3360],
// '23': [3360, 10816],
// '3': [7817],
// '32': [12270, 7817],
// '3m': [7412, 7817],
// '4': [7898],
// '4/m': [7898, 3198],
// '4/mmm': [3360, 7898, 3198],
// '422': [3360, 7898],
// '432': [10816, 7898],
// '4mm': [7898, 16322],
// '6': [3200, 7817],
// '6/m': [3200, 3198, 7817],
// '6/mmm': [3200, 12270, 7817, 3198],
// '622': [3200, 12270, 7817],
// '6mm': [3200, 7412, 7817],
// 'm': [16322],
// 'm-3': [3360, 10816, 3198],
// 'm-3m': [10816, 7898, 3198, 12272],
// 'mm2': [3200, 16322],
// 'mmm': [3200, 3360, 3198]}

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
            (11, "2".to_string()),
            (12, "2/m".to_string()),
            (13, "222".to_string()),
            (14, "23".to_string()),
            (15, "3".to_string()),
            (16, "32".to_string()),
            (17, "3m".to_string()),
            (18, "4".to_string()),
            (19, "4/m".to_string()),
            (20, "4/mmm".to_string()),
            (21, "422".to_string()),
            (22, "432".to_string()),
            (23, "4mm".to_string()),
            (24, "6".to_string()),
            (25, "6/m".to_string()),
            (26, "6/mmm".to_string()),
            (27, "622".to_string()),
            (28, "6mm".to_string()),
            (29, "m".to_string()),
            (30, "m-3".to_string()),
            (31, "m-3m".to_string()),
            (32, "mm2".to_string()),
            (33, "mmm".to_string()),
        ])
    };

    // Holohedry point group number lookup from Bravais symbol
    pub static ref BRAVAIS_TO_HOLOHEDRY: HashMap<BravaisType, u8> = {
        HashMap::from([
            (BravaisType::aP, 2),
            (BravaisType::mP, 12),
            (BravaisType::mC, 12),
            (BravaisType::mI, 12),
            (BravaisType::oP, 33),
            (BravaisType::oC, 33),
            (BravaisType::oI, 33),
            (BravaisType::oF, 33),
            (BravaisType::tP, 20),
            (BravaisType::tI, 20),
            (BravaisType::hR, 4),
            (BravaisType::hP, 26),
            (BravaisType::cP, 31),
            (BravaisType::cI, 31),
            (BravaisType::cF, 31)
        ])
    };

    // Encoded point group generator lookup from number
    pub static ref PG_NUM_TO_HOLOHEDRY: HashMap<u8, Vec<u16>> = {
        HashMap::from([
            (1, Vec::from([3198])),
            (2, Vec::from([3198, 7817])),
            (3, Vec::from([7410, 3198, 7817])),
            (4, Vec::from([11784])),
            (5, Vec::from([3360, 11784])),
            (6, Vec::from([10816, 11784, 12272])),
            (7, Vec::from([16482, 7817])),
            (8, Vec::from([16482, 7412, 7817])),
            (10, Vec::from([3360])),
            (11, Vec::from([3360, 3198])),
            (12, Vec::from([3200, 3360])),
            (13, Vec::from([3360, 10816])),
            (14, Vec::from([7817])),
            (15, Vec::from([12270, 7817])),
            (16, Vec::from([7412, 7817])),
            (17, Vec::from([7898])),
            (18, Vec::from([7898, 3198])),
            (19, Vec::from([3360, 7898, 3198])),
            (20, Vec::from([3360, 7898])),
            (21, Vec::from([10816, 7898])),
            (22, Vec::from([7898, 16322])),
            (23, Vec::from([3200, 7817])),
            (24, Vec::from([3200, 3198, 7817])),
            (25, Vec::from([3200, 12270, 7817, 3198])),
            (26, Vec::from([3200, 12270, 7817])),
            (27, Vec::from([3200, 7412, 7817])),
            (28, Vec::from([16322])),
            (29, Vec::from([3360, 10816, 3198])),
            (30, Vec::from([10816, 7898, 3198, 12272])),
            (31, Vec::from([3200, 16322])),
            (32, Vec::from([3200, 3360, 3198])),
        ])
    };
}
