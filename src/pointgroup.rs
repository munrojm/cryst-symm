// PointGroup struct
// pub fn from_number
// fn generate_operations(generators)
// symbol, operations
use crate::data::pointgroup::{PG_NUM_TO_HOLOHEDRY_GEN, PG_NUM_TO_SYMBOL};
use crate::utils::decode;
use itertools::iproduct;
use nalgebra::Matrix3;
use std::collections::HashSet;
use std::iter::FromIterator;
use std::string::String;

pub struct PointGroup {
    pub symbol: String,
    pub operations: Vec<Matrix3<i8>>,
}

impl PointGroup {
    pub fn from_number(pg_num: &u8) -> Self {
        let symbol = PG_NUM_TO_SYMBOL.get(&pg_num).unwrap();
        let encoded_matrices = PG_NUM_TO_HOLOHEDRY_GEN.get(&pg_num).unwrap();

        let mut decoded_matrices: Vec<Matrix3<i8>> = Vec::new();

        for encoding in encoded_matrices.to_owned() {
            let matrix_vec = decode(encoding, 3, 1, 9);

            let matrix: Matrix3<i8> = Matrix3::from_iterator(matrix_vec).transpose();

            decoded_matrices.push(matrix);
        }

        let operations = Self::generate_operations(decoded_matrices);

        return Self {
            symbol: symbol.to_owned(),
            operations: operations,
        };
    }

    fn generate_operations(generators: Vec<Matrix3<i8>>) -> Vec<Matrix3<i8>> {
        let mut symm_ops: HashSet<Matrix3<i8>> = HashSet::from_iter(generators.clone());
        let mut new_ops: HashSet<Matrix3<i8>> = HashSet::from_iter(generators.clone());

        while new_ops.len() > 0 {
            let mut generated: HashSet<Matrix3<i8>> = HashSet::new();

            for (m1, m2) in iproduct!(symm_ops.to_owned(), new_ops.iter()) {
                let op = m1 * m2;
                if !(symm_ops.contains(&op)) {
                    symm_ops.insert(op);
                    generated.insert(op);
                }
            }

            new_ops = generated;
        }

        let final_vec: Vec<Matrix3<i8>> = Vec::from_iter(symm_ops);

        return final_vec;
    }
}
