use crate::data::pointgroup::{PG_NUM_TO_GENERATOR_MATRICES, PG_NUM_TO_SYMBOL};
use crate::utils::decode;
use itertools::iproduct;
use nalgebra::{try_invert_to, Matrix3};
use std::collections::HashSet;
use std::iter::FromIterator;

pub struct PointGroup {
    pub symbol: &'static str,
    pub operations: Vec<Matrix3<i8>>,
    pub generators: Vec<Matrix3<i8>>,
}

impl PointGroup {
    pub fn from_number(pg_num: &u8) -> Self {
        let symbol = PG_NUM_TO_SYMBOL.get(pg_num).unwrap();
        let encoded_matrices = PG_NUM_TO_GENERATOR_MATRICES.get(pg_num).unwrap();

        let mut decoded_matrices: Vec<Matrix3<i8>> = Vec::new();

        for encoding in encoded_matrices {
            let matrix_vec = decode(*encoding, 3, 1, 9);

            let matrix: Matrix3<i8> = Matrix3::from_iterator(matrix_vec).transpose();

            decoded_matrices.push(matrix);
        }

        let operations = Self::generate_operations(&decoded_matrices);

        Self {
            symbol,
            operations,
            generators: decoded_matrices,
        }
    }

    /// Transform symmetry operations `R -> R'` according to `T^-1 R T = R'`.
    pub fn apply_transformation(&mut self, transformation_matrix: &Matrix3<f64>) {
        let mut inv_trans_mat: Matrix3<f64> = Matrix3::identity();
        let inverted = try_invert_to(*transformation_matrix, &mut inv_trans_mat);

        if !inverted {
            panic!("Transformation matrix is not invertible!");
        }

        for op in self.operations.iter_mut() {
            let mut float_op = op.cast::<f64>();
            float_op = inv_trans_mat * float_op * transformation_matrix;
            *op = Matrix3::from_iterator(float_op.iter().map(|&x| x.round() as i8));
        }

        for op in self.generators.iter_mut() {
            let mut float_op = op.cast::<f64>();
            float_op = inv_trans_mat * float_op * transformation_matrix;
            *op = Matrix3::from_iterator(float_op.iter().map(|&x| x.round() as i8));
        }
    }

    fn generate_operations(generators: &[Matrix3<i8>]) -> Vec<Matrix3<i8>> {
        let mut symm_ops: HashSet<Matrix3<i8>> = HashSet::from_iter(generators.to_owned());
        let mut new_ops: HashSet<Matrix3<i8>> = HashSet::from_iter(generators.to_owned());

        while !new_ops.is_empty() {
            let mut generated: HashSet<Matrix3<i8>> = HashSet::new();

            for (m1, m2) in iproduct!(symm_ops.clone().iter(), new_ops.iter()) {
                let op = m1 * m2;
                if !(symm_ops.contains(&op)) {
                    symm_ops.insert(op);
                    generated.insert(op);
                }
            }

            new_ops = generated;
        }

        let operations_vec: Vec<Matrix3<i8>> = Vec::from_iter(symm_ops);

        operations_vec
    }
}
