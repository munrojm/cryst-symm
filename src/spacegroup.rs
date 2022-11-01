use crate::symmop::SymmOp;
use itertools::iproduct;
use nalgebra::Vector3;
use nalgebra::{try_invert_to, Matrix3};

pub struct SpaceGroup {
    pub operations: Vec<SymmOp>,
    pub generators: Vec<SymmOp>,
}

impl SpaceGroup {
    pub fn from_generators(generators: &Vec<SymmOp>, frac_tols: &Vector3<f64>) -> Self {
        let mut symm_ops: Vec<SymmOp> = generators.clone();
        let mut new_ops: Vec<SymmOp> = generators.clone();

        while new_ops.len() > 0 {
            let mut generated: Vec<SymmOp> = Vec::new();

            for (m1, m2) in iproduct!(symm_ops.clone(), new_ops) {
                let mut op = m1 * m2;

                if !(symm_ops
                    .iter()
                    .any(|symm_op| symm_op.is_approx_eq(&op, frac_tols)))
                {
                    op.normalize(frac_tols);
                    symm_ops.push(op.clone());
                    generated.push(op);
                }
            }

            new_ops = generated;
        }

        return Self {
            operations: symm_ops,
            generators: generators.clone(),
        };
    }

    /// Transform symmetry operations `{R|t} -> {R'|t'}` according to `T^-1 R T = R'` and `T^-1 t = t'`.
    pub fn apply_transformation(&mut self, transformation_matrix: &Matrix3<f64>) {
        let mut inv_trans_mat: Matrix3<f64> = Matrix3::identity();
        let inverted = try_invert_to(transformation_matrix.clone(), &mut inv_trans_mat);

        if !inverted {
            panic!("Transformation matrix is not invertible!");
        }

        for op in self.operations.iter_mut() {
            let mut float_rot = op.rotation.cast::<f64>();

            float_rot = inv_trans_mat * float_rot * transformation_matrix;
            op.translation = inv_trans_mat * op.translation;

            op.rotation = Matrix3::from_iterator(float_rot.iter().map(|&x| x as i8));
        }

        for op in self.generators.iter_mut() {
            let mut float_rot = op.rotation.cast::<f64>();

            float_rot = inv_trans_mat * float_rot * transformation_matrix;
            op.translation = inv_trans_mat * op.translation;

            op.rotation = Matrix3::from_iterator(float_rot.iter().map(|&x| x as i8));
        }
    }
}
