use crate::data::core::ZERO_TOL;
use crate::symmop::SymmOp;
use crate::utils::normalize_frac_vectors;
use itertools::iproduct;
use nalgebra::Vector3;
use nalgebra::{try_invert_to, Matrix3};

#[derive(Debug)]
pub struct SpaceGroup {
    pub operations: Vec<SymmOp>,
    pub generators: Vec<SymmOp>,
}

impl SpaceGroup {
    pub fn from_generators(generators: &Vec<SymmOp>, frac_tols: &Vector3<f64>) -> Self {
        let symm_ops = Self::ops_from_generators(generators, frac_tols);

        return Self {
            operations: symm_ops,
            generators: generators.clone(),
        };
    }

    /// Transform symmetry operations `{R|t} -> {R'|t'}` according to `T^-1 R T = R'` and `T^-1 t = t'`.
    pub fn apply_transformation(
        &mut self,
        transformation_matrix: &Matrix3<f64>,
        frac_tols: &Vector3<f64>,
    ) {
        let mut inv_trans_mat: Matrix3<f64> = Matrix3::identity();
        let inverted = try_invert_to(transformation_matrix.clone(), &mut inv_trans_mat);

        if !inverted {
            panic!("Transformation matrix is not invertible!");
        }

        for op in self.generators.iter_mut() {
            let mut float_rot = op.rotation.cast::<f64>();

            float_rot = inv_trans_mat * float_rot * transformation_matrix;

            let mut float_trans = vec![inv_trans_mat * op.translation];
            normalize_frac_vectors(&mut float_trans, frac_tols);

            op.translation = float_trans[0];
            op.rotation = Matrix3::from_iterator(float_rot.iter().map(|&x| x.round() as i8));
        }

        //Find extra unit translations included in new cell if it is bigger (this is broken...)

        let volume_ratio = transformation_matrix.determinant().abs().round();

        let mut new_generators = self.generators.clone();

        while (new_generators.len() as f64 / self.generators.len() as f64)
            < (volume_ratio - ZERO_TOL)
        {
            let translation_vecs: Vec<Vector3<f64>> = vec![
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.0, 1.0, 0.0),
                Vector3::new(0.0, 0.0, 1.0),
                Vector3::new(1.0, 1.0, 0.0),
                Vector3::new(0.0, 1.0, 1.0),
                Vector3::new(1.0, 0.0, 1.0),
                Vector3::new(1.0, 1.0, 1.0),
            ];

            let original_ops = self.generators.clone();

            for op in original_ops.iter() {
                let mut found_translation_vecs: Vec<Vector3<f64>> = vec![op.translation];

                for vec in translation_vecs.iter() {
                    let mut trans_vec = vec![inv_trans_mat * (vec + op.translation)];
                    let mut unique = true;

                    for found_trans_vec in found_translation_vecs.iter() {
                        let mut diff = vec![trans_vec[0] - found_trans_vec];
                        normalize_frac_vectors(&mut diff, &frac_tols);

                        let translation_eq: bool = diff[0]
                            .iter()
                            .enumerate()
                            .all(|(i, val)| val.abs() < frac_tols[i]);

                        if translation_eq {
                            unique = false;
                        }
                    }

                    if unique {
                        normalize_frac_vectors(&mut trans_vec, &frac_tols);
                        let is_zero: bool = trans_vec[0]
                            .iter()
                            .enumerate()
                            .all(|(i, val)| val.abs() < frac_tols[i]);

                        if !is_zero {
                            found_translation_vecs.push(trans_vec[0]);

                            // Add new ops to group
                            let mut new_op = op.clone();
                            new_op.translation = trans_vec[0];
                            new_generators.push(new_op)
                        }
                    }
                }
            }
        }
        self.operations = Self::ops_from_generators(&new_generators, frac_tols);
    }

    fn ops_from_generators(generators: &Vec<SymmOp>, frac_tols: &Vector3<f64>) -> Vec<SymmOp> {
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
        return symm_ops;
    }
}
