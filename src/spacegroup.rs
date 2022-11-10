use std::collections::HashMap;

use crate::data::core::ZERO_TOL;
use crate::symmop::SymmOp;
use crate::utils::normalize_frac_vectors;
use itertools::iproduct;
use nalgebra::Vector3;
use nalgebra::{try_invert_to, Matrix3};

#[derive(Debug)]
/// Represenation of a space group and its operations.
pub struct SpaceGroup {
    pub operations: Vec<SymmOp>,
}

impl SpaceGroup {
    /// Constructs a space group from a list of generators.
    ///
    /// # Arguments
    ///
    /// * `generators` - List of generator symmetry operations
    pub fn from_generators(generators: &[SymmOp], frac_tols: &Vector3<f64>) -> Self {
        let symm_ops = Self::ops_from_generators(generators, frac_tols);

        Self {
            operations: symm_ops,
        }
    }

    /// Transform symmetry operations `{R|t} -> {R'|t'}` according to `T^-1 R T = R'` and `T^-1 t = t'`.
    ///
    /// If a transformation results in a larger periodic cell, new operations will be searched for and added.
    ///
    /// # Arguments
    ///
    /// * `transformation_matrix` - Transformation matrix `T`
    /// * `frac_tols` - List of fractional tolerance values to use for each translation component comparison.
    pub fn apply_transformation(
        &mut self,
        transformation_matrix: &Matrix3<f64>,
        frac_tols: &Vector3<f64>,
    ) {
        //TODO: Use operations to smartly choose which translation vectors to use
        let mut inv_trans_mat: Matrix3<f64> = Matrix3::identity();
        let inverted = try_invert_to(*transformation_matrix, &mut inv_trans_mat);

        if !inverted {
            panic!("Transformation matrix is not invertible!");
        }

        // Find extra unit translations included in new cell if it is bigger

        let volume_ratio = (transformation_matrix.determinant().abs() * 100.0).round() / 100.0;

        let mut new_operations: Vec<SymmOp> = Vec::new();

        // Set up found vector mapping
        let mut found_vectors_mapping: HashMap<Matrix3<i8>, Vec<Vector3<f64>>> = HashMap::new();

        let mut initial_num_ops = self.operations.len() as f64;

        let mut mult = 2;
        while (new_operations.len() as f64 / initial_num_ops) < volume_ratio {
            let mut translation_vecs: Vec<Vector3<f64>> = Vec::new();

            // Generate unit translation vectors
            for (i, j, k) in iproduct!(0..mult, 0..mult, 0..mult) {
                translation_vecs.push(Vector3::new(i as f64, j as f64, k as f64))
            }

            for op in self.operations.iter() {
                // Calculate transformed rotation matrix
                let mut float_rot = op.rotation.cast::<f64>();
                float_rot = inv_trans_mat * float_rot * transformation_matrix;

                // Throw out non-integer rotation matrices
                if !float_rot
                    .iter()
                    .any(|val| (val.abs() < ZERO_TOL) || ((val.abs() - 1.0).abs() < ZERO_TOL))
                {
                    initial_num_ops -= 1.0;
                    continue;
                }

                let trans_rotation =
                    Matrix3::from_iterator(float_rot.iter().map(|&x| x.round() as i8));

                // Iterate through candidate transformed translation vector for a specific symm op
                // and check if it has been seen before. If not, add it to the list of new operations.
                for vec in translation_vecs.iter() {
                    let mut trans_vec = vec![inv_trans_mat * (op.translation + vec)];
                    normalize_frac_vectors(&mut trans_vec, frac_tols);

                    let mut unique = true;

                    for found_trans_vec in found_vectors_mapping
                        .entry(trans_rotation)
                        .or_default()
                        .iter()
                    {
                        let mut diff = vec![trans_vec[0] - found_trans_vec];
                        normalize_frac_vectors(&mut diff, frac_tols);

                        let translation_eq: bool = diff[0]
                            .iter()
                            .enumerate()
                            .all(|(i, val)| val.abs() < frac_tols[i]);

                        if translation_eq {
                            unique = false;
                            break;
                        }
                    }

                    if unique {
                        normalize_frac_vectors(&mut trans_vec, frac_tols);
                        found_vectors_mapping
                            .entry(trans_rotation)
                            .or_insert(Vec::new())
                            .push(trans_vec[0]);

                        // Add new ops to group
                        let mut new_op = op.clone();
                        new_op.translation = trans_vec[0];
                        new_op.rotation = trans_rotation;
                        new_operations.push(new_op)
                    }
                }
            }
            mult += 1;
        }
        self.operations = new_operations;
    }

    fn ops_from_generators(generators: &[SymmOp], frac_tols: &Vector3<f64>) -> Vec<SymmOp> {
        let mut symm_ops: Vec<SymmOp> = generators.to_vec();
        let mut new_ops: Vec<SymmOp> = generators.to_vec();
        while !new_ops.is_empty() {
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
        symm_ops
    }
}

#[cfg(test)]

mod tests {
    use super::*;
    use crate::{symmop::SymmOp, utils::approx_equal_iter};
    use nalgebra::Matrix3;

    #[test]
    fn test_from_generators() {
        let gen_rotations: [Matrix3<i8>; 3] = [
            Matrix3::identity(),
            Matrix3::new(1, 0, 0, 0, -1, 0, 0, 0, 1),
            Matrix3::identity(),
        ];

        let gen_translations: [Vector3<f64>; 3] = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.5),
            Vector3::new(0.5, 0.5, 0.0),
        ];

        let generators: Vec<SymmOp> = gen_rotations
            .iter()
            .zip(gen_translations)
            .map(|(r, t)| SymmOp {
                rotation: *r,
                translation: t,
            })
            .collect();

        let frac_tols = Vector3::new(1e-3, 1e-3, 1e-3);

        let sg = SpaceGroup::from_generators(&generators, &frac_tols);

        let correct_ops = [
            SymmOp {
                rotation: Matrix3::identity(),
                translation: Vector3::new(0.0, 0.0, 0.0),
            },
            SymmOp {
                rotation: Matrix3::identity(),
                translation: Vector3::new(0.5, 0.5, 0.0),
            },
            SymmOp {
                rotation: Matrix3::new(1, 0, 0, 0, -1, 0, 0, 0, 1),
                translation: Vector3::new(0.0, 0.0, 0.5),
            },
            SymmOp {
                rotation: Matrix3::new(1, 0, 0, 0, -1, 0, 0, 0, 1),
                translation: Vector3::new(0.5, 0.5, 0.5),
            },
        ];

        for op in correct_ops {
            let mut found = false;

            for sg_op in sg.operations.iter() {
                if sg_op.rotation == op.rotation
                    && approx_equal_iter(sg_op.translation.iter(), op.translation.iter(), &ZERO_TOL)
                {
                    found = true;
                }
            }

            assert!(found);
        }
    }

    #[test]
    fn test_apply_transformation() {
        let gen_rotations: [Matrix3<i8>; 3] = [
            Matrix3::identity(),
            Matrix3::new(1, 0, 0, 0, -1, 0, 0, 0, 1),
            Matrix3::identity(),
        ];

        let gen_translations: [Vector3<f64>; 3] = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.5),
            Vector3::new(0.5, 0.5, 0.0),
        ];

        let generators: Vec<SymmOp> = gen_rotations
            .iter()
            .zip(gen_translations)
            .map(|(r, t)| SymmOp {
                rotation: *r,
                translation: t,
            })
            .collect();

        let frac_tols = Vector3::new(1e-3, 1e-3, 1e-3);

        let mut sg = SpaceGroup::from_generators(&generators, &frac_tols);

        let transformation_matrix: Matrix3<f64> =
            Matrix3::new(2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

        sg.apply_transformation(&transformation_matrix, &frac_tols);

        let correct_transformed = [
            SymmOp {
                rotation: Matrix3::identity(),
                translation: Vector3::new(0.0, 0.0, 0.0),
            },
            SymmOp {
                rotation: Matrix3::identity(),
                translation: Vector3::new(0.25, 0.5, 0.0),
            },
            SymmOp {
                rotation: Matrix3::identity(),
                translation: Vector3::new(0.5, 0.0, 0.0),
            },
            SymmOp {
                rotation: Matrix3::identity(),
                translation: Vector3::new(0.75, 0.5, 0.0),
            },
            SymmOp {
                rotation: Matrix3::new(1, 0, 0, 0, -1, 0, 0, 0, 1),
                translation: Vector3::new(0.0, 0.0, 0.5),
            },
            SymmOp {
                rotation: Matrix3::new(1, 0, 0, 0, -1, 0, 0, 0, 1),
                translation: Vector3::new(0.75, 0.5, 0.5),
            },
            SymmOp {
                rotation: Matrix3::new(1, 0, 0, 0, -1, 0, 0, 0, 1),
                translation: Vector3::new(0.5, 0.0, 0.5),
            },
            SymmOp {
                rotation: Matrix3::new(1, 0, 0, 0, -1, 0, 0, 0, 1),
                translation: Vector3::new(0.75, 0.5, 0.5),
            },
        ];

        for symm_op in correct_transformed {
            let mut found = false;
            println!("{:?}", symm_op);

            for sg_op in sg.operations.iter() {
                if sg_op.rotation == symm_op.rotation
                    && approx_equal_iter(
                        sg_op.translation.iter(),
                        symm_op.translation.iter(),
                        &ZERO_TOL,
                    )
                {
                    found = true;
                    break;
                }
            }

            assert!(found);
        }
    }
}
