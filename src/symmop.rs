use crate::utils::normalize_frac_vectors;
use nalgebra::{Matrix3, Vector3};
use std::ops::Mul;

#[derive(Debug, Clone)]
pub struct SymmOp {
    pub rotation: Matrix3<i8>,
    pub translation: Vector3<f64>,
}

impl SymmOp {
    pub fn is_approx_eq(&self, op: &Self, frac_tols: &Vector3<f64>) -> bool {
        let mut diff = vec![self.translation - op.translation];
        normalize_frac_vectors(&mut diff, frac_tols);

        let translation_eq: bool = diff[0]
            .iter()
            .enumerate()
            .all(|(i, val)| val.abs() < frac_tols[i]);

        (self.rotation == op.rotation) && translation_eq
    }

    pub fn normalize(&mut self, frac_tols: &Vector3<f64>) {
        let mut new_translation = vec![self.translation];
        normalize_frac_vectors(&mut new_translation, frac_tols);

        self.translation = new_translation[0];
    }
}

impl Mul for SymmOp {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let new_rotation = self.rotation * rhs.rotation;

        let float_rotation: Matrix3<f64> = rhs.rotation.cast::<f64>();

        let new_translation = (float_rotation * self.translation) + rhs.translation;

        Self {
            rotation: new_rotation,
            translation: new_translation,
        }
    }
}

impl Mul<Vector3<f64>> for SymmOp {
    type Output = Vector3<f64>;

    fn mul(self, rhs: Vector3<f64>) -> Vector3<f64> {
        let float_rotation: Matrix3<f64> = self.rotation.cast::<f64>();

        (float_rotation * rhs) + self.translation
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::approx_equal_iter;
    use nalgebra::{Matrix3, Vector3};

    fn generate_symm_ops() -> (SymmOp, SymmOp) {
        let rotation: Matrix3<i8> = Matrix3::new(0, 0, 1, 0, -1, 0, 1, 0, 0);
        let translation: Vector3<f64> = Vector3::new(0.0, 0.5, -0.75);
        let symm_op_1 = SymmOp {
            rotation: rotation,
            translation: translation,
        };

        let rotation: Matrix3<i8> = Matrix3::new(0, 1, 0, 1, 0, 0, 0, 0, 1);
        let translation: Vector3<f64> = Vector3::new(0.0, 0.50001, -0.74999);
        let symm_op_2 = SymmOp {
            rotation: rotation,
            translation: translation,
        };

        (symm_op_1, symm_op_2)
    }

    #[test]
    fn test_symm_ops() {
        let (mut symm_op_1, symm_op_2) = generate_symm_ops();

        let frac_tols = Vector3::new(1e-3, 1e-3, 1e-3);

        symm_op_1.normalize(&frac_tols);

        assert!(approx_equal_iter(
            symm_op_1.translation.iter(),
            Vector3::new(0.0, 0.5, 0.25).iter(),
            &1e-4
        ));

        assert!(symm_op_1.is_approx_eq(&symm_op_2, &frac_tols));
    }

    #[test]
    fn test_mult() {
        let (symm_op_1, symm_op_2) = generate_symm_ops();

        let symm_op_3 = symm_op_1.clone();

        let symm_op_mult = symm_op_1 * symm_op_2;

        assert!(approx_equal_iter(
            symm_op_mult.translation.iter(),
            Vector3::new(0.5, 0.5, -1.5).iter(),
            &1e-4
        ));

        let test_mat: Matrix3<i8> = Matrix3::new(0, 0, 1, -1, 0, 0, 0, 1, 0);

        assert_eq!(test_mat, symm_op_mult.rotation);

        let mult_vec = symm_op_3 * Vector3::new(0.5, 0.25, -0.75);

        println!("{:?}", mult_vec);

        assert!(approx_equal_iter(
            mult_vec.iter(),
            Vector3::new(-0.75, 0.25, -0.25).iter(),
            &1e-4
        ));
    }
}
