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
        normalize_frac_vectors(&mut diff, &frac_tols);

        let translation_eq: bool = diff[0]
            .iter()
            .enumerate()
            .all(|(i, val)| val.abs() < frac_tols[i]);

        return (self.rotation == op.rotation) && translation_eq;
    }

    pub fn normalize(&mut self, frac_tols: &Vector3<f64>) {
        let mut new_translation = vec![self.translation];
        normalize_frac_vectors(&mut new_translation, &frac_tols);

        self.translation = new_translation[0];
    }
}

impl Mul for SymmOp {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let new_rotation = self.rotation * rhs.rotation;

        let float_rotation: Matrix3<f64> =
            Matrix3::from_iterator(rhs.rotation.iter().map(|&x| x as f64));

        let new_translation = (float_rotation * self.translation) + rhs.translation;

        return Self {
            rotation: new_rotation,
            translation: new_translation,
        };
    }
}

impl Mul<Vector3<f64>> for SymmOp {
    type Output = Vector3<f64>;

    fn mul(self, rhs: Vector3<f64>) -> Vector3<f64> {
        let float_rotation: Matrix3<f64> =
            Matrix3::from_iterator(self.rotation.iter().map(|&x| x as f64));
        return (float_rotation * rhs) + self.translation;
    }
}
