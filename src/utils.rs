use nalgebra::{Matrix3, Vector3};
use std::f64::consts::PI;

///Normalizes fraction vectors components such that they are within `(-frac_tol, 1-frac_tol)`
pub fn normalize_frac_vectors(vecs: &mut Vec<Vector3<f64>>, frac_tols: &Vector3<f64>) {
    for vec in vecs {
        for i in 0..3 {
            vec[i] %= 1.0;
            if vec[i] < (-1.0 * frac_tols[i]) {
                vec[i] += 1.0;
            } else if vec[i] > (1.0 - frac_tols[i]) {
                vec[i] -= 1.0;
            };
        }
    }
}

///Calculate the uncertainty in the dot product of two vectors.
pub fn calculate_dot_uncertainty(
    v1: &Vector3<f64>,
    v2: &Vector3<f64>,
    dv1: &f64,
    dv2: &f64,
    dtheta: &f64,
) -> f64 {
    let v1_norm = v1.magnitude();
    let v2_norm = v2.magnitude();
    let theta = (v1.dot(v2) / (v1_norm * v2_norm)).acos();

    ((v2_norm * theta.cos() * dv1).powf(2.0)
        + (v1_norm * theta.cos() * dv2).powf(2.0)
        + (v1_norm * v2_norm * theta.sin() * dtheta * (PI / 180.0)).powf(2.0))
    .sqrt()
}

///Count the number of negative or zero components in a vector within a tolerance `epsilon`
pub fn num_negative_zero(vec: &Vec<f64>, epsilon: &f64) -> (i8, i8) {
    let mut num_negative = 0;
    let mut num_zero = 0;

    for val in vec {
        if val < &(-1.0 * epsilon) {
            num_negative += 1;
        } else if val < epsilon {
            num_zero += 1;
        }
    }
    (num_negative, num_zero)
}

///Custom equals evaluation to avoid floating point operations
pub fn cust_eq(a: &f64, b: &f64, epsilon: &f64) -> bool {
    !((a < &(b - epsilon)) || (b < &(a - epsilon)))
}

/// Decode a vector of integers which is encoded in a particular base.
/// Outputs data into the vector passed to `out`.
pub fn decode(e: u16, base: u8, sub: i8, len: i8) -> Vec<i8> {
    let mut out: Vec<i8> = Vec::new();
    let mut new_e = e;
    for _ in 0..len {
        out.insert(0, (new_e % base as u16) as i8 - sub);
        new_e /= base as u16
    }
    out
}

/// Function to compare two float rotation matrices element-wise with some tolerance
/// defined for each entry.
pub fn compare_matrices(a: &Matrix3<f64>, b: &Matrix3<f64>, tols: &[f64]) -> bool {
    for n in 0..9 {
        if !cust_eq(&a[n], &b[n], &tols[n]) {
            return false;
        }
    }

    true
}

/// Function to compare two float values within some tolerance
pub fn approx_equal(a: &f64, b: &f64, tol: &f64) -> bool {
    (a - b).abs() < *tol
}

pub fn approx_equal_iter(
    a: impl ExactSizeIterator<Item = f64>,
    b: impl ExactSizeIterator<Item = f64>,
    tol: &f64,
) -> bool {
    let a = a.into_iter();
    let b = b.into_iter();

    if a.len() != b.len() {
        panic!("Iterators provided for comparison do not have equal lengths")
    }
    a.zip(b).all(|(a_i, b_i)| approx_equal(&a_i, &b_i, tol))
}
