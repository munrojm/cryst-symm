use nalgebra::Vector3;
use std::f32::consts::PI;

pub fn normalize_frac_vectors(vecs: &mut Vec<Vector3<f32>>, frac_tols: &Vector3<f32>) {
    for vec in vecs {
        for i in 0..3 {
            vec[i] = vec[i] % 1.0;
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
    v1: &Vector3<f32>,
    v2: &Vector3<f32>,
    dv1: &f32,
    dv2: &f32,
    dtheta: &f32,
) -> f32 {
    let v1_norm = v1.magnitude();
    let v2_norm = v2.magnitude();
    let theta = (v1.dot(v2) / (v1_norm * v2_norm)).acos();

    let dd = ((v2_norm * theta.cos() * dv1).powf(2.0)
        + (v1_norm * theta.cos() * dv2).powf(2.0)
        + (v1_norm * v2_norm * theta.sin() * dtheta * (PI / 180.0)).powf(2.0))
    .sqrt();

    return dd;
}

pub fn num_negative_zero(vec: &Vec<f32>, epsilon: &f32) -> (i8, i8) {
    let mut num_negative = 0;
    let mut num_zero = 0;

    for val in vec {
        if val < &(-1.0 * epsilon) {
            num_negative += 1;
        } else if val < epsilon {
            num_zero += 1;
        }
    }
    return (num_negative, num_zero);
}

pub fn cust_eq(a: &f32, b: &f32, epsilon: &f32) -> bool {
    return !((a < &(b - epsilon)) || (b < &(a - epsilon)));
}
