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
