use nalgebra::Vector3;

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
