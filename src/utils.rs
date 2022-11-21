use nalgebra::{Matrix3, Vector3};

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
pub fn decode(e: u32, base: u8, sub: i8, len: i8) -> Vec<i8> {
    let mut out: Vec<i8> = Vec::new();
    let mut new_e = e;
    for _ in 0..len {
        out.insert(0, (new_e % base as u32) as i8 - sub);
        new_e /= base as u32
    }
    out
}

pub fn decode_spg_op(num: u32) -> (Vec<i8>, Vec<i8>) {
    let base8_translation_map = [0, 2, 3, 4, 6, 8, 9, 10];

    let num_rot = num % (3_u32.pow(9));
    let num_trans = num / (3_u32.pow(9));

    let dec_rot = decode(num_rot, 3, 1, 9);
    let dec_trans: Vec<i8> = decode(num_trans, 8, 0, 3)
        .iter()
        .map(|&i| base8_translation_map[i as usize])
        .collect();

    (dec_rot, dec_trans)
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

pub fn approx_equal_iter<'a, T: ExactSizeIterator<Item = &'a f64>>(a: T, b: T, tol: &f64) -> bool {
    let a = a.into_iter();
    let b = b.into_iter();

    if a.len() != b.len() {
        panic!("Iterators provided for comparison do not have equal lengths")
    }
    a.zip(b).all(|(a_i, b_i)| approx_equal(&a_i, &b_i, tol))
}

#[cfg(test)]

mod tests {

    use super::*;
    use nalgebra::{Matrix3, Vector3};

    #[test]
    fn test_approx_equal() {
        let base = [0.1, 0.2, 0.3];
        let good = [0.10001, 0.20001, 0.30001];
        let bad = [0.11, 0.21, 0.31];

        assert!(!approx_equal_iter(base.iter(), bad.iter(), &1e-3));
        assert!(approx_equal_iter(base.iter(), good.iter(), &1e-3));
    }

    #[test]
    fn test_num_negative() {
        let vec = vec![0.1, 0.2, 0.3, 1e-7, -0.1];

        let num = num_negative_zero(&vec, &1e-6);

        assert_eq!(num, (1, 1));
    }

    #[test]
    fn test_decode() {
        let num = 3360;
        let mat: Vec<i8> = vec![-1, 0, 0, 0, 1, 0, 0, 0, -1];

        let decode_mat = decode(num, 3, 1, 9);

        assert_eq!(mat, decode_mat);
    }

    #[test]
    fn test_custom_compare() {
        let mat1: Matrix3<f64> = Matrix3::new(1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0);
        let mat2: Matrix3<f64> =
            Matrix3::new(1.00001, 0.0, 1e-6, 0.0, -1.000001, 0.0, 0.0, 0.0, 1.0);
        let mat3: Matrix3<f64> = Matrix3::new(1.01, 0.0, 0.0, 0.0, -1.01, 0.0, 0.0, 0.0, 1.01);

        let tols: [f64; 9] = [
            0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
        ];

        assert!(!compare_matrices(&mat1, &mat3, &tols));
        assert!(compare_matrices(&mat1, &mat2, &tols));
    }
    #[test]
    fn test_normalize() {
        let mut vecs: Vec<Vector3<f64>> = vec![
            Vector3::new(-1.25, 1.000001, 1.47),
            Vector3::new(-1e-5, 0.47, -0.78),
        ];
        let frac_tols: Vector3<f64> = Vector3::new(1e-3, 1e-3, 1e-3);
        normalize_frac_vectors(&mut vecs, &frac_tols);

        assert!(approx_equal_iter(
            vecs[0].iter(),
            Vector3::new(0.75, 1.000001, 0.47).iter(),
            &1e-3
        ));
        assert!(approx_equal_iter(
            vecs[1].iter(),
            Vector3::new(-1e-5, 0.47, 0.22).iter(),
            &1e-3
        ));
    }
}
