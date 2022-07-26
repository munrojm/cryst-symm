use crate::data::core::{
    CENTERING_TO_PRIM_TRANS, LATTICE_CHAR_TO_BRAVAIS, LATTICE_CHAR_TO_CONV_TRANS,
};
use crate::reduce::Reducer;
use crate::structure::Structure;
use crate::utils::{cust_eq, num_negative_zero};
use nalgebra::{Matrix3, Matrix4, Vector3};

#[derive(Debug, Clone)]
pub struct SymmetryAnalyzer {
    pub dtol: f32,
    pub atol: f32,
}

impl SymmetryAnalyzer {
    /// Obtains the crystallographic primitive crystal structure
    pub fn get_standard_primitive_structure(&self, structure: &Structure) -> Structure {
        let reducer = Reducer {
            dtol: self.dtol,
            atol: self.atol,
        };

        let prim_structure = reducer.find_primitive_cell(structure);

        let mut reduced_structure = reducer.niggli_reduce(&prim_structure, &1e-5);

        let lattice_character = self.get_lattice_character(&reduced_structure, &1e-5);

        let trans_mat = LATTICE_CHAR_TO_CONV_TRANS.get(&lattice_character);

        match trans_mat {
            Some(mat) => {
                let float_mat = Matrix3::from_iterator(mat.iter().map(|&x| x as f32));
                reduced_structure.apply_transformation(&float_mat.transpose(), &self.dtol);
            }
            None => panic!("Could not find the appropriate conventional transformation matrix!"),
        }

        let bravais_symbol = LATTICE_CHAR_TO_BRAVAIS.get(&lattice_character);

        match bravais_symbol {
            Some(symbol) => {
                let prim_trans_mat = CENTERING_TO_PRIM_TRANS.get(&symbol.centering());
                match prim_trans_mat {
                    Some(mat) => {
                        let float_mat: Matrix4<f32> =
                            Matrix4::from_iterator(mat.iter().map(|&x| x as f32));

                        let prim_trans_mat: Matrix3<f32> =
                            Matrix3::from(float_mat.fixed_slice::<3, 3>(0, 0)) / float_mat.m44;
                        reduced_structure.apply_transformation(&prim_trans_mat, &self.dtol);
                    }
                    None => {
                        panic!("Could not find the appropriate primitive transformation matrix!")
                    }
                }
            }
            None => panic!("Could not find the appropriate lattice character or bravais symbol!"),
        }

        return reduced_structure;
    }

    /// Obtains the crystallographic standard conventional structure
    /// TODO: Deal with I-centered monoclinc
    pub fn get_standard_conventional_structure(&self, structure: &Structure) -> Structure {
        let reducer = Reducer {
            dtol: self.dtol,
            atol: self.atol,
        };

        let prim_structure = reducer.find_primitive_cell(structure);
        let mut reduced_structure = reducer.niggli_reduce(&prim_structure, &1e-5);

        let lattice_character = self.get_lattice_character(&reduced_structure, &1e-5);
        let trans_mat = LATTICE_CHAR_TO_CONV_TRANS.get(&lattice_character);

        match trans_mat {
            Some(mat) => {
                let float_mat = Matrix3::from_iterator(mat.iter().map(|&x| x as f32));
                reduced_structure.apply_transformation(&float_mat.transpose(), &self.dtol);
            }
            None => panic!("Could not find the appropriate conventional transformation matrix!"),
        }

        return reduced_structure;
    }

    /// Classify the reduced structure accoring to one of the 44 lattice characters.
    fn get_lattice_character(&self, reduced_structure: &Structure, tol: &f32) -> u8 {
        let epsilon = tol * reduced_structure.volume().powf(1.0 / 3.0);

        let a_vec = Vector3::from(reduced_structure.lattice.column(0));
        let b_vec = Vector3::from(reduced_structure.lattice.column(1));
        let c_vec = Vector3::from(reduced_structure.lattice.column(2));

        let a = a_vec.magnitude_squared();
        let b = b_vec.magnitude_squared();
        let c = c_vec.magnitude_squared();

        let d = b_vec.dot(&c_vec);
        let e = a_vec.dot(&c_vec);
        let f = a_vec.dot(&b_vec);

        // Count negatives to avoid floating point multiplication
        let (num_negative, num_zero) = num_negative_zero(&vec![d, e, f], &epsilon);

        let type_i = if (num_negative % 2 == 0) && (num_zero == 0) {
            true
        } else {
            false
        };

        let def = Vector3::new(d, e, f);

        let mut char_num = 0;

        let vecs: Vec<(Vector3<f32>, u8)>;

        if cust_eq(&a, &b, &epsilon) && cust_eq(&a, &c, &epsilon) {
            vecs = Self::get_first_set(type_i, a, b, d, e, f, epsilon);
        } else if cust_eq(&a, &b, &epsilon) {
            vecs = Self::get_second_set(type_i, a, b, d, e, f, epsilon);
        } else if cust_eq(&b, &c, &epsilon) {
            vecs = Self::get_third_set(type_i, a, b, d, e, f, epsilon);
        } else {
            vecs = Self::get_fourth_set(type_i, a, b, d, e, f);
        }

        for (vec, num) in vecs {
            if vec
                .iter()
                .enumerate()
                .all(|(i, val)| cust_eq(val, &def[i], &epsilon))
            {
                char_num = num;
                break;
            }
        }

        if char_num == 14 && Self::extra_eval_a(a, b, d, e, f, epsilon) {
            char_num += 2;
        } else if char_num == 44 && Self::extra_eval_b(a, b, d, e, f, epsilon) {
            char_num -= 1;
        }

        return if char_num != 0 {
            char_num
        } else {
            panic!("Cannot find lattice character of input structure!")
        };
    }

    fn extra_eval_a(a: f32, b: f32, d: f32, e: f32, f: f32, epsilon: f32) -> bool {
        return cust_eq(&((d + e + f).abs() * 2.0), &(a + b), &epsilon);
    }

    fn extra_eval_b(a: f32, b: f32, d: f32, e: f32, f: f32, epsilon: f32) -> bool {
        let c1 = cust_eq(&((d + e + f).abs() * 2.0), &(a + b), &epsilon);
        let c2 = cust_eq(&((2.0 * d + f).abs()), &b, &epsilon);
        return c1 && c2;
    }

    /// Generate first set of lattice character DEF-vectors corresponding to A=B=C
    fn get_first_set(
        type_i: bool,
        a: f32,
        b: f32,
        d: f32,
        e: f32,
        f: f32,
        epsilon: f32,
    ) -> Vec<(Vector3<f32>, u8)> {
        let mut vecs = Vec::new();
        if type_i {
            vecs.push((Vector3::new(a, a, a) / 2.0, 1));
            vecs.push((Vector3::new(d, d, d), 2));
        } else {
            vecs.push((Vector3::new(0.0, 0.0, 0.0), 3));
            vecs.push((Vector3::new(a, a, a) / (-3.0), 5));
            vecs.push((Vector3::new(d, d, d), 4));
            if Self::extra_eval_a(a, b, d, e, f, epsilon) {
                vecs.push((Vector3::new(d, d, f), 6));
                vecs.push((Vector3::new(d, e, e), 7));
                vecs.push((Vector3::new(d, e, f), 8));
            }
        }

        return vecs;
    }

    /// Generate second set of lattice character DEF-vectors corresponding to A=B and no conditions on C
    fn get_second_set(
        type_i: bool,
        a: f32,
        b: f32,
        d: f32,
        e: f32,
        f: f32,
        epsilon: f32,
    ) -> Vec<(Vector3<f32>, u8)> {
        let mut vecs = Vec::new();
        if type_i {
            vecs.push((Vector3::new(a, a, a) / 2.0, 9));
            vecs.push((Vector3::new(d, d, f), 10));
        } else {
            vecs.push((Vector3::new(0.0, 0.0, 0.0), 11));
            vecs.push((Vector3::new(0.0, 0.0, a) / -2.0, 12));
            vecs.push((Vector3::new(0.0, 0.0, f), 13));
            vecs.push((Vector3::new(a, a, 0.0) / -2.0, 15));
            vecs.push((Vector3::new(d, d, f), 14));
            if Self::extra_eval_a(a, b, d, e, f, epsilon) {
                vecs.push((Vector3::new(d, e, f), 17));
            }
        }

        return vecs;
    }

    /// Generate third set of lattice character DEF-vectors corresponding to B=C and no conditions on A
    fn get_third_set(
        type_i: bool,
        a: f32,
        b: f32,
        d: f32,
        e: f32,
        f: f32,
        epsilon: f32,
    ) -> Vec<(Vector3<f32>, u8)> {
        let mut vecs = Vec::new();
        if type_i {
            vecs.push((Vector3::new(a / 4.0, a / 2.0, a / 2.0), 18));
            vecs.push((Vector3::new(d, a / 2.0, a / 2.0), 19));
            vecs.push((Vector3::new(d, e, e), 20));
        } else {
            vecs.push((Vector3::new(0.0, 0.0, 0.0), 21));
            vecs.push((Vector3::new(b, 0.0, 0.0) / -2.0, 22));
            vecs.push((Vector3::new(d, 0.0, 0.0), 23));
            if Self::extra_eval_a(a, b, d, e, f, epsilon) {
                vecs.push((Vector3::new(d, a / -3.0, a / -3.0), 24));
            }
            vecs.push((Vector3::new(d, e, e), 25));
        }

        return vecs;
    }

    /// Generate fourth set of lattice character DEF-vectors corresponding to no conditions on A, B, C
    fn get_fourth_set(
        type_i: bool,
        a: f32,
        b: f32,
        d: f32,
        e: f32,
        f: f32,
    ) -> Vec<(Vector3<f32>, u8)> {
        let mut vecs = Vec::new();
        if type_i {
            vecs.push((Vector3::new(a / 4.0, a / 2.0, a / 2.0), 26));
            vecs.push((Vector3::new(d, a / 2.0, a / 2.0), 27));
            vecs.push((Vector3::new(d, a / 2.0, d * 2.0), 28));
            vecs.push((Vector3::new(d, d * 2.0, a / 2.0), 29));
            vecs.push((Vector3::new(b / 2.0, e, e * 2.0), 30));
            vecs.push((Vector3::new(d, e, f), 31));
        } else {
            vecs.push((Vector3::new(0.0, 0.0, 0.0), 32));
            vecs.push((Vector3::new(b, 0.0, 0.0) / -2.0, 40));
            vecs.push((Vector3::new(d, 0.0, 0.0), 35));
            vecs.push((Vector3::new(0.0, a, 0.0) / -2.0, 36));
            vecs.push((Vector3::new(0.0, e, 0.0), 33));
            vecs.push((Vector3::new(0.0, 0.0, a) / -2.0, 38));
            vecs.push((Vector3::new(0.0, 0.0, f), 34));
            vecs.push((Vector3::new(b, a, 0.0) / -2.0, 42));
            vecs.push((Vector3::new(b / -2.0, e, 0.0), 41));
            vecs.push((Vector3::new(d, a / -2.0, 0.0), 37));
            vecs.push((Vector3::new(d, 0.0, a / -2.0), 39));

            vecs.push((Vector3::new(d, e, f), 44));
        }

        return vecs;
    }
}
