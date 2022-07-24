use crate::data::{LATTICE_CHAR_TO_CONV_TRANS, ZERO_TOL};
use crate::reduce::Reducer;
use crate::structure::Structure;
use crate::utils::calculate_dot_uncertainty;
use nalgebra::{Matrix3, Vector3};

#[derive(Debug, Copy, Clone)]
pub struct SymmetryAnalyzer {
    pub dtol: f32,
    pub atol: f32,
}

impl SymmetryAnalyzer {
    /// Obtains the conventional structure
    /// TODO: Deal with I-centered monoclinc
    pub fn get_standard_conventional_structure(&self, structure: &Structure) -> Structure {
        let reducer = Reducer {
            dtol: self.dtol,
            atol: self.atol,
        };

        let prim_structure = reducer.find_primitive_cell(structure);
        let mut reduced_structure = reducer.niggli_reduce(&prim_structure, &1e-5);

        let lattice_character = self.get_lattice_character(&reduced_structure);
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
    fn get_lattice_character(&self, reduced_structure: &Structure) -> u8 {
        let a_vec = Vector3::from(reduced_structure.lattice.column(0));
        let b_vec = Vector3::from(reduced_structure.lattice.column(1));
        let c_vec = Vector3::from(reduced_structure.lattice.column(2));

        let a = a_vec.magnitude_squared();
        let da = 2.0_f32.sqrt() * a_vec.magnitude() * self.dtol;

        let b = b_vec.magnitude_squared();
        let db = 2.0_f32.sqrt() * b_vec.magnitude() * self.dtol;

        let c = c_vec.magnitude_squared();
        let dc = 2.0_f32.sqrt() * c_vec.magnitude() * self.dtol;

        let d = b_vec.dot(&c_vec);
        let dd = calculate_dot_uncertainty(&b_vec, &c_vec, &self.dtol, &self.dtol, &self.atol);

        let e = a_vec.dot(&c_vec);
        let de = calculate_dot_uncertainty(&a_vec, &c_vec, &self.dtol, &self.dtol, &self.atol);

        let f = a_vec.dot(&b_vec);
        let df = calculate_dot_uncertainty(&a_vec, &b_vec, &self.dtol, &self.dtol, &self.atol);

        let type_i = if (d * e * f) > ZERO_TOL { true } else { false };

        let def = Vector3::new(d, e, f);
        let d_def = Vector3::new(dd, de, df);

        let mut char_num = 0;

        let vecs: Vec<(Vector3<f32>, Vector3<f32>, u8)>;

        if (a - b).abs() < (da + db) && (a - c).abs() < (da + dc) {
            vecs = Self::get_first_set(type_i, a, da, b, db, d, dd, e, de, f, df);
        } else if (a - b).abs() < (da + db) {
            vecs = Self::get_second_set(type_i, a, da, b, db, d, dd, e, de, f, df);
        } else if (b - c).abs() < (db + dc) {
            vecs = Self::get_third_set(type_i, a, da, b, db, d, dd, e, de, f, df);
        } else {
            vecs = Self::get_fourth_set(type_i, a, da, b, db, d, dd, e, de, f, df);
        }

        for (vec, d_vec, num) in vecs {
            if vec
                .iter()
                .enumerate()
                .all(|(i, val)| (val - def[i]).abs() < (d_def[i] + d_vec[i]))
            {
                char_num = num;
                break;
            }
        }

        if char_num == 14 && Self::extra_eval_a(a, da, b, db, d, dd, e, de, f, df) {
            char_num += 2;
        } else if char_num == 44 && Self::extra_eval_b(a, da, b, db, d, dd, e, de, f, df) {
            char_num -= 1;
        }

        return if char_num != 0 {
            char_num
        } else {
            panic!("Cannot find lattice character of input structure!")
        };
    }

    fn extra_eval_a(
        a: f32,
        da: f32,
        b: f32,
        db: f32,
        d: f32,
        dd: f32,
        e: f32,
        de: f32,
        f: f32,
        df: f32,
    ) -> bool {
        return (((d + e + f).abs() * 2.0) - (a + b)).abs() < (((dd + de + df) * 2.0) + da + db);
    }

    fn extra_eval_b(
        a: f32,
        da: f32,
        b: f32,
        db: f32,
        d: f32,
        dd: f32,
        e: f32,
        de: f32,
        f: f32,
        df: f32,
    ) -> bool {
        let c1 = (((d + e + f).abs() * 2.0) - (a + b)).abs() < (((dd + de + df) * 2.0) + da + db);
        let c2 = ((2.0 * d + f).abs() - b).abs() < (2.0 * dd + df + db);
        return c1 && c2;
    }

    /// Generate first set of lattice character DEF-vectors corresponding to A=B=C
    fn get_first_set(
        type_i: bool,
        a: f32,
        da: f32,
        b: f32,
        db: f32,
        d: f32,
        dd: f32,
        e: f32,
        de: f32,
        f: f32,
        df: f32,
    ) -> Vec<(Vector3<f32>, Vector3<f32>, u8)> {
        let mut vecs = Vec::new();
        if type_i {
            vecs.push((
                Vector3::new(a, a, a) / 2.0,
                Vector3::new(da, da, da) / 2.0,
                1,
            ));
            vecs.push((Vector3::new(d, d, d), Vector3::new(dd, dd, dd), 2));
        } else {
            vecs.push((Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0), 3));
            vecs.push((
                Vector3::new(a, a, a) / (-3.0),
                Vector3::new(da, da, da) / 3.0,
                5,
            ));
            vecs.push((Vector3::new(d, d, d), Vector3::new(dd, dd, dd), 4));
            if Self::extra_eval_a(a, da, b, db, d, dd, e, de, f, df) {
                vecs.push((Vector3::new(d, d, f), Vector3::new(dd, dd, df), 6));
                vecs.push((Vector3::new(d, e, e), Vector3::new(dd, de, de), 7));
                vecs.push((Vector3::new(d, e, f), Vector3::new(dd, de, df), 8));
            }
        }

        return vecs;
    }

    /// Generate second set of lattice character DEF-vectors corresponding to A=B and no conditions on C
    fn get_second_set(
        type_i: bool,
        a: f32,
        da: f32,
        b: f32,
        db: f32,
        d: f32,
        dd: f32,
        e: f32,
        de: f32,
        f: f32,
        df: f32,
    ) -> Vec<(Vector3<f32>, Vector3<f32>, u8)> {
        let mut vecs = Vec::new();
        if type_i {
            vecs.push((
                Vector3::new(a, a, a) / 2.0,
                Vector3::new(da, da, da) / 2.0,
                9,
            ));
            vecs.push((Vector3::new(d, d, f), Vector3::new(dd, dd, df), 10));
        } else {
            vecs.push((Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0), 11));
            vecs.push((
                Vector3::new(0.0, 0.0, a) / -2.0,
                Vector3::new(0.0, 0.0, da) / 2.0,
                12,
            ));
            vecs.push((Vector3::new(0.0, 0.0, f), Vector3::new(0.0, 0.0, df), 13));
            vecs.push((
                Vector3::new(a, a, 0.0) / -2.0,
                Vector3::new(da, da, 0.0) / 2.0,
                15,
            ));
            vecs.push((Vector3::new(d, d, f), Vector3::new(dd, dd, df), 14));
            if (((d + e + f).abs() * 2.0) - (a + b)).abs() < (((dd + de + df) * 2.0) + da + db) {
                vecs.push((Vector3::new(d, e, f), Vector3::new(dd, de, df), 17));
            }
        }

        return vecs;
    }

    /// Generate third set of lattice character DEF-vectors corresponding to B=C and no conditions on A
    fn get_third_set(
        type_i: bool,
        a: f32,
        da: f32,
        b: f32,
        db: f32,
        d: f32,
        dd: f32,
        e: f32,
        de: f32,
        f: f32,
        df: f32,
    ) -> Vec<(Vector3<f32>, Vector3<f32>, u8)> {
        let mut vecs = Vec::new();
        if type_i {
            vecs.push((
                Vector3::new(a / 4.0, a / 2.0, a / 2.0),
                Vector3::new(da / 4.0, da / 2.0, da / 2.0),
                18,
            ));
            vecs.push((
                Vector3::new(d, a / 2.0, a / 2.0),
                Vector3::new(dd, da / 2.0, da / 2.0),
                19,
            ));
            vecs.push((Vector3::new(d, e, e), Vector3::new(dd, de, de), 20));
        } else {
            vecs.push((Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0), 21));
            vecs.push((
                Vector3::new(b, 0.0, 0.0) / -2.0,
                Vector3::new(db, 0.0, 0.0) / 2.0,
                22,
            ));
            vecs.push((Vector3::new(d, 0.0, 0.0), Vector3::new(dd, 0.0, 0.0), 23));
            if (((d + e + f).abs() * 2.0) - (a + b)).abs() < (((dd + de + df) * 2.0) + da + db) {
                vecs.push((
                    Vector3::new(d, a / -3.0, a / -3.0),
                    Vector3::new(dd, da / 3.0, da / 3.0),
                    24,
                ));
            }
            vecs.push((Vector3::new(d, e, e), Vector3::new(dd, de, de), 25));
        }

        return vecs;
    }

    /// Generate fourth set of lattice character DEF-vectors corresponding to no conditions on A, B, C
    fn get_fourth_set(
        type_i: bool,
        a: f32,
        da: f32,
        b: f32,
        db: f32,
        d: f32,
        dd: f32,
        e: f32,
        de: f32,
        f: f32,
        df: f32,
    ) -> Vec<(Vector3<f32>, Vector3<f32>, u8)> {
        let mut vecs = Vec::new();
        if type_i {
            vecs.push((
                Vector3::new(a / 4.0, a / 2.0, a / 2.0),
                Vector3::new(da / 4.0, da / 2.0, da / 2.0),
                26,
            ));
            vecs.push((
                Vector3::new(d, a / 2.0, a / 2.0),
                Vector3::new(dd, da / 2.0, da / 2.0),
                27,
            ));
            vecs.push((
                Vector3::new(d, a / 2.0, d * 2.0),
                Vector3::new(dd, da / 2.0, dd * 2.0),
                28,
            ));
            vecs.push((
                Vector3::new(d, d * 2.0, a / 2.0),
                Vector3::new(dd, dd * 2.0, da / 2.0),
                29,
            ));
            vecs.push((
                Vector3::new(b / 2.0, e, e * 2.0),
                Vector3::new(db / 2.0, de, de * 2.0),
                30,
            ));
            vecs.push((Vector3::new(d, e, f), Vector3::new(dd, de, df), 31));
        } else {
            vecs.push((Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0), 32));
            vecs.push((
                Vector3::new(b, 0.0, 0.0) / -2.0,
                Vector3::new(db, 0.0, 0.0) / 2.0,
                40,
            ));
            vecs.push((Vector3::new(d, 0.0, 0.0), Vector3::new(dd, 0.0, 0.0), 35));
            vecs.push((
                Vector3::new(0.0, a, 0.0) / -2.0,
                Vector3::new(0.0, da, 0.0) / 2.0,
                36,
            ));
            vecs.push((Vector3::new(0.0, e, 0.0), Vector3::new(0.0, de, 0.0), 33));
            vecs.push((
                Vector3::new(0.0, 0.0, a) / -2.0,
                Vector3::new(0.0, 0.0, da) / 2.0,
                38,
            ));
            vecs.push((Vector3::new(0.0, 0.0, f), Vector3::new(0.0, 0.0, df), 34));
            vecs.push((
                Vector3::new(b, a, 0.0) / -2.0,
                Vector3::new(db, da, 0.0) / 2.0,
                42,
            ));
            vecs.push((
                Vector3::new(b / -2.0, e, 0.0),
                Vector3::new(db / 2.0, de, 0.0),
                41,
            ));
            vecs.push((
                Vector3::new(d, a / -2.0, 0.0),
                Vector3::new(dd, da / 2.0, 0.0),
                37,
            ));
            vecs.push((
                Vector3::new(d, 0.0, a / -2.0),
                Vector3::new(dd, 0.0, da / 2.0),
                39,
            ));

            vecs.push((Vector3::new(d, e, f), Vector3::new(dd, de, df), 44));
        }

        return vecs;
    }
}
