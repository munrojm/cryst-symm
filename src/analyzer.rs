use crate::data::core::{
    BravaisType, CENTERING_TO_PRIM_TRANS, INT_CONVERSION_MULT, LATTICE_CHAR_TO_BRAVAIS,
    LATTICE_CHAR_TO_CONV_TRANS,
};
use crate::data::pointgroup::BRAVAIS_TO_HOLOHEDRY_NUM;
use crate::pointgroup::PointGroup;
use crate::reduce::Reducer;
use crate::spacegroup::SpaceGroup;
use crate::structure::Structure;
use crate::symmop::SymmOp;
use crate::utils::{cust_eq, normalize_frac_vectors, num_negative_zero};
use nalgebra::{Matrix3, Vector3};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct SymmetryAnalyzer {
    pub dtol: f64,
    pub atol: f64,
}

impl SymmetryAnalyzer {
    pub fn get_space_group_operations(&self, structure: &Structure) -> SpaceGroup {
        // Get conventional structure, classify, and get space group generators
        let reducer = Reducer {
            dtol: self.dtol,
            atol: self.atol,
        };

        let prim_structure = reducer.find_primitive_cell(structure);

        let mut reduced_structure = reducer.niggli_reduce(&prim_structure, &1e-5);

        let (conv_trans_mat, lattice_character) =
            self.conventional_transformation_matrix(&reduced_structure);

        reduced_structure.apply_transformation(&conv_trans_mat.transpose(), &self.dtol); // Standard conventional structure

        let (prim_trans_mat, bravais_symbol) =
            self.primitive_transformation_matrix(&lattice_character);

        reduced_structure.apply_transformation(&prim_trans_mat, &self.dtol); // Standard primitive structure

        let sg_generators = self.get_primitive_space_group_ops(&reduced_structure, bravais_symbol);

        // Generate full coset representatives of space group with primitive generators and centering
        let frac_tols = Vector3::from_iterator(
            reduced_structure
                .lattice
                .column_iter()
                .map(|col| self.dtol / col.magnitude()),
        );

        let mut sg = SpaceGroup::from_generators(&sg_generators, &frac_tols);

        // Put operations back into the basis of the input structure
        // TODO: Need to fix transforming to the larger cell.
        let niggli_trans_mat = prim_structure.get_transformation_matrix(structure);

        let total_trans_mat = niggli_trans_mat * conv_trans_mat * prim_trans_mat;

        sg.apply_transformation(&total_trans_mat);

        return sg;
    }

    fn get_primitive_space_group_ops(
        &self,
        prim_structure: &Structure,
        bravais_symbol: &BravaisType,
    ) -> Vec<SymmOp> {
        let frac_tols = Vector3::from_iterator(
            prim_structure
                .lattice
                .column_iter()
                .map(|col| self.dtol / col.magnitude()),
        );

        // Get holohedral point group in primitive basis
        let holohedry_num = BRAVAIS_TO_HOLOHEDRY_NUM.get(bravais_symbol).unwrap();
        let mut holohedry_pg = PointGroup::from_number(holohedry_num);

        let int_trans = CENTERING_TO_PRIM_TRANS
            .get(&bravais_symbol.centering())
            .unwrap();

        let float_trans = int_trans.cast::<f64>();

        let prim_transformation: Matrix3<f64> =
            Matrix3::from(float_trans.fixed_slice::<3, 3>(0, 0)) / float_trans.m44;

        holohedry_pg.apply_transformation(&prim_transformation);

        // 1. Apply holohedral point group operation
        // 2. Generate all candidate translation vectors to other sites
        // 3. If the occurance number of a particular vector is equal to number of sites,
        //    choose it as canonical translation for the operation.
        // NOTE: These translations include both instrinsic translation and origin shift.
        let mut op_trans_vecs: HashMap<Matrix3<i8>, Vec<Vector3<f64>>> = HashMap::new();

        let (_, ele_inds, ele_counts) = prim_structure.get_min_element();

        for op in holohedry_pg.operations {
            let float_op = Matrix3::from_iterator(op.iter().map(|&x| x as f64));
            let mut translation_vectors: Vec<Vector3<f64>> = Vec::new();

            for (ele, ind) in ele_inds.iter() {
                let start = ind.to_owned();
                let end = start + ele_counts.get(ele).unwrap().to_owned();

                for coord_ind_a in start..end {
                    let rotated_coord = float_op * prim_structure.frac_coords[coord_ind_a as usize];

                    for coord_ind_b in start..end {
                        let mut translation_vec =
                            vec![prim_structure.frac_coords[coord_ind_b as usize] - rotated_coord];

                        normalize_frac_vectors(&mut translation_vec, &frac_tols);

                        translation_vectors.push(translation_vec[0]);
                    }
                }
            }

            op_trans_vecs.insert(op, translation_vectors);
        }

        // Sort through candidate vectors and get counts
        let mut sg_operations: Vec<SymmOp> = Vec::new();
        let mut matched;

        for (op, translation_vectors) in op_trans_vecs.iter() {
            let mut translation_counts: HashMap<Vector3<i64>, u64> = HashMap::new();

            for translation in translation_vectors.iter() {
                let int_translation: Vector3<i64> = Vector3::from_iterator(
                    translation
                        .iter()
                        .map(|&x| (x * INT_CONVERSION_MULT) as i64),
                );

                matched = false;

                for (count_vector, count) in translation_counts.iter_mut() {
                    let float_count_translation: Vector3<f64> = Vector3::from_iterator(
                        count_vector.iter().map(|&x| x as f64 / INT_CONVERSION_MULT),
                    );

                    let mut diff = vec![translation - float_count_translation];
                    normalize_frac_vectors(&mut diff, &frac_tols);

                    let cart_coord_delta =
                        Structure::get_cart_coords(&prim_structure.lattice, &diff);

                    if cart_coord_delta[0].magnitude().abs() <= (2.0 * self.dtol) {
                        *count += 1;
                        matched = true;
                        break;
                    }
                }

                if !matched {
                    translation_counts.insert(int_translation, 1);
                }
            }

            // Check for counts that are equal to the number of sites and choose final translation
            // vector for the given rotation operation.
            for (int_translation_vector, count) in translation_counts.into_iter() {
                if count as usize == prim_structure.num_sites() {
                    let float_translation: Vector3<f64> = Vector3::from_iterator(
                        int_translation_vector
                            .iter()
                            .map(|&x| x as f64 / INT_CONVERSION_MULT),
                    );

                    sg_operations.push(SymmOp {
                        rotation: op.clone(),
                        translation: float_translation,
                    });
                }
            }
        }

        return sg_operations;
    }

    /// Obtains the crystallographic primitive crystal structure
    pub fn get_standard_primitive_structure(&self, structure: &Structure) -> Structure {
        let reducer = Reducer {
            dtol: self.dtol,
            atol: self.atol,
        };

        let prim_structure = reducer.find_primitive_cell(structure);

        let mut reduced_structure = reducer.niggli_reduce(&prim_structure, &1e-5);

        let (trans_mat, lattice_character) =
            self.conventional_transformation_matrix(&reduced_structure);

        reduced_structure.apply_transformation(&trans_mat.transpose(), &self.dtol);

        let (prim_trans_mat, _) = self.primitive_transformation_matrix(&lattice_character);

        reduced_structure.apply_transformation(&prim_trans_mat, &self.dtol);

        return reduced_structure;
    }

    fn primitive_transformation_matrix(
        &self,
        lattice_character: &u8,
    ) -> (Matrix3<f64>, &BravaisType) {
        let bravais_symbol = LATTICE_CHAR_TO_BRAVAIS
            .get(&lattice_character)
            .expect("Could not find the appropriate lattice character or bravais symbol!");

        let prim_trans_mat_int = CENTERING_TO_PRIM_TRANS
            .get(&bravais_symbol.centering())
            .expect("Could not find the appropriate primitive transformation matrix!");

        let float_mat = prim_trans_mat_int.cast::<f64>();

        let prim_trans_mat: Matrix3<f64> =
            Matrix3::from(float_mat.fixed_slice::<3, 3>(0, 0)) / float_mat.m44;

        return (prim_trans_mat, bravais_symbol);
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

        let (trans_mat, _) = self.conventional_transformation_matrix(&reduced_structure);

        reduced_structure.apply_transformation(&trans_mat.transpose(), &self.dtol);

        return reduced_structure;
    }

    fn conventional_transformation_matrix(
        &self,
        reduced_structure: &Structure,
    ) -> (Matrix3<f64>, u8) {
        let lattice_character = self.get_lattice_character(&reduced_structure, &1e-5);

        let trans_mat = LATTICE_CHAR_TO_CONV_TRANS
            .get(&lattice_character)
            .expect("Could not find the appropriate conventional transformation matrix!");

        let float_mat = trans_mat.cast::<f64>();

        return (float_mat, lattice_character);
    }

    /// Classify the reduced structure accoring to one of the 44 lattice characters.
    fn get_lattice_character(&self, reduced_structure: &Structure, tol: &f64) -> u8 {
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

        let vecs: Vec<(Vector3<f64>, u8)>;

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

    fn extra_eval_a(a: f64, b: f64, d: f64, e: f64, f: f64, epsilon: f64) -> bool {
        return cust_eq(&((d + e + f).abs() * 2.0), &(a + b), &epsilon);
    }

    fn extra_eval_b(a: f64, b: f64, d: f64, e: f64, f: f64, epsilon: f64) -> bool {
        let c1 = cust_eq(&((d + e + f).abs() * 2.0), &(a + b), &epsilon);
        let c2 = cust_eq(&((2.0 * d + f).abs()), &b, &epsilon);
        return c1 && c2;
    }

    /// Generate first set of lattice character DEF-vectors corresponding to A=B=C
    fn get_first_set(
        type_i: bool,
        a: f64,
        b: f64,
        d: f64,
        e: f64,
        f: f64,
        epsilon: f64,
    ) -> Vec<(Vector3<f64>, u8)> {
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
        a: f64,
        b: f64,
        d: f64,
        e: f64,
        f: f64,
        epsilon: f64,
    ) -> Vec<(Vector3<f64>, u8)> {
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
        a: f64,
        b: f64,
        d: f64,
        e: f64,
        f: f64,
        epsilon: f64,
    ) -> Vec<(Vector3<f64>, u8)> {
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
        a: f64,
        b: f64,
        d: f64,
        e: f64,
        f: f64,
    ) -> Vec<(Vector3<f64>, u8)> {
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
