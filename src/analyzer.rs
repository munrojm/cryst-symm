use crate::data::{SELLING_TO_BRAVAIS, SELLING_TO_CONV_TRANS, ZERO_TOL};
use crate::reduce::Reducer;
use crate::structure::Structure;
use itertools::Itertools;
use nalgebra::{Matrix3, Matrix3x4, Matrix4, Vector4, Vector6};
use std::collections::HashMap;
use std::option::Option;
use std::string::String;

#[derive(Debug, Copy, Clone)]
pub struct SymmetryAnalyzer {
    pub dtol: f32,
}

impl SymmetryAnalyzer {
    pub fn get_conventional_structure(&self, structure: &Structure) -> Option<Structure> {
        let reducer = Reducer { dtol: self.dtol };

        let prim_structure = reducer.find_primitive_cell(structure);
        let del_prim_structure = reducer.delaunay_reduce(&prim_structure);

        let mut delaunay_mat = del_prim_structure.compute_delaunay_mat();

        let _ = self.delaunay_to_bravais(&mut delaunay_mat);

        let selling_vec = Reducer::get_scalar_prods(&delaunay_mat);

        let simplified_selling = self.get_simplified_selling(&selling_vec);

        let trans_mat = SELLING_TO_CONV_TRANS.get(&simplified_selling);

        // Need to alter lattice with new delaunay mat..
        // First get transformation from old -> new delaunay then transform del_prim
        let new_lattice = Matrix3::from(delaunay_mat.fixed_slice::<3, 3>(0, 0));

        let mut new_structure = Structure::new(
            new_lattice,
            del_prim_structure.species,
            del_prim_structure.coords,
            true,
        );

        match trans_mat {
            Some(mat) => new_structure.apply_transformation(&mat, &self.dtol),
            None => return Option::None,
        }

        return Option::Some(new_structure);
    }

    fn delaunay_to_bravais(&self, delaunay_mat: &mut Matrix3x4<f32>) -> Option<&String> {
        let trans_mats = self.generate_delaunay_permutation_mats();
        let mut bravais: Option<&String> = Option::None;
        let original_delaunay = delaunay_mat.clone();
        for mat in trans_mats.iter() {
            *delaunay_mat = original_delaunay * mat;
            let selling_vec = Reducer::get_scalar_prods(&delaunay_mat);
            let simplified_selling = self.get_simplified_selling(&selling_vec);
            bravais = SELLING_TO_BRAVAIS.get(&simplified_selling);
            match bravais {
                Some(_) => break,
                None => (),
            }
        }
        return bravais;
    }

    fn get_simplified_selling(&self, &selling_vec: &Vector6<f32>) -> Vector6<usize> {
        let mut value_map: HashMap<usize, usize> = HashMap::from([(0, 0)]);
        let mut count = 1;
        let mut simplified_selling = Vector6::new(0, 0, 0, 0, 0, 0);
        for (i, j) in selling_vec.iter().enumerate() {
            let mut rounded_val = (j.abs() * 100.0) as usize; // May need to match with tolerance here later.
            if j.abs() <= ZERO_TOL {
                rounded_val = 0
            }
            let new_val: usize;
            match value_map.get(&rounded_val) {
                Some(t) => new_val = t.clone(),
                None => {
                    value_map.insert(rounded_val, count);
                    new_val = count;
                    count += 1;
                }
            }
            simplified_selling[i] = new_val;
        }
        return simplified_selling;
    }

    fn generate_delaunay_permutation_mats(&self) -> Vec<Matrix4<f32>> {
        let base_vecs: Vec<Vector4<f32>> = vec![
            Vector4::new(1.0, 0.0, 0.0, 0.0),
            Vector4::new(0.0, 1.0, 0.0, 0.0),
            Vector4::new(0.0, 0.0, 1.0, 0.0),
            Vector4::new(0.0, 0.0, 0.0, 1.0),
        ];

        let mut mats: Vec<Matrix4<f32>> = Vec::new();

        for vec in base_vecs.into_iter().permutations(4) {
            let mat = Matrix4::from_columns(&vec[..]);
            mats.push(mat);
        }

        return mats;
    }
}
