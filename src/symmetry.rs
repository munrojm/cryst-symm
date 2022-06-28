use crate::data::SELLING_TO_BRAVAIS;
use crate::reduce::Reducer;
use crate::structure::Structure;
use nalgebra::{Matrix3, Matrix4, Vector6};
use std::collections::HashMap;
use std::option::Option;
use std::string::String;

#[derive(Debug, Copy, Clone)]
pub struct SymmetryAnalyzer {
    pub dtol: f32,
}

impl SymmetryAnalyzer {
    pub fn get_crystal_system(
        self,
        structure: &Structure,
    ) -> Option<&(String, String, Matrix3<isize>)> {
        let reducer = Reducer { dtol: self.dtol };

        let prim_structure = reducer.find_primitive_cell(structure);
        let del_prim_structure = reducer.delaunay_reduce(&prim_structure);

        let trans_mats: Vec<Matrix4<f32>> = vec![
            Matrix4::identity(),
            Matrix4::new(
                0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ),
            Matrix4::new(
                0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ),
            Matrix4::new(
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ),
        ];
        let delaunay_mat = Reducer::compute_delaunay_mat(&del_prim_structure.lattice);

        let mut lookup_data: Option<&(String, String, Matrix3<isize>)> = Option::None;

        for mat in trans_mats.iter() {
            let selling_vec = Reducer::get_scalar_prods(&(delaunay_mat * mat));
            let simplified_selling = self.get_simplified_selling(&selling_vec);
            println!("{:?}", (selling_vec, simplified_selling));
            lookup_data = SELLING_TO_BRAVAIS.get(&simplified_selling);
            match lookup_data {
                Some(_) => break,
                None => (),
            }
        }

        return lookup_data;
    }

    fn get_simplified_selling(self, &selling_vec: &Vector6<f32>) -> Vector6<usize> {
        let mut value_map: HashMap<usize, usize> = HashMap::from([(0, 0)]);
        let mut count = 1;
        let mut simplified_selling = Vector6::new(0, 0, 0, 0, 0, 0);
        for (i, j) in selling_vec.iter().enumerate() {
            let mut rounded_val = (j.abs() / self.dtol) as usize; // Will need to match with tolerance here later.
            if j.abs() <= self.dtol / 100.0 {
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
}
