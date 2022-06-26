use crate::structure::Structure;
use nalgebra::{Matrix3, Vector3};
use std::collections::HashMap;
use std::string::String;

#[derive(Debug)]
pub struct Reducer {
    pub dtol: f32,
}

impl Reducer {
    pub fn find_primitive_cell(self, structure: &Structure) -> Structure {
        let mut temp_structure = structure.clone();
        let frac_tols = Vector3::from_iterator(
            temp_structure
                .lattice
                .column_iter()
                .map(|col| self.dtol / col.magnitude()),
        );

        //
        // 1.) Find atom type with the fewest sites
        //
        let mut type_count: HashMap<&String, u8> = HashMap::new();
        let mut ele_inds: HashMap<&String, u8> = HashMap::new();

        for (ele_ind, ele) in structure.species.iter().enumerate() {
            *type_count.entry(ele).or_insert(0) += 1;
            ele_inds.entry(ele).or_insert(ele_ind as u8);
        }

        let min_ele = type_count
            .iter()
            .min_by_key(|entry| entry.1)
            .unwrap()
            .0
            .clone();
        //
        // 2.) Shift origin to first atom of min_ele type and normalize coords
        //
        let new_origin = temp_structure.frac_coords[*ele_inds.get(&min_ele).unwrap() as usize];
        temp_structure.set_origin(new_origin);
        temp_structure.normalize_coords(self.dtol);

        //
        // 3.) Get all potential new unit translation lattice vectors
        //
        let mut candidate_vecs: Vec<Vector3<f32>> = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];

        let base_site = temp_structure.frac_coords[*ele_inds.get(min_ele).unwrap() as usize];

        for (coord_ind, coord) in temp_structure.frac_coords.iter().enumerate() {
            if (&temp_structure.species[coord_ind] == min_ele)
                && (coord_ind != *ele_inds.get(min_ele).unwrap() as usize)
            {
                candidate_vecs.push(coord - &base_site);
            }
        }

        let mut translation_vecs: Vec<Vector3<f32>> = Vec::new();
        let mut all_matched: bool;
        let mut matched: bool;

        for vec in candidate_vecs.iter() {
            all_matched = true;
            for (coord, specie) in temp_structure
                .frac_coords
                .iter()
                .zip(temp_structure.species.iter())
            {
                let tr_coord = coord + vec;
                matched = false;
                for (coord_p, specie_p) in temp_structure
                    .frac_coords
                    .iter()
                    .zip(temp_structure.species.iter())
                {
                    if specie == specie_p {
                        let mut coord_delta = tr_coord - coord_p;
                        Self::normalize_frac_vector(&mut coord_delta, &(2.0 * &frac_tols));

                        let cart_coord_delta =
                            Structure::get_cart_coords(&temp_structure.lattice, &vec![coord_delta]);

                        if cart_coord_delta.first().unwrap().magnitude().abs() <= (2.0 * self.dtol)
                        {
                            matched = true;
                            break;
                        }
                    }
                }
                if !matched {
                    all_matched = false;
                    break;
                }
            }
            if all_matched {
                translation_vecs.push(*vec);
            }
        }

        //
        // 3.) Sort potential vectors from smallest to largest and find three which span the lattice
        //

        let mut cart_translation_vecs: Vec<Vector3<f32>> =
            Structure::get_cart_coords(&temp_structure.lattice, &translation_vecs);

        cart_translation_vecs.sort_by(|a, b| a.magnitude().partial_cmp(&b.magnitude()).unwrap());

        let first = cart_translation_vecs[0];
        let mut second = cart_translation_vecs[1];
        let mut third = cart_translation_vecs[2];
        let cross_vec = first.cross(&second);
        let mut second_ind = 1;

        for (vec_num, vec) in cart_translation_vecs[1..].iter().enumerate() {
            if cross_vec.magnitude().abs() <= self.dtol {
                second = vec.clone();
                second_ind = vec_num;
                break;
            }
        }

        for vec in cart_translation_vecs[(second_ind + 1)..].iter() {
            if !((cross_vec.dot(vec)).abs() <= self.dtol) {
                third = vec.clone();
                break;
            }
        }
        //
        // 5.) Create a new structure object with new lattice and coords.
        //     TODO: Average positions which fold into one another in the new cell.
        //

        let new_lattice = Matrix3::from_columns(&[first, second, third]);

        let temp_frac_coords = Structure::get_frac_coords(&new_lattice, &temp_structure.coords);

        let mut new_frac_coords: Vec<Vector3<f32>> = Vec::new();
        let mut new_species: Vec<String> = Vec::new();

        for (coord, specie) in temp_frac_coords.iter().zip(&structure.species) {
            let mut eq = false;
            for new_coord in new_frac_coords.iter() {
                let mut coord_diff = coord - new_coord;
                Self::normalize_frac_vector(&mut coord_diff, &frac_tols);
                let cart_coord_delta =
                    Structure::get_cart_coords(&temp_structure.lattice, &vec![coord_diff]);

                if cart_coord_delta[0].magnitude().abs() <= (self.dtol * 2.0) {
                    eq = true;
                    break;
                }
            }
            if !eq {
                let mut norm_coord = coord.clone();
                Self::normalize_frac_vector(&mut norm_coord, &frac_tols);
                new_frac_coords.push(norm_coord);
                new_species.push(specie.clone());
            }
        }

        let mut prim_structure = Structure::new(new_lattice, new_species, new_frac_coords, false);
        prim_structure.normalize_coords(self.dtol);
        return prim_structure;
    }

    fn normalize_frac_vector(vec: &mut Vector3<f32>, frac_tols: &Vector3<f32>) {
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
