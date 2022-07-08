use crate::structure::Structure;
use nalgebra::{Matrix3, Matrix3x4, Vector3, Vector6};
use std::collections::HashMap;
use std::string::String;

#[derive(Debug, Copy, Clone)]
pub struct Reducer {
    pub dtol: f32,
}

impl Reducer {
    pub fn find_primitive_cell(&self, structure: &Structure) -> Structure {
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

        let frac_translation_vecs = self.get_shortest_translation_vecs(&mut cart_translation_vecs);

        //
        // 5.) Add combinations of shortest vectors
        //

        let smallest_frac_vecs =
            Structure::get_frac_coords(&temp_structure.lattice, &frac_translation_vecs);

        let frac_matrix = Matrix3::from_columns(&smallest_frac_vecs);

        let base_comb_vecs: Vec<Vector3<f32>> = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ];

        let mut combination_vecs: Vec<Vector3<f32>> = Vec::new();

        for comb_vec in base_comb_vecs {
            let transformed_vec = frac_matrix * comb_vec;
            combination_vecs.push(transformed_vec);
            combination_vecs.push(-1.0 * transformed_vec);

            for i in 0..3 {
                let mut new_vec = transformed_vec.clone();
                new_vec[i] = -1.0 * new_vec[i];
                combination_vecs.push(new_vec);
            }
        }

        let mut cart_translation_vecs: Vec<Vector3<f32>> =
            Structure::get_cart_coords(&temp_structure.lattice, &combination_vecs);

        let final_cart_vecs = self.get_shortest_translation_vecs(&mut cart_translation_vecs);

        //
        // 6.) Create a new structure object with new lattice and coords.
        //     TODO: Average positions which fold into one another in the new cell.
        //

        let new_lattice = Matrix3::from_columns(&final_cart_vecs);

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

    pub fn delaunay_reduce(&self, structure: &Structure) -> Structure {
        let mut new_structure = structure.clone();
        let pair_map: HashMap<usize, (u8, u8)> = HashMap::from([
            (0, (0, 1)),
            (1, (0, 2)),
            (2, (0, 3)),
            (3, (1, 2)),
            (4, (1, 3)),
            (5, (2, 3)),
        ]);

        let mut delaunay_mat = structure.compute_delaunay_mat();

        let mut scalar_prods = Self::get_scalar_prods(&delaunay_mat);

        let mut count = 0;

        while scalar_prods.iter().any(|val| val > &0.0001) && count <= 100 {
            let mut pos_pair = &(0, 0);

            for (i, entry) in scalar_prods.iter().enumerate() {
                if entry > &0.0001 {
                    pos_pair = pair_map.get(&i).unwrap();
                }
            }

            for i in 0..4 {
                if !(pos_pair.0 == i || pos_pair.1 == i) {
                    let new_col = &delaunay_mat.column(pos_pair.0 as usize)
                        + &delaunay_mat.column(i as usize);
                    delaunay_mat.set_column(i as usize, &new_col);
                };
            }

            let neg_column = &delaunay_mat.column(pos_pair.0 as usize) * (-1.0);

            delaunay_mat.set_column(pos_pair.0 as usize, &neg_column);
            scalar_prods = Self::get_scalar_prods(&delaunay_mat);
            count += 1;
        }

        if count == 100 {
            panic!("Reached max number of iterations without reduction!")
        }

        let new_lattice = Matrix3::from(delaunay_mat.fixed_columns::<3>(0));
        let new_frac_coords = Structure::get_frac_coords(&new_lattice, &new_structure.coords);

        new_structure = Structure::new(new_lattice, new_structure.species, new_frac_coords, false);
        new_structure.normalize_coords(self.dtol);

        return new_structure;
    }

    fn get_shortest_translation_vecs(
        &self,
        cart_vecs: &mut Vec<Vector3<f32>>,
    ) -> Vec<Vector3<f32>> {
        cart_vecs.sort_by(|a, b| a.magnitude().partial_cmp(&b.magnitude()).unwrap());

        let first = cart_vecs[0];
        let mut cross_vec: Vector3<f32> = Vector3::new(0.0, 0.0, 0.0);
        let mut second_ind = 1;
        let mut third_ind = 2;

        for (vec_num, cart_vec) in cart_vecs[1..].iter().enumerate() {
            cross_vec = first.cross(cart_vec);

            if cross_vec.magnitude().abs() > self.dtol {
                second_ind = vec_num + 1;
                break;
            }
        }

        for (vec_num, cart_vec) in cart_vecs[(second_ind + 1)..].iter().enumerate() {
            if !((cross_vec.dot(cart_vec)).abs() <= self.dtol) {
                third_ind = vec_num + second_ind + 1;
                break;
            }
        }

        let final_vecs = vec![cart_vecs[0], cart_vecs[second_ind], cart_vecs[third_ind]];

        return final_vecs;
    }

    pub fn get_scalar_prods(delaunay_mat: &Matrix3x4<f32>) -> Vector6<f32> {
        let mut scalar_prods: Vector6<f32> = Vector6::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let mut ind = 0;
        for i in 0..3 {
            for j in (i + 1)..4 {
                scalar_prods[ind] = delaunay_mat.column(i).dot(&delaunay_mat.column(j));
                ind += 1;
            }
        }
        return scalar_prods;
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
