use crate::data::ZERO_TOL;
use crate::structure::Structure;
use crate::utils::normalize_frac_vectors;
use nalgebra::{Matrix3, Matrix3x4, Vector3, Vector6};
use std::collections::HashMap;
use std::iter::FromIterator;
use std::string::String;

#[derive(Debug, Clone)]
pub struct Reducer {
    pub dtol: f32,
    pub atol: f32,
}

impl Reducer {
    /// Produces a primitive structure using the following steps:
    ///
    /// 1. Find atom type with fewest sites.
    /// 2. Get all new potential unit translation vectors.
    /// 3. Select three shortest vectors which span the lattice.
    /// 4. Create new structure and fold atomic sites into new smaller cell.
    ///
    pub fn find_primitive_cell(&self, structure: &Structure) -> Structure {
        let temp_structure = structure.clone();
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

        // Sorting ensures determinism when getting min
        let mut sorted_pairs: Vec<(&String, u8)> = Vec::from_iter(type_count.into_iter());

        sorted_pairs.sort_by(|a, b| a.0.cmp(b.0));

        let min_ele = &sorted_pairs
            .iter()
            .min_by_key(|entry| entry.1)
            .unwrap()
            .0
            .clone();
        //
        // 2.) Shift origin to first atom of min_ele type and normalize coords (not needed)
        //
        // let new_origin = temp_structure.frac_coords[*ele_inds.get(&min_ele).unwrap() as usize];
        // temp_structure.set_origin(new_origin);
        // temp_structure.normalize_coords(&self.dtol);

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
                        let mut coord_delta = vec![tr_coord - coord_p];
                        normalize_frac_vectors(&mut coord_delta, &(2.0 * &frac_tols));

                        let cart_coord_delta =
                            Structure::get_cart_coords(&temp_structure.lattice, &coord_delta);

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
            Vector3::new(-1.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 1.0),
            Vector3::new(0.0, -1.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(-1.0, 1.0, 1.0),
            Vector3::new(1.0, -1.0, 1.0),
            Vector3::new(1.0, 1.0, -1.0),
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
                let mut coord_diff = vec![coord - new_coord];
                normalize_frac_vectors(&mut coord_diff, &frac_tols);
                let cart_coord_delta =
                    Structure::get_cart_coords(&temp_structure.lattice, &coord_diff);

                if cart_coord_delta[0].magnitude().abs() <= (self.dtol * 2.0) {
                    eq = true;
                    break;
                }
            }
            if !eq {
                let norm_coord = coord.clone();
                normalize_frac_vectors(&mut vec![norm_coord], &frac_tols);
                new_frac_coords.push(norm_coord);
                new_species.push(specie.clone());
            }
        }

        let mut prim_structure = Structure::new(new_lattice, new_species, new_frac_coords, false);
        prim_structure.normalize_coords(&self.dtol);

        return prim_structure;
    }

    /// Produce a new structure which is Niggli reduced
    /// using the algorithm by Grosse-Kuntsleve et al.
    /// [Acta Cryst. (2004). A60, 1-6](https://doi.org/10.1107/S010876730302186X)
    pub fn niggli_reduce(&self, structure: &Structure, tol: &f32) -> Structure {
        let epsilon = tol * structure.volume().powf(1.0 / 3.0);

        let mut metric_tensor = structure.metric_tensor();

        let mut transformation: Matrix3<f32> = Matrix3::identity();

        fn cust_eq(a: f32, b: f32, epsilon: f32) -> bool {
            return !((a < (b - epsilon)) || (b < (a - epsilon)));
        }

        let (mut a, mut b, mut c, mut xi, mut eta, mut zeta): (f32, f32, f32, f32, f32, f32);

        for iteration in 0..101 {
            if iteration == 100 {
                panic!("Could not find Niggli reduced structure!")
            }

            (a, b, c, xi, eta, zeta) = (
                metric_tensor.m11,
                metric_tensor.m22,
                metric_tensor.m33,
                2.0 * metric_tensor.m23,
                2.0 * metric_tensor.m13,
                2.0 * metric_tensor.m12,
            );

            // A1
            let c1: Matrix3<f32> = Matrix3::new(0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0);
            if ((a - epsilon) > b) || (cust_eq(a, b, epsilon) && ((xi.abs() - epsilon) > eta.abs()))
            {
                metric_tensor = c1.transpose() * metric_tensor * c1;
                transformation = transformation * c1;
            }

            (a, b, c, xi, eta, zeta) = (
                metric_tensor.m11,
                metric_tensor.m22,
                metric_tensor.m33,
                2.0 * metric_tensor.m23,
                2.0 * metric_tensor.m13,
                2.0 * metric_tensor.m12,
            );

            // A2
            let c2: Matrix3<f32> = Matrix3::new(-1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0);
            if ((b - epsilon) > c)
                || ((cust_eq(b, c, epsilon)) && ((eta.abs() - epsilon) > zeta.abs()))
            {
                metric_tensor = c2.transpose() * metric_tensor * c2;
                transformation = transformation * c2;
                continue;
            }

            // A3 + A4

            let mut num_negative = 0;
            let mut num_zero = 0;

            for val in [xi, eta, zeta] {
                if val < (-1.0 * epsilon) {
                    num_negative += 1;
                } else if val < epsilon {
                    num_zero += 1;
                }
            }

            let mut i = 1.0;
            let mut j = 1.0;
            let mut k = 1.0;

            if (num_negative % 2 == 0) && (num_zero == 0) {
                if xi < -epsilon {
                    i = -1.0;
                };
                if eta < -epsilon {
                    j = -1.0;
                };
                if zeta < -epsilon {
                    k = -1.0;
                };
            } else {
                let mut p: &mut f32 = &mut 1.0;
                let mut prod = 1.0;

                if (xi - epsilon) > 0.0 {
                    i = -1.0;
                    prod *= -1.0;
                } else if !(xi < -epsilon) {
                    p = &mut i
                };

                if (eta - epsilon) > 0.0 {
                    j = -1.0;
                    prod *= -1.0;
                } else if !(eta < -epsilon) {
                    p = &mut j
                };

                if (zeta - epsilon) > 0.0 {
                    k = -1.0;
                    prod *= -1.0;
                } else if !(zeta < -epsilon) {
                    p = &mut k
                };

                if prod < -epsilon {
                    *p = -1.0;
                };
            };

            let c3_4: Matrix3<f32> = Matrix3::new(i, 0.0, 0.0, 0.0, j, 0.0, 0.0, 0.0, k);
            metric_tensor = c3_4.transpose() * metric_tensor * c3_4;
            transformation = transformation * c3_4;

            (a, b, c, xi, eta, zeta) = (
                metric_tensor.m11,
                metric_tensor.m22,
                metric_tensor.m33,
                2.0 * metric_tensor.m23,
                2.0 * metric_tensor.m13,
                2.0 * metric_tensor.m12,
            );

            // A5

            if ((xi.abs() - epsilon) > b)
                || ((cust_eq(xi, b, epsilon)) && ((zeta - epsilon) > (2.0 * eta)))
                || ((cust_eq(xi, -1.0 * b, epsilon)) && (zeta < -epsilon))
            {
                let c5: Matrix3<f32> =
                    Matrix3::new(1.0, 0.0, 0.0, 0.0, 1.0, -1.0 * xi.signum(), 0.0, 0.0, 1.0);
                metric_tensor = c5.transpose() * metric_tensor * c5;
                transformation = transformation * c5;
                continue;
            }

            // A6

            if ((eta.abs() - epsilon) > a)
                || ((cust_eq(eta, a, epsilon)) && ((zeta - epsilon) > (2.0 * xi)))
                || ((cust_eq(eta, -1.0 * a, epsilon)) && (zeta < -epsilon))
            {
                let c6: Matrix3<f32> =
                    Matrix3::new(1.0, 0.0, -1.0 * eta.signum(), 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
                metric_tensor = c6.transpose() * metric_tensor * c6;
                transformation = transformation * c6;
                continue;
            }

            // A7

            if ((zeta.abs() - epsilon) > a)
                || ((cust_eq(zeta, a, epsilon)) && ((eta - epsilon) > (2.0 * xi)))
                || ((cust_eq(zeta, -1.0 * a, epsilon)) && (eta < -epsilon))
            {
                let c7: Matrix3<f32> =
                    Matrix3::new(1.0, -1.0 * zeta.signum(), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
                metric_tensor = c7.transpose() * metric_tensor * c7;
                transformation = transformation * c7;
                continue;
            }

            // A8

            let bool_sum = xi + eta + zeta + a + b;

            if (bool_sum < -epsilon)
                || ((cust_eq(bool_sum, 0.0, epsilon))
                    && (((2.0 * (a + eta)) + zeta - epsilon) > 0.0))
            {
                let c8: Matrix3<f32> = Matrix3::new(1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0);
                metric_tensor = c8.transpose() * metric_tensor * c8;
                transformation = transformation * c8;
                continue;
            }

            break;
        }

        (a, b, c, xi, eta, zeta) = (
            metric_tensor.m11,
            metric_tensor.m22,
            metric_tensor.m33,
            2.0 * metric_tensor.m23,
            2.0 * metric_tensor.m13,
            2.0 * metric_tensor.m12,
        );

        // Construct new structure
        let new_lattice = structure.lattice * transformation;

        // Does coord folding need to be considered?
        let new_frac_coords = Structure::get_frac_coords(&new_lattice, &structure.coords);

        let mut new_structure = Structure::new(
            new_lattice,
            structure.species.clone(),
            new_frac_coords,
            false,
        );

        new_structure.normalize_coords(&self.dtol);

        return new_structure;
    }

    /// Produce a new structure which is Delaunay reduced.
    pub fn delaunay_reduce(&self, structure: &Structure, full_reduction: bool) -> Structure {
        let mut new_structure = structure.clone();

        // Defined pair map for scalar prod vector
        let pair_map: HashMap<usize, (u8, u8)> = HashMap::from([
            (0, (0, 1)),
            (1, (0, 2)),
            (2, (0, 3)),
            (3, (1, 2)),
            (4, (1, 3)),
            (5, (2, 3)),
        ]);

        let mut delaunay_mat = structure.delaunay_matrix();

        let mut scalar_prods = Self::get_scalar_prods(&delaunay_mat);

        let mut count = 0;

        // Iterate until all scalar prods are <= zero
        while scalar_prods.iter().any(|val| val > &ZERO_TOL) && count <= 100 {
            let mut max_scalar_prod = &0.0;
            let mut max_index = 0 as usize;
            for (i, entry) in scalar_prods.iter().enumerate() {
                if entry > &ZERO_TOL && entry > &max_scalar_prod {
                    max_scalar_prod = entry;
                    max_index = i;
                }
            }

            let pos_pair = pair_map.get(&max_index).unwrap();

            // Delaunay transformation
            // e.g. If a*b > 0, the transformation is {a -> a, b -> -b, c -> c+B, d -> d+b}
            for i in 0..4 {
                if !(pos_pair.0 == i || pos_pair.1 == i) {
                    let new_col = &delaunay_mat.column(pos_pair.1 as usize)
                        + &delaunay_mat.column(i as usize);
                    delaunay_mat.set_column(i as usize, &new_col);
                };
            }
            let neg_column = &delaunay_mat.column(pos_pair.1 as usize) * (-1.0);
            delaunay_mat.set_column(pos_pair.1 as usize, &neg_column);

            scalar_prods = Self::get_scalar_prods(&delaunay_mat);

            count += 1;
        }

        if count == 100 {
            panic!("Reached max number of iterations without reduction!")
        }

        let new_lattice: Matrix3<f32>;

        if full_reduction {
            // New lattice vectors are chosen as shortest from {a, b, c, d, a+b, b+c, c+a}

            let mut base_vecs: Vec<Vector3<f32>> = delaunay_mat
                .column_iter()
                .map(|vec| Vector3::from(vec))
                .collect();

            // Sort quadruple to get three shortest vectors
            base_vecs.sort_by(|a, b| a.magnitude().partial_cmp(&b.magnitude()).unwrap());

            let mut first = base_vecs[0];
            let mut second = base_vecs[1];
            let mut third = base_vecs[2];
            let fourth = base_vecs[3];

            let g12 = first.dot(&second);
            let g13 = first.dot(&third);
            let g14 = first.dot(&fourth);
            let g23 = second.dot(&third);
            let g24 = second.dot(&fourth);

            // Check if tri-acute conditions are met and apply relevant transformations

            if -g12 < -(g13 + g14) {
                second = -1.0 * (first + second);
                first *= -1.0;
            } else if -g23 < -(g12 + g24) {
                third = -1.0 * (second + third);
                second *= -1.0;
            } else if -g13 < -(g12 + g14) {
                third = -1.0 * (first + third);
                first *= -1.0;
            }

            new_lattice = Matrix3::from_columns(&[first, second, third]);
        } else {
            new_lattice = Matrix3::from(delaunay_mat.fixed_slice::<3, 3>(0, 0));
        }

        // Does coord folding need to be considered?
        let new_frac_coords = Structure::get_frac_coords(&new_lattice, &new_structure.coords);

        new_structure = Structure::new(new_lattice, new_structure.species, new_frac_coords, false);
        new_structure.normalize_coords(&self.dtol);

        return new_structure;
    }

    /// Function to get three shortest (right-handed) cartestian translation vectors which still span the lattice.
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

            if cross_vec.magnitude() > self.dtol && first.dot(&cart_vec) > self.dtol {
                second_ind = vec_num + 1;
                break;
            }
        }

        for (vec_num, cart_vec) in cart_vecs[(second_ind + 1)..].iter().enumerate() {
            if cross_vec.dot(cart_vec) > self.dtol
                && cart_vec.dot(&first) > self.dtol
                && cart_vec.dot(&cart_vecs[second_ind]) > self.dtol
            {
                third_ind = vec_num + second_ind + 1;
                break;
            }
        }

        let final_vecs = vec![cart_vecs[0], cart_vecs[second_ind], cart_vecs[third_ind]];

        return final_vecs;
    }

    fn get_scalar_prods(delaunay_mat: &Matrix3x4<f32>) -> Vector6<f32> {
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
}
