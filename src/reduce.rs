#![allow(clippy::neg_cmp_op_on_partial_ord)]
use crate::data::core::ZERO_TOL;
use crate::structure::Structure;
use crate::utils::{cust_eq, normalize_frac_vectors, num_negative_zero};
use nalgebra::{Matrix3, Matrix3x4, Vector3, Vector6};
use std::collections::HashMap;
use std::string::String;

#[derive(Debug, Clone)]
pub struct Reducer {
    pub dtol: f64,
    pub atol: f64,
}

impl Reducer {
    /// Produces a primitive structure using the following steps:
    ///
    /// 1. Find atom type with fewest sites.
    /// 2. Get all new potential unit translation vectors.
    /// 3. Select three shortest vectors which span the lattice.
    /// 4. Create new structure and fold atomic sites into new smaller cell.
    ///
    /// # Arguments
    ///
    /// * `structure` - Input structure to reduce.
    pub fn find_primitive_cell(&self, structure: &Structure) -> Structure {
        let mut temp_structure = structure.clone();

        let frac_tols = Vector3::from_iterator(
            temp_structure
                .lattice
                .column_iter()
                .map(|col| self.dtol / col.magnitude()),
        );

        temp_structure.normalize_coords(&self.dtol);

        //
        // 1.) Find atom type with the fewest sites
        //

        let (min_ele, ele_inds, _) = temp_structure.get_min_element();

        //
        // 2.) Shift origin to first atom of min_ele type and normalize coords (deprecated)
        //

        //let new_origin = temp_structure.frac_coords[*ele_inds.get(&min_ele).unwrap() as usize];
        //temp_structure.set_origin(new_origin);

        //
        // 3.) Get all potential new unit translation lattice vectors
        //
        let mut candidate_vecs: Vec<Vector3<f64>> = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];

        let base_site = temp_structure.frac_coords[*ele_inds.get(&min_ele).unwrap() as usize];

        for (coord_ind, coord) in temp_structure.frac_coords.iter().enumerate() {
            if (temp_structure.species[coord_ind] == min_ele)
                && (coord_ind != *ele_inds.get(&min_ele).unwrap() as usize)
            {
                let mut frac_diff = vec![coord - base_site];
                normalize_frac_vectors(&mut frac_diff, &frac_tols);
                candidate_vecs.push(frac_diff[0]);
            }
        }

        let mut translation_vecs: Vec<Vector3<f64>> = Vec::new();
        let mut all_matched: bool;
        let mut matched: bool;

        // TODO: Ensure inner loop is only over sites with the same specie type to speed up
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
                        normalize_frac_vectors(&mut coord_delta, &frac_tols);

                        let cart_coord_delta =
                            Structure::get_cart_coords(&temp_structure.lattice, &coord_delta);

                        if cart_coord_delta[0].magnitude().abs() <= (2.0 * self.dtol) {
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

        let mut cart_translation_vecs: Vec<Vector3<f64>> =
            Structure::get_cart_coords(&temp_structure.lattice, &translation_vecs);

        let frac_translation_vecs = self.get_shortest_translation_vecs(&mut cart_translation_vecs);

        //
        // 5.) Add combinations of shortest vectors
        //

        let smallest_frac_vecs =
            Structure::get_frac_coords(&temp_structure.lattice, &frac_translation_vecs);

        let frac_matrix = Matrix3::from_columns(&smallest_frac_vecs);

        let base_comb_vecs: Vec<Vector3<f64>> = vec![
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

        let mut combination_vecs: Vec<Vector3<f64>> = Vec::new();

        for comb_vec in base_comb_vecs {
            let transformed_vec = frac_matrix * comb_vec;
            combination_vecs.push(transformed_vec);
            combination_vecs.push(-1.0 * transformed_vec);

            let mut new_vec = transformed_vec;
            new_vec = -1.0 * new_vec;
            if new_vec.magnitude().abs() > ZERO_TOL {
                combination_vecs.push(new_vec);
            }
        }

        let mut cart_translation_vecs: Vec<Vector3<f64>> =
            Structure::get_cart_coords(&temp_structure.lattice, &combination_vecs);

        let final_cart_vecs = self.get_shortest_translation_vecs(&mut cart_translation_vecs);

        //
        // 6.) Create a new structure object with new lattice and coords.
        //     TODO: Average positions which fold into one another in the new cell.
        //

        let new_lattice = Matrix3::from_columns(&final_cart_vecs);

        let temp_frac_coords = Structure::get_frac_coords(&new_lattice, &temp_structure.coords);

        let mut new_frac_coords: Vec<Vector3<f64>> = Vec::new();

        let mut new_species: Vec<String> = Vec::new();

        for (coord, specie) in temp_frac_coords.iter().zip(&temp_structure.species) {
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
                let norm_coord = *coord;
                normalize_frac_vectors(&mut vec![norm_coord], &frac_tols);
                new_frac_coords.push(norm_coord);
                new_species.push(specie.clone());
            }
        }

        let mut prim_structure = Structure::new(new_lattice, new_species, new_frac_coords, false);
        prim_structure.normalize_coords(&self.dtol);

        prim_structure
    }

    /// Produce a new structure which is Niggli reduced
    /// using the algorithm by Grosse-Kuntsleve et al.
    /// [Acta Cryst. (2004). A60, 1-6](https://doi.org/10.1107/S010876730302186X)
    /// # Arguments
    ///
    /// * `structure` - Input structure to reduce.
    /// * `tol` - Value used to calculate tolerance (eps) as `eps = tol*(volume^(1/3))`
    pub fn niggli_reduce(&self, structure: &Structure, tol: &f64) -> Structure {
        let epsilon = tol * structure.volume().powf(1.0 / 3.0);

        let mut metric_tensor = structure.metric_tensor();

        let mut transformation: Matrix3<f64> = Matrix3::identity();

        let (mut a, mut b, mut c, mut xi, mut eta, mut zeta): (f64, f64, f64, f64, f64, f64);

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
            let c1: Matrix3<f64> = Matrix3::new(0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0);
            if ((a - epsilon) > b)
                || (cust_eq(&a, &b, &epsilon) && ((xi.abs() - epsilon) > eta.abs()))
            {
                metric_tensor = c1.transpose() * metric_tensor * c1;
                transformation *= c1;

                (_, b, c, xi, eta, zeta) = (
                    metric_tensor.m11,
                    metric_tensor.m22,
                    metric_tensor.m33,
                    2.0 * metric_tensor.m23,
                    2.0 * metric_tensor.m13,
                    2.0 * metric_tensor.m12,
                );
            }

            // A2
            let c2: Matrix3<f64> = Matrix3::new(-1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0);
            if ((b - epsilon) > c)
                || ((cust_eq(&b, &c, &epsilon)) && ((eta.abs() - epsilon) > zeta.abs()))
            {
                metric_tensor = c2.transpose() * metric_tensor * c2;
                transformation *= c2;
                continue;
            }

            // A3 + A4

            let (num_negative, num_zero) = num_negative_zero(&vec![xi, eta, zeta], &epsilon);

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
                let mut p: &mut f64 = &mut 1.0;
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

            let c3_4: Matrix3<f64> = Matrix3::new(i, 0.0, 0.0, 0.0, j, 0.0, 0.0, 0.0, k);
            metric_tensor = c3_4.transpose() * metric_tensor * c3_4;
            transformation *= c3_4;

            (a, b, _, xi, eta, zeta) = (
                metric_tensor.m11,
                metric_tensor.m22,
                metric_tensor.m33,
                2.0 * metric_tensor.m23,
                2.0 * metric_tensor.m13,
                2.0 * metric_tensor.m12,
            );

            // A5

            if ((xi.abs() - epsilon) > b)
                || ((cust_eq(&xi, &b, &epsilon)) && ((zeta - epsilon) > (2.0 * eta)))
                || ((cust_eq(&xi, &(-1.0 * b), &epsilon)) && (zeta < -epsilon))
            {
                let c5: Matrix3<f64> =
                    Matrix3::new(1.0, 0.0, 0.0, 0.0, 1.0, -1.0 * xi.signum(), 0.0, 0.0, 1.0);
                metric_tensor = c5.transpose() * metric_tensor * c5;
                transformation *= c5;
                continue;
            }

            // A6

            if ((eta.abs() - epsilon) > a)
                || ((cust_eq(&eta, &a, &epsilon)) && ((zeta - epsilon) > (2.0 * xi)))
                || ((cust_eq(&eta, &(-1.0 * a), &epsilon)) && (zeta < -epsilon))
            {
                let c6: Matrix3<f64> =
                    Matrix3::new(1.0, 0.0, -1.0 * eta.signum(), 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
                metric_tensor = c6.transpose() * metric_tensor * c6;
                transformation *= c6;
                continue;
            }

            // A7

            if ((zeta.abs() - epsilon) > a)
                || ((cust_eq(&zeta, &a, &epsilon)) && ((eta - epsilon) > (2.0 * xi)))
                || ((cust_eq(&zeta, &(-1.0 * a), &epsilon)) && (eta < -epsilon))
            {
                let c7: Matrix3<f64> =
                    Matrix3::new(1.0, -1.0 * zeta.signum(), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
                metric_tensor = c7.transpose() * metric_tensor * c7;
                transformation *= c7;
                continue;
            }

            // A8

            let bool_sum = xi + eta + zeta + a + b;

            if (bool_sum < -epsilon)
                || ((cust_eq(&bool_sum, &0.0, &epsilon))
                    && (((2.0 * (a + eta)) + zeta - epsilon) > 0.0))
            {
                let c8: Matrix3<f64> = Matrix3::new(1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0);
                metric_tensor = c8.transpose() * metric_tensor * c8;
                transformation *= c8;
                continue;
            }

            break;
        }

        // Construct new structure
        let new_lattice = structure.lattice * transformation;

        // Should coord folding (averaging) need to be considered?
        let new_frac_coords = Structure::get_frac_coords(&new_lattice, &structure.coords);

        let mut new_structure = Structure::new(
            new_lattice,
            structure.species.clone(),
            new_frac_coords,
            false,
        );

        new_structure.normalize_coords(&self.dtol);

        new_structure
    }

    /// Produce a new structure which is Delaunay reduced.
    ///
    /// **NOTE**: This method may need some changes. Use with caution.
    ///
    /// # Arguments
    ///
    /// * `structure` - Input structure to reduce.
    /// * `full_reduction` - Whether to choose the shortest vectors after Delaunay reduction for the new structure.
    /// If not, outputs the lattice directly after the Delaunay algorithm.
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
            let mut max_index = 0_usize;
            for (i, entry) in scalar_prods.iter().enumerate() {
                if entry > &ZERO_TOL && entry > max_scalar_prod {
                    max_scalar_prod = entry;
                    max_index = i;
                }
            }

            let pos_pair = pair_map.get(&max_index).unwrap();

            // Delaunay transformation
            // e.g. If a*b > 0, the transformation is {a -> a, b -> -b, c -> c+B, d -> d+b}
            for i in 0..4 {
                if !(pos_pair.0 == i || pos_pair.1 == i) {
                    let new_col =
                        delaunay_mat.column(pos_pair.1 as usize) + delaunay_mat.column(i as usize);
                    delaunay_mat.set_column(i as usize, &new_col);
                };
            }
            let neg_column = delaunay_mat.column(pos_pair.1 as usize) * (-1.0);
            delaunay_mat.set_column(pos_pair.1 as usize, &neg_column);

            scalar_prods = Self::get_scalar_prods(&delaunay_mat);

            count += 1;
        }

        if count == 100 {
            panic!("Reached max number of iterations without reduction!")
        }

        let new_lattice: Matrix3<f64>;

        if full_reduction {
            // New lattice vectors are chosen as shortest from {a, b, c, d, a+b, b+c, c+a}

            let mut base_vecs: Vec<Vector3<f64>> =
                delaunay_mat.column_iter().map(Vector3::from).collect();

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

        new_structure
    }

    /// Function to get three shortest (right-handed) cartestian translation vectors which still span the lattice.
    fn get_shortest_translation_vecs(&self, cart_vecs: &mut [Vector3<f64>]) -> Vec<Vector3<f64>> {
        cart_vecs.sort_by(|a, b| a.magnitude().partial_cmp(&b.magnitude()).unwrap());

        let first = cart_vecs[0];
        let mut cross_vec: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);
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
            if cross_vec.dot(cart_vec) > self.dtol {
                third_ind = vec_num + second_ind + 1;
                break;
            }
        }

        let final_vecs = vec![cart_vecs[0], cart_vecs[second_ind], cart_vecs[third_ind]];

        final_vecs
    }

    /// Function to get scalar products from the Delaunay matrix
    fn get_scalar_prods(delaunay_mat: &Matrix3x4<f64>) -> Vector6<f64> {
        let mut scalar_prods: Vector6<f64> = Vector6::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let mut ind = 0;
        for i in 0..3 {
            for j in (i + 1)..4 {
                scalar_prods[ind] = delaunay_mat.column(i).dot(&delaunay_mat.column(j));
                ind += 1;
            }
        }
        scalar_prods
    }
}

impl Default for Reducer {
    fn default() -> Reducer {
        Self {
            dtol: 0.05,
            atol: 5.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::approx_equal_iter;

    fn generate_structure() -> Structure {
        let a = Vector3::new(-3.748244, 0.0, 0.0);
        let b = Vector3::new(1.874122, -4.750729, 0.0);
        let c = Vector3::new(0.0, 1.5562529, 6.25794799);

        let lattice = Matrix3::from_columns(&[a, b, c]);

        let elements = ["Au", "Au", "Se", "Se"];
        let species: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
        let coords = vec![
            Vector3::new(-1.8741, 0.7781, 3.1290),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(-1.8741, -3.0213, 1.7273),
            Vector3::new(0.0, -0.1732, 4.5307),
        ];

        Structure::new(lattice, species, coords, true)
    }

    fn generate_reducer() -> Reducer {
        Reducer {
            ..Reducer::default()
        }
    }

    #[test]
    fn test_find_primitive() {
        let s = generate_structure();
        let mut supercell = s.clone();
        let transformation_matrix = Matrix3::new(2.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0);

        let r = generate_reducer();

        supercell.apply_transformation(&transformation_matrix, &0.05);

        let prim_s = r.find_primitive_cell(&supercell);

        assert!(approx_equal_iter(
            s.lattice.iter(),
            prim_s.lattice.iter(),
            &ZERO_TOL
        ));
    }

    #[test]
    fn test_niggli_reduce() {
        let s = generate_structure();
        let r = generate_reducer();

        let prim_s = r.find_primitive_cell(&s);
        let niggli_s = r.niggli_reduce(&prim_s, &1e-5);

        assert!(approx_equal_iter(
            niggli_s.lattice.iter(),
            Matrix3::new(-3.748244, 1.874122, 0.0, 0.0, 4.750729, -1.556253, 0.0, 0.0, -6.257948)
                .iter(),
            &ZERO_TOL
        ));
    }

    #[test]
    fn test_delaunay_reduce() {
        let s = generate_structure();
        let r = generate_reducer();

        let prim_s = r.find_primitive_cell(&s);
        let delaunay_s = r.delaunay_reduce(&prim_s, true);

        assert!(approx_equal_iter(
            delaunay_s.lattice.iter(),
            Matrix3::new(
                -3.748244, -1.874122, -1.874122, 0.0, 4.750729, 3.1944761, 0.0, 0.0, -6.257948
            )
            .iter(),
            &ZERO_TOL
        ));
    }
}
