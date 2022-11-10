use crate::data::core::ZERO_TOL;
use crate::utils::normalize_frac_vectors;
use itertools::Itertools;
use itertools::{iproduct, izip};
use nalgebra::{try_invert_to, Matrix3, Matrix3x4, Vector3};
use std::collections::HashMap;
use std::f64::consts::PI;
use std::iter::FromIterator;
use std::string::String;

#[derive(Debug, Clone)]
/// Representation of a periodic unit cell of a crystal structure.
pub struct Structure {
    pub lattice: Matrix3<f64>,
    pub species: Vec<String>,
    pub coords: Vec<Vector3<f64>>,
    pub frac_coords: Vec<Vector3<f64>>,
}

impl Structure {
    /// Returns an instance of a crystal structure defined with a lattice, species, and coordinates.
    ///
    /// # Arguments
    ///
    /// * `lattice` - A 3x3 matrix with columns representing the lattice vectors.
    /// * `species` - A list of strings representing the elements at each atomic site.
    /// * `coords` - A list of cartesian (Å) or fractional vector coordinates of the atomic sites.
    /// * `coords_are_cart` - Whether the input coordinates are in a cartesian basis.
    ///
    /// # Examples
    ///
    /// ```
    /// # use crystsymm::structure::Structure;
    /// use nalgebra::{Matrix3, Vector3};
    ///
    /// let a = Vector3::new(-3.748244, 0.0, 0.0);
    /// let b = Vector3::new(1.874122, -4.750729, 0.0);
    /// let c = Vector3::new(0.0, 1.5562529, 6.25794799);
    ///
    /// let lattice = Matrix3::from_columns(&[a, b, c]);
    ///
    /// let elements = ["Au", "Au", "Se", "Se"];
    /// let species: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
    /// let coords = vec![
    ///     Vector3::new(-1.8741, 0.7781, 3.1290),
    ///     Vector3::new(0.0, 0.0, 0.0),
    ///     Vector3::new(-1.8741, -3.0213, 1.7273),
    ///     Vector3::new(0.0, -0.1732, 4.5307),
    /// ];
    ///
    /// let s = Structure::new(lattice, species, coords, true);
    /// ```
    pub fn new(
        lattice: Matrix3<f64>,
        species: Vec<String>,
        coords: Vec<Vector3<f64>>,
        coords_are_cart: bool,
    ) -> Structure {
        let frac_coords: Vec<Vector3<f64>>;
        let cart_coords: Vec<Vector3<f64>>;

        if coords_are_cart {
            frac_coords = Structure::get_frac_coords(&lattice, &coords);
            cart_coords = coords;
        } else {
            cart_coords = Structure::get_cart_coords(&lattice, &coords);
            frac_coords = coords;
        }

        Self {
            lattice,
            species,
            coords: cart_coords,
            frac_coords,
        }
    }

    /// Length of cell vector `a` in angstroms
    pub fn num_sites(&self) -> usize {
        self.species.len()
    }

    /// Length of cell vector `a` in angstroms
    pub fn a(&self) -> f64 {
        return self.lattice.column(0).magnitude();
    }

    /// Length of cell vector `b` in angstroms
    pub fn b(&self) -> f64 {
        return self.lattice.column(1).magnitude();
    }

    /// Length of cell vector `c` in angstroms
    pub fn c(&self) -> f64 {
        return self.lattice.column(2).magnitude();
    }

    /// Angle `alpha` in degrees
    pub fn alpha(&self) -> f64 {
        return ((self.lattice.column(1).dot(&self.lattice.column(2)))
            / (self.b().abs() * self.c().abs()))
        .acos()
            * (180.0 / PI);
    }

    /// Angle `beta` in degrees
    pub fn beta(&self) -> f64 {
        return ((self.lattice.column(0).dot(&self.lattice.column(2)))
            / (self.a().abs() * self.c().abs()))
        .acos()
            * (180.0 / PI);
    }

    /// Angle `gamma` in degrees
    pub fn gamma(&self) -> f64 {
        return ((self.lattice.column(0).dot(&self.lattice.column(1)))
            / (self.a().abs() * self.b().abs()))
        .acos()
            * (180.0 / PI);
    }

    /// Volume of the unit cell in angstroms
    pub fn volume(&self) -> f64 {
        self.lattice.determinant()
    }

    /// Metric tensor of lattice matrix
    pub fn metric_tensor(&self) -> Matrix3<f64> {
        self.lattice.transpose() * self.lattice
    }

    /// Applied a transformation to the structure as: `lattice * trans_mat = lattice'`
    ///
    /// If the transformation results in a volume increase, then new atomic sites are searched
    /// for and added to the transformed cell.
    ///
    /// # Arguments
    ///
    /// * `transformation_matrix` - A 3x3 matrix representing the transformation.
    /// * `pos_tol` - Position tolerance in angstroms to use when comparing atomic sites.
    pub fn apply_transformation(&mut self, transformation_matrix: &Matrix3<f64>, pos_tol: &f64) {
        let frac_tols = Vector3::from_iterator(
            self.lattice
                .column_iter()
                .map(|col| pos_tol / col.magnitude()),
        );

        let mut inv_trans_mat: Matrix3<f64> = Matrix3::identity();

        let inverted = try_invert_to(*transformation_matrix, &mut inv_trans_mat);

        if !inverted {
            panic!("Transformation matrix is not invertible!");
        }

        let new_lattice = self.lattice * transformation_matrix;

        let mut new_frac_coords: Vec<Vector3<f64>> = Vec::new();
        let mut new_species: Vec<String> = Vec::new();

        let volume_ratio = (transformation_matrix.determinant().abs() * 100.0).round() / 100.0;

        // Search for new sites in transformed cell using volume ratio to determine proper number of sites
        let mut mult = 2;
        while (new_frac_coords.len() as f64 / self.coords.len() as f64) < volume_ratio {
            let mut translation_vecs: Vec<Vector3<f64>> = Vec::new();

            for (i, j, k) in iproduct!(0..mult, 0..mult, 0..mult) {
                translation_vecs.push(Vector3::new(i as f64, j as f64, k as f64))
            }

            for (coord, specie) in self.frac_coords.iter().zip(&self.species) {
                for translation_vec in translation_vecs.iter() {
                    let mut frac_coord = *coord + (translation_vec * mult as f64);

                    frac_coord = inv_trans_mat * frac_coord;

                    let mut eq = false;
                    for new_coord in new_frac_coords.iter() {
                        let mut coord_diff = vec![frac_coord - new_coord];
                        normalize_frac_vectors(&mut coord_diff, &frac_tols);
                        let cart_coord_delta = Self::get_cart_coords(&self.lattice, &coord_diff);

                        if cart_coord_delta[0].magnitude().abs() <= (pos_tol * 2.0) {
                            eq = true;
                            break;
                        }
                    }
                    if !eq {
                        let norm_coord = frac_coord;
                        normalize_frac_vectors(&mut vec![norm_coord], &frac_tols);
                        new_frac_coords.push(norm_coord);
                        new_species.push(specie.clone());
                    }
                }
            }
            mult += 1;
        }

        self.lattice = new_lattice;
        self.frac_coords = new_frac_coords;
        self.coords = Self::get_cart_coords(&self.lattice, &self.frac_coords);
        self.species = new_species;
        self.normalize_coords(pos_tol);
    }

    /// Computes the Delaunay matrix from the lattice matrix of the structure.
    pub fn delaunay_matrix(&self) -> Matrix3x4<f64> {
        let d4 = -1.0 * (self.lattice.column(0) + self.lattice.column(1) + self.lattice.column(2));

        let cols: Vec<Vector3<f64>> = self
            .lattice
            .column_iter()
            .map(|vec| Vector3::new(vec[0], vec[1], vec[2]))
            .chain([d4])
            .collect();

        Matrix3x4::from_columns(&cols)
    }

    /// Normalizes fractional atomic positions to ensure all sites are within the unit cell.
    ///
    /// # Arguments
    ///
    /// * `pos_tol` - Position tolerance in angstroms to use when comparing atomic sites.
    pub fn normalize_coords(&mut self, pos_tol: &f64) {
        let frac_tols = Vector3::from_iterator(
            self.lattice
                .column_iter()
                .map(|col| pos_tol / col.magnitude()),
        );

        normalize_frac_vectors(&mut self.frac_coords, &frac_tols);

        let cart_coords = Self::get_cart_coords(&self.lattice, &self.frac_coords);

        // Sort coords according to species type to keep organized
        let mut zipped_coords: Vec<(&String, &Vector3<f64>, &Vector3<f64>)> =
            izip!(&self.species, &cart_coords, &self.frac_coords).collect();

        zipped_coords.sort_by(|a, b| a.0.cmp(b.0));

        let new_species: Vec<String> = zipped_coords.iter().map(|tuple| tuple.0.clone()).collect();

        let new_cart_coords: Vec<Vector3<f64>> =
            zipped_coords.iter().map(|tuple| *tuple.1).collect();

        let new_frac_coords: Vec<Vector3<f64>> =
            zipped_coords.iter().map(|tuple| *tuple.2).collect();

        self.species = new_species;
        self.coords = new_cart_coords;
        self.frac_coords = new_frac_coords;
    }

    /// Compute coordinates in a cartesian basis from the provided lattice and fractional coordinate list.
    ///
    /// # Arguments
    ///
    /// * `lattice` - A 3x3 matrix with columns representing the lattice vectors.
    /// * `frac_coords` - A list of fractional vector coordinates of the atomic sites.
    pub fn get_cart_coords(
        lattice: &Matrix3<f64>,
        frac_coords: &[Vector3<f64>],
    ) -> Vec<Vector3<f64>> {
        let mut new_coords = frac_coords.to_owned();

        for (i, coord) in frac_coords.iter().enumerate() {
            new_coords[i] = lattice * coord;
        }
        new_coords
    }

    /// Compute coordinates in a fractional basis from the provided lattice and cartesian coordinate list.
    ///
    /// # Arguments
    ///
    /// * `lattice` - A 3x3 matrix with columns representing the lattice vectors.
    /// * `cart_coords` - A list of cartesian (Å) vector coordinates of the atomic sites.
    pub fn get_frac_coords(
        lattice: &Matrix3<f64>,
        cart_coords: &[Vector3<f64>],
    ) -> Vec<Vector3<f64>> {
        let inverted_lattice = (*lattice)
            .pseudo_inverse(ZERO_TOL)
            .expect("Crystal lattice is not invertible!");

        let mut new_coords = cart_coords.to_owned();

        for (i, coord) in cart_coords.iter().enumerate() {
            new_coords[i] = inverted_lattice * coord;
        }
        new_coords
    }

    /// Reciprocal lattice of the structure
    pub fn reciprocal_lattice(&self) -> Matrix3<f64> {
        let mut reciprocal_lattice = Matrix3::identity();
        let inverted = try_invert_to(self.lattice, &mut reciprocal_lattice);

        if !inverted {
            panic!("Crystal lattice is not invertible!");
        }

        reciprocal_lattice.transpose() * 2.0 * PI
    }

    /// Get the element with the smallest amount of sites, alongside element-site index and element-count maps.
    pub fn get_min_element(&self) -> (String, HashMap<String, u16>, HashMap<String, u16>) {
        let mut type_count: HashMap<String, u16> = HashMap::new();
        let mut ele_inds: HashMap<String, u16> = HashMap::new();
        let temp_ele = self.species.clone();

        for (ele_ind, ele) in temp_ele.iter().enumerate() {
            *type_count.entry(ele.clone()).or_insert(0) += 1;
            ele_inds.entry(ele.clone()).or_insert(ele_ind as u16);
        }

        // Sorting ensures determinism when getting min
        let mut sorted_pairs: Vec<(String, u16)> = Vec::from_iter(type_count.clone().into_iter());

        sorted_pairs.sort_by(|a, b| a.0.cmp(&b.0));

        let min_ele = &sorted_pairs
            .iter()
            .min_by_key(|entry| entry.1)
            .unwrap()
            .0
            .clone();

        (min_ele.clone(), ele_inds, type_count)
    }

    /// Formula and reduced formula of the structure.
    pub fn formulas(&self) -> (String, String) {
        let mut species_tally = HashMap::<&String, i8>::new();

        let mut max_count = 0;

        for specie in self.species.iter() {
            let count = species_tally.entry(specie).or_insert(0);
            *count += 1;

            if count > &mut max_count {
                max_count = *count;
            }
        }
        let mut gcf: i8 = 1;

        for div in 2..(max_count + 1) {
            if species_tally.values().all(|&tally| (tally % div) == 0) {
                gcf = div;
            }
        }

        let mut formula = String::new();
        let mut reduced_formula = String::new();

        species_tally
            .into_iter()
            .sorted()
            .for_each(|(element, count)| {
                formula.push_str(element);
                reduced_formula.push_str(element);

                if count != 1 {
                    formula.push_str(&count.to_string());
                }

                let new_count = count / gcf;

                if new_count != 1 {
                    reduced_formula.push_str(&new_count.to_string());
                }
            });

        (formula, reduced_formula)
    }

    /// Attempts to get the matrix that transforms this structure into the one provided.
    ///
    /// # Arguments
    ///
    /// * `structure` - A transformed structure
    pub fn get_transformation_matrix(&self, structure: &Structure) -> Matrix3<f64> {
        let inv_lattice: Matrix3<f64> = self
            .lattice
            .pseudo_inverse(ZERO_TOL)
            .expect("Crystal lattice is not invertible!");

        inv_lattice * structure.lattice
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

    #[test]
    fn test_scalar_attributes() {
        let s = generate_structure();

        let real_scalars = [
            4.0,
            3.748244,
            5.107030380008033,
            6.44855306023863,
            102.97327735350085,
            90.0,
            111.52881069906073,
            111.43460086012757,
        ];
        let test_scalars = [
            s.num_sites() as f64,
            s.a(),
            s.b(),
            s.c(),
            s.alpha(),
            s.beta(),
            s.gamma(),
            s.volume(),
        ];

        assert!(approx_equal_iter(
            real_scalars.into_iter(),
            test_scalars.into_iter(),
            &ZERO_TOL
        ));
    }

    #[test]
    fn test_formulas() {
        let mut s = generate_structure();
        let (formula, reduced_formula) = s.formulas();
        assert_eq!(formula, "Au2Se2".to_string());
        assert_eq!(reduced_formula, "AuSe".to_string());

        let elements = [
            "Si", "Si", "Si", "Si", "Si", "Si", "O", "O", "P", "P", "P", "P",
        ];
        s.species = elements.iter().map(|s| s.to_string()).collect();

        let (formula, reduced_formula) = s.formulas();
        assert_eq!(formula, "O2P4Si6".to_string());
        assert_eq!(reduced_formula, "OP2Si3".to_string());
    }

    #[test]
    fn test_min_element() {
        let s = generate_structure();

        let ele_map = HashMap::from([("Au".to_string(), 2), ("Se".to_string(), 2)]);
        let ind_map = HashMap::from([("Au".to_string(), 0), ("Se".to_string(), 2)]);

        let (min_ele, ele_inds, type_count) = s.get_min_element();

        assert_eq!(min_ele, "Au");
        assert_eq!(type_count, ele_map);
        assert_eq!(ele_inds, ind_map);
    }

    #[test]
    fn test_metric_tensor() {
        let s = generate_structure();

        assert!(approx_equal_iter(
            s.metric_tensor().into_iter().copied(),
            Matrix3::new(
                14.04933308,
                -7.02466654,
                0.0,
                -7.02466654,
                26.0817593,
                -7.39333626,
                0.0,
                -7.39333626,
                41.58383657
            )
            .into_iter()
            .copied(),
            &ZERO_TOL
        ));
    }

    #[test]
    fn test_reciprocal_lattice() {
        let s = generate_structure();

        assert!(approx_equal_iter(
            s.reciprocal_lattice().into_iter().copied(),
            Matrix3::new(
                -1.67630104,
                -0.66128644,
                0.16445151,
                0.0,
                -1.32257287,
                0.32890302,
                0.0,
                0.0,
                1.00403284
            )
            .transpose()
            .into_iter()
            .copied(),
            &ZERO_TOL
        ));
    }

    #[test]
    fn test_delaunay() {
        let s = generate_structure();

        assert!(approx_equal_iter(
            s.delaunay_matrix().into_iter().copied(),
            Matrix3x4::new(
                -3.748244,
                1.874122,
                0.0,
                1.874122,
                0.0,
                -4.750729,
                1.5562529,
                3.1944760999999997,
                0.0,
                0.0,
                6.25794799,
                -6.25794799
            )
            .into_iter()
            .copied(),
            &ZERO_TOL
        ));
    }

    #[test]
    fn test_get_cart_coords() {
        let s = generate_structure();

        let cart_coords = Structure::get_cart_coords(&s.lattice, &s.frac_coords);

        let test_coords: Vec<Vector3<f64>> = vec![
            Vector3::new(-1.8741, 0.7781, 3.1290),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(-1.8741, -3.0213, 1.7273),
            Vector3::new(0.0, -0.1732, 4.5307),
        ];

        for (real, test) in cart_coords.into_iter().zip(test_coords) {
            assert!(approx_equal_iter(
                real.iter().copied(),
                test.iter().copied(),
                &1e-4
            ));
        }
    }

    #[test]
    fn test_get_frac_coords() {
        let s = generate_structure();

        let frac_coords = Structure::get_frac_coords(&s.lattice, &s.coords);

        let test_coords: Vec<Vector3<f64>> = vec![
            Vector3::new(0.4999, 0.0000, 0.5),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.863192, 0.726384, 0.276014),
            Vector3::new(0.136808, 0.273616, 0.723986),
        ];

        for (real, test) in frac_coords.into_iter().zip(test_coords) {
            assert!(approx_equal_iter(
                real.iter().copied(),
                test.iter().copied(),
                &1e-4
            ));
        }
    }

    #[test]
    fn test_normalize_coords() {
        let mut s = generate_structure();

        let test_coords = vec![
            Vector3::new(1.05, 2.6, 1.5),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(-0.75, 0.5, 0.5),
            Vector3::new(-1.05, -2.61, -1.5),
        ];

        let real_coords = vec![
            Vector3::new(0.05, 0.6, 0.5),
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.25, 0.5, 0.5),
            Vector3::new(0.95, 0.39, 0.5),
        ];

        s.frac_coords = test_coords;

        s.normalize_coords(&0.1);

        for (real, test) in s.frac_coords.into_iter().zip(real_coords) {
            assert!(approx_equal_iter(
                real.iter().copied(),
                test.iter().copied(),
                &1e-4
            ));
        }
    }

    #[test]
    fn test_apply_transformation() {
        let mut s = generate_structure();
        println!("{:?}", s.lattice);
        let trans_mat = Matrix3::new(2.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

        s.apply_transformation(&trans_mat, &0.01);

        assert!(approx_equal_iter(
            s.lattice.into_iter().copied(),
            Matrix3::new(-7.496488, -1.874122, 0.0, 0.0, -4.750729, 1.556253, 0.0, 0.0, 6.257948)
                .into_iter()
                .copied(),
            &ZERO_TOL
        ));
    }

    #[test]
    fn test_get_transformation_matrix() {
        let s = generate_structure();
        let mut trans_s = s.clone();

        let trans_mat = Matrix3::new(2.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

        trans_s.apply_transformation(&trans_mat, &0.1);

        let recovered_mat = s.get_transformation_matrix(&trans_s);

        assert!(approx_equal_iter(
            recovered_mat.into_iter().copied(),
            trans_mat.into_iter().copied(),
            &ZERO_TOL
        ));
    }
}
