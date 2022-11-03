use crate::data::core::ZERO_TOL;
use crate::utils::normalize_frac_vectors;
use itertools::izip;
use nalgebra::{try_invert_to, Matrix3, Matrix3x4, Vector3};
use std::collections::HashMap;
use std::f64::consts::PI;
use std::iter::FromIterator;
use std::string::String;

#[derive(Debug, Clone)]
///Representation of a periodic unit cell of a crystal structure.
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
    /// * `lattice` - A 3x3 matrix with columns containing the lattice vectors.
    /// * `species` - A list of strings representing the elements at each atomic site.
    /// * `coords` - A list of cartesian (Å) or fractional vector coordinates of the atomic sites.
    /// * `coords_are_cart` - Whether the input coordinates are in a cartesian basis.
    ///
    /// # Examples
    ///
    /// ```
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
        return ((self.lattice.column(0).dot(&self.lattice.column(1)))
            / (self.a().abs() * self.b().abs()))
        .acos()
            * (180.0 / PI);
    }

    /// Angle `beta` in degrees
    pub fn beta(&self) -> f64 {
        return ((self.lattice.column(1).dot(&self.lattice.column(2)))
            / (self.b().abs() * self.c().abs()))
        .acos()
            * (180.0 / PI);
    }

    /// Angle `gamma` in degrees
    pub fn gamma(&self) -> f64 {
        return ((self.lattice.column(0).dot(&self.lattice.column(2)))
            / (self.a().abs() * self.c().abs()))
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
    /// If the transformation results in a volume increase, then new atomic sites are searched for and added to the transformed cell.
    pub fn apply_transformation(&mut self, trans_mat: &Matrix3<f64>, pos_tol: &f64) {
        let frac_tols = Vector3::from_iterator(
            self.lattice
                .column_iter()
                .map(|col| pos_tol / col.magnitude()),
        );

        let mut inv_trans_mat: Matrix3<f64> = Matrix3::identity();

        let inverted = try_invert_to(*trans_mat, &mut inv_trans_mat);

        if !inverted {
            panic!("Transformation matrix is not invertible!");
        }

        let new_lattice = self.lattice * trans_mat;

        let mut new_frac_coords: Vec<Vector3<f64>> = Vec::new();
        let mut new_species: Vec<String> = Vec::new();

        let translation_vecs: Vec<Vector3<f64>> = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
        ];

        let volume_ratio = trans_mat.determinant().abs();

        // Search for new sites in transformed cell using volume ratio to determine proper number of sites
        let mut mult = 0;
        while (new_frac_coords.len() as f64 / self.coords.len() as f64) < volume_ratio && mult <= 10
        {
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

        if mult == 10 {
            panic!("Could not find all atomic sites in larger transformed cell.");
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

    /// Sets a new origin point for the unit cell.
    pub fn set_origin(&mut self, fractional_vector: Vector3<f64>) {
        for vec in self.frac_coords.iter_mut() {
            *vec -= fractional_vector;
        }
        self.coords = Self::get_cart_coords(&self.lattice, &self.frac_coords);
    }

    /// Normalizes fractional atomic positions to ensure all sites are within the unit cell.
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

    pub fn get_cart_coords(lattice: &Matrix3<f64>, coords: &[Vector3<f64>]) -> Vec<Vector3<f64>> {
        let mut new_coords = coords.to_owned();

        for (i, coord) in coords.iter().enumerate() {
            new_coords[i] = lattice * coord;
        }
        new_coords
    }

    pub fn get_frac_coords(lattice: &Matrix3<f64>, coords: &[Vector3<f64>]) -> Vec<Vector3<f64>> {
        let inverted_lattice = (*lattice)
            .pseudo_inverse(ZERO_TOL)
            .expect("Crystal lattice is not invertible!");

        let mut new_coords = coords.to_owned();

        for (i, coord) in coords.iter().enumerate() {
            new_coords[i] = inverted_lattice * coord;
        }
        new_coords
    }

    pub fn reciprocal_lattice(&self) -> Matrix3<f64> {
        let mut reciprocal_lattice = Matrix3::identity();
        let inverted = try_invert_to(self.lattice, &mut reciprocal_lattice);

        if !inverted {
            panic!("Crystal lattice is not invertible!");
        }

        reciprocal_lattice.transpose() * 2.0 * PI
    }

    /// Get the element with the smallest amount of sites, alongside an element-site-number map.
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

        for (element, count) in species_tally.into_iter() {
            formula.push_str(element);
            reduced_formula.push_str(element);

            if count != 1 {
                formula.push_str(&count.to_string());
            }

            let new_count = count / gcf;

            if new_count != 1 {
                reduced_formula.push_str(&new_count.to_string());
            }
        }
        (formula, reduced_formula)
    }

    /// Attempts to get the matrix that transforms this structure into the one provided.
    pub fn get_transformation_matrix(&self, structure: &Self) -> Matrix3<f64> {
        let inv_lattice: Matrix3<f64> = self
            .lattice
            .pseudo_inverse(ZERO_TOL)
            .expect("Crystal lattice is not invertible!");

        inv_lattice * structure.lattice
    }
}
