use nalgebra::{try_invert_to, Matrix3, Vector3};
use std::collections::HashMap;
use std::f32::consts::PI;
use std::string::String;

#[derive(Debug)]
///Representation of a periodic unit cell of a crystal structure.
pub struct Structure {
    pub lattice: Matrix3<f32>,
    pub reciprocal_lattice: Matrix3<f32>,
    pub species: Vec<String>,
    pub coords: Vec<Vector3<f32>>,
    pub frac_coords: Vec<Vector3<f32>>,
    pub formula: String,
    pub reduced_formula: String,
}

impl Structure {
    /// Returns an instance of a crystal structure defined with a lattice, species, and coordinates.
    ///
    /// # Arguments
    ///
    /// * `lattice` - A 3x3 matrix with columns containing the lattice vectors.
    /// * `species` - A list of strings representing the elements at each atomic site.
    /// * `coords` - A list of cartesiean vector coordinates of the atomic sites.
    ///
    /// # Examples
    ///
    /// ```
    /// use nalgebra::{Matrix3, Vector3};
    ///
    /// let a1 = Vector3::new(-3.748244, 0.0, 0.0);
    /// let a2 = Vector3::new(1.874122, -4.750729, 0.0);
    /// let a3 = Vector3::new(0.0, 1.5562529, 6.25794799);
    ///
    /// let lattice = Matrix3::from_columns(&[a1, a2, a3]);
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
    /// let s = Structure::new(lattice, species, coords);
    /// ```
    pub fn new(
        lattice: Matrix3<f32>,
        species: Vec<String>,
        coords: Vec<Vector3<f32>>,
    ) -> Structure {
        let (formula, reduced_formula) = Structure::get_formulas(species.clone());
        let reciprocal_lattice = Structure::get_reciprocal_lattice(lattice.clone());
        let frac_coords = Structure::get_frac_coords(lattice.clone(), coords.clone());
        Self {
            lattice: lattice,
            reciprocal_lattice: reciprocal_lattice,
            species: species,
            coords: coords,
            frac_coords: frac_coords,
            formula: formula,
            reduced_formula: reduced_formula,
        }
    }

    pub fn set_origin(&mut self, fractional_vector: Vector3<f32>) {
        for vec in self.frac_coords.iter_mut() {
            *vec -= fractional_vector;
        }

        self.normalize_coords()
    }

    pub fn normalize_coords(&mut self) {
        for vec in self.frac_coords.iter_mut() {
            for ind in vec {
                *ind = ((*ind % 1.0) + 1.0) % 1.0;
            }
        }

        let new_cart_coords = Self::get_cart_coords(self.lattice, self.coords.clone());

        self.coords = new_cart_coords;
    }

    fn get_cart_coords(lattice: Matrix3<f32>, mut coords: Vec<Vector3<f32>>) -> Vec<Vector3<f32>> {
        for coord in coords.iter_mut() {
            *coord = lattice * *coord;
        }
        return coords;
    }

    fn get_frac_coords(lattice: Matrix3<f32>, mut coords: Vec<Vector3<f32>>) -> Vec<Vector3<f32>> {
        let mut inverted_lattice = Matrix3::identity();
        let inverted = try_invert_to(lattice, &mut inverted_lattice);

        if !inverted {
            panic!("Crystal lattice is not invertible!");
        }

        for coord in coords.iter_mut() {
            *coord = inverted_lattice * *coord;
        }
        return coords;
    }

    fn _species_coords_map(&self) -> HashMap<&String, &Vector3<f32>> {
        let mut map = HashMap::new();

        for pair in self.species.iter().zip(&self.coords) {
            map.insert(pair.0, pair.1);
        }

        return map;
    }

    fn get_reciprocal_lattice(lattice: Matrix3<f32>) -> Matrix3<f32> {
        let mut reciprocal_lattice = Matrix3::identity();
        let inverted = try_invert_to(lattice, &mut reciprocal_lattice);

        if !inverted {
            panic!("Crystal lattice is not invertible!");
        }

        return reciprocal_lattice.transpose() * 2.0 * PI;
    }

    fn get_formulas(species: Vec<String>) -> (String, String) {
        let mut species_tally = HashMap::<&String, i8>::new();

        let mut max_count = 0;

        for specie in species.iter() {
            let count = species_tally.entry(specie).or_insert(0);
            *count += 1;

            if count > &mut max_count {
                max_count = count.clone();
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

            let new_count = &count / &gcf;

            if new_count != 1 {
                reduced_formula.push_str(&new_count.to_string());
            }
        }
        return (formula, reduced_formula);
    }
}

// impl fmt::Display for Structure {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         write!(
//             f,
//             "
// Structure {:?} ({:?})

// Lattice:
// {:?}
// {:?}
// {:?}

// Sites:
// {:?}
// ",
//             self.formula,
//             self.reduced_formula,
//             self.lattice[0],
//             self.lattice[1],
//             self.lattice[2],
//             self.species_coords_map().iter()
//         )
//     }
// }
