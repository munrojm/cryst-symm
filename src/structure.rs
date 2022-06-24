use nalgebra::{try_invert_to, Matrix3, Vector3};
use std::collections::HashMap;
use std::f32::consts::PI;
use std::string::String;

#[derive(Debug)]
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
