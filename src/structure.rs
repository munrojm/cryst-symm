use std::collections::HashMap;
use std::fmt;
use std::string::String;

pub struct Structure {
    pub lattice: Vec<Vec<f32>>,
    pub species: Vec<String>,
    pub coords: Vec<Vec<f32>>,
    pub formula: String,
    pub reduced_formula: String,
}

impl Structure {
    pub fn new(lattice: Vec<Vec<f32>>, species: Vec<String>, coords: Vec<Vec<f32>>) -> Self {
        let (formula, reduced_formula) = Structure::get_formulas(species.clone());

        Self {
            lattice: lattice,
            species: species,
            coords: coords,
            formula: formula,
            reduced_formula: reduced_formula,
        }
    }

    fn species_coords_map(&self) -> HashMap<&String, &Vec<f32>> {
        let mut map = HashMap::new();

        for pair in self.species.iter().zip(&self.coords) {
            map.insert(pair.0, pair.1);
        }

        return map;
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

impl fmt::Display for Structure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "
Structure {:?} ({:?})

Lattice:
{:?},
{:?},
{:?}

Sites
{:?}",
            self.formula,
            self.reduced_formula,
            self.lattice[0],
            self.lattice[1],
            self.lattice[2],
            self.species_coords_map()
        )
    }
}
