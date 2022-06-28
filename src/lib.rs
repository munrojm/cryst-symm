pub mod data;
pub mod reduce;
pub mod structure;
pub mod symmetry;

use nalgebra::{Matrix3, Vector3};
use structure::Structure;
use symmetry::SymmetryAnalyzer;

use pyo3::prelude::*;

#[pyfunction]
fn get_crystal_system(
    lattice: Vec<f32>,
    species: Vec<String>,
    coords: Vec<Vec<f32>>,
    coords_are_cart: bool,
) -> PyResult<(String, String)> {
    let formatted_lattice = Matrix3::from_iterator(lattice.into_iter());
    let formatted_coords: Vec<Vector3<f32>> = coords
        .iter()
        .map(|vec| Vector3::new(vec[0], vec[1], vec[2]))
        .collect();

    let structure = Structure::new(
        formatted_lattice,
        species,
        formatted_coords,
        coords_are_cart,
    );
    let sa = SymmetryAnalyzer { dtol: 0.1 };
    let result = sa.get_crystal_system(&structure).unwrap();
    Ok((result.0.clone(), result.1.clone()))
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn crystsymm(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_crystal_system, m)?)?;
    Ok(())
}
