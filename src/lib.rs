pub mod analyzer;
pub mod data;
pub mod reduce;
pub mod structure;
pub mod utils;

use analyzer::SymmetryAnalyzer;
use nalgebra::{Matrix3, Vector3};
use reduce::Reducer;
use structure::Structure;

use pyo3::prelude::*;

#[pyfunction]
fn niggli_reduce(
    lattice: Vec<f32>,
    species: Vec<String>,
    coords: Vec<Vec<f32>>,
    coords_are_cart: bool,
    dtol: f32,
    atol: f32,
    tol: f32,
) -> PyResult<(Vec<f32>, Vec<String>, Vec<Vec<f32>>)> {
    let formatted_lattice = Matrix3::from_iterator(lattice.into_iter()); // Lattice needs to be column-major iterator
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
    let reducer = Reducer {
        dtol: dtol,
        atol: atol,
    };

    let niggli_structure = reducer.niggli_reduce(&structure, &tol);

    let lattice_vec: Vec<f32> = niggli_structure.lattice.iter().map(|x| *x).collect();
    let mut coords_vec: Vec<Vec<f32>> = Vec::new();
    for vec in niggli_structure.frac_coords.iter() {
        let new_vec: Vec<f32> = vec.iter().map(|x| *x).collect();
        coords_vec.push(new_vec);
    }
    Ok((lattice_vec, niggli_structure.species, coords_vec))
}

#[pyfunction]
fn get_standard_conventional_structure(
    lattice: Vec<f32>,
    species: Vec<String>,
    coords: Vec<Vec<f32>>,
    coords_are_cart: bool,
    dtol: f32,
    atol: f32,
) -> PyResult<(Vec<f32>, Vec<String>, Vec<Vec<f32>>)> {
    let formatted_lattice = Matrix3::from_iterator(lattice.into_iter()); // Lattice needs to be column-major iterator

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
    let sa = SymmetryAnalyzer {
        dtol: dtol,
        atol: atol,
    };
    let conv_struct = sa.get_standard_conventional_structure(&structure);

    let lattice_vec: Vec<f32> = conv_struct.lattice.iter().map(|x| *x).collect();

    let mut coords_vec: Vec<Vec<f32>> = Vec::new();
    for vec in conv_struct.frac_coords.iter() {
        let new_vec: Vec<f32> = vec.iter().map(|x| *x).collect();
        coords_vec.push(new_vec);
    }
    Ok((lattice_vec, conv_struct.species, coords_vec))
}

#[pyfunction]
fn get_standard_primitive_structure(
    lattice: Vec<f32>,
    species: Vec<String>,
    coords: Vec<Vec<f32>>,
    coords_are_cart: bool,
    dtol: f32,
    atol: f32,
) -> PyResult<(Vec<f32>, Vec<String>, Vec<Vec<f32>>)> {
    let formatted_lattice = Matrix3::from_iterator(lattice.into_iter()); // Lattice needs to be column-major iterator

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
    let sa = SymmetryAnalyzer {
        dtol: dtol,
        atol: atol,
    };
    let conv_struct = sa.get_standard_primitive_structure(&structure);

    let lattice_vec: Vec<f32> = conv_struct.lattice.iter().map(|x| *x).collect();

    let mut coords_vec: Vec<Vec<f32>> = Vec::new();
    for vec in conv_struct.frac_coords.iter() {
        let new_vec: Vec<f32> = vec.iter().map(|x| *x).collect();
        coords_vec.push(new_vec);
    }
    Ok((lattice_vec, conv_struct.species, coords_vec))
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn crystsymm(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(niggli_reduce, m)?)?;
    m.add_function(wrap_pyfunction!(get_standard_conventional_structure, m)?)?;
    m.add_function(wrap_pyfunction!(get_standard_primitive_structure, m)?)?;
    Ok(())
}
