pub mod analyzer;
pub mod data;
pub mod pointgroup;
pub mod reduce;
pub mod spacegroup;
pub mod structure;
pub mod symmop;
pub mod utils;

use analyzer::SymmetryAnalyzer;
use nalgebra::{Matrix3, Vector3};
use reduce::Reducer;
use structure::Structure;

use pyo3::prelude::*;

fn generate_input_structure(
    lattice: Vec<f64>,
    species: Vec<String>,
    coords: Vec<Vec<f64>>,
    coords_are_cart: bool,
) -> Structure {
    let formatted_lattice = Matrix3::from_iterator(lattice.into_iter()); // Lattice needs to be column-major iterator

    let formatted_coords: Vec<Vector3<f64>> = coords
        .iter()
        .map(|vec| Vector3::new(vec[0], vec[1], vec[2]))
        .collect();

    Structure::new(
        formatted_lattice,
        species,
        formatted_coords,
        coords_are_cart,
    )
}

fn generate_output_structure_data(structure: Structure) -> (Vec<f64>, Vec<String>, Vec<Vec<f64>>) {
    let lattice_vec: Vec<f64> = structure.lattice.iter().copied().collect();

    let mut coords_vecs: Vec<Vec<f64>> = Vec::new();

    for vec in structure.frac_coords.iter() {
        let new_vec: Vec<f64> = Vec::from_iter(vec.iter().copied());
        coords_vecs.push(new_vec);
    }

    (lattice_vec, structure.species, coords_vecs)
}

#[pyfunction]
fn find_primitive(
    lattice: Vec<f64>,
    species: Vec<String>,
    coords: Vec<Vec<f64>>,
    coords_are_cart: bool,
    dtol: f64,
    atol: f64,
) -> PyResult<(Vec<f64>, Vec<String>, Vec<Vec<f64>>)> {
    let structure = generate_input_structure(lattice, species, coords, coords_are_cart);

    let reducer = Reducer { dtol, atol };

    let prim_structure = reducer.find_primitive_cell(&structure);

    let output_data = generate_output_structure_data(prim_structure);

    Ok(output_data)
}

#[pyfunction]
fn niggli_reduce(
    lattice: Vec<f64>,
    species: Vec<String>,
    coords: Vec<Vec<f64>>,
    coords_are_cart: bool,
    tol: f64,
) -> PyResult<(Vec<f64>, Vec<String>, Vec<Vec<f64>>)> {
    let structure = generate_input_structure(lattice, species, coords, coords_are_cart);

    let reducer = Reducer {
        ..Default::default()
    };

    let niggli_structure = reducer.niggli_reduce(&structure, &tol);

    let output_data = generate_output_structure_data(niggli_structure);

    Ok(output_data)
}

#[pyfunction]
fn get_space_group_info(
    lattice: Vec<f64>,
    species: Vec<String>,
    coords: Vec<Vec<f64>>,
    coords_are_cart: bool,
    dtol: f64,
    atol: f64,
) -> PyResult<(u8, String, Vec<(Vec<i8>, Vec<f64>)>)> {
    let structure = generate_input_structure(lattice, species, coords, coords_are_cart);

    let sa = SymmetryAnalyzer { dtol, atol };

    let (space_group_num, space_group_symbol, space_group) = sa.get_space_group_info(&structure);

    let mut op_list: Vec<(Vec<i8>, Vec<f64>)> = Vec::new();
    for symm_op in space_group.operations.iter() {
        let mat_list: Vec<i8> = symm_op.rotation.transpose().iter().copied().collect();
        let vec_list: Vec<f64> = symm_op.translation.iter().copied().collect();
        op_list.push((mat_list, vec_list));
    }

    Ok((space_group_num, space_group_symbol, op_list))
}

#[pyfunction]
fn get_standard_conventional_structure(
    lattice: Vec<f64>,
    species: Vec<String>,
    coords: Vec<Vec<f64>>,
    coords_are_cart: bool,
    dtol: f64,
    atol: f64,
) -> PyResult<(Vec<f64>, Vec<String>, Vec<Vec<f64>>)> {
    let structure = generate_input_structure(lattice, species, coords, coords_are_cart);

    let sa = SymmetryAnalyzer { dtol, atol };
    let conv_struct = sa.get_standard_conventional_structure(&structure);

    let output_data = generate_output_structure_data(conv_struct);

    Ok(output_data)
}

#[pyfunction]
fn get_standard_primitive_structure(
    lattice: Vec<f64>,
    species: Vec<String>,
    coords: Vec<Vec<f64>>,
    coords_are_cart: bool,
    dtol: f64,
    atol: f64,
) -> PyResult<(Vec<f64>, Vec<String>, Vec<Vec<f64>>)> {
    let structure = generate_input_structure(lattice, species, coords, coords_are_cart);

    let sa = SymmetryAnalyzer { dtol, atol };
    let conv_struct = sa.get_standard_primitive_structure(&structure);

    let output_data = generate_output_structure_data(conv_struct);

    Ok(output_data)
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn crystsymm(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(niggli_reduce, m)?)?;
    m.add_function(wrap_pyfunction!(get_standard_conventional_structure, m)?)?;
    m.add_function(wrap_pyfunction!(get_standard_primitive_structure, m)?)?;
    m.add_function(wrap_pyfunction!(get_space_group_info, m)?)?;
    m.add_function(wrap_pyfunction!(find_primitive, m)?)?;

    Ok(())
}
