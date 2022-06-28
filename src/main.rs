use crystsymm::reduce::Reducer;
use crystsymm::structure::Structure;
use crystsymm::symmetry::SymmetryAnalyzer;
use nalgebra::{Matrix3, Matrix3x4, Vector3};

fn main() {
    let a1 = Vector3::new(1.4214491371745352, 2.4620221259612376, 0.0);
    let a2 = Vector3::new(2.8428982743490705, 0.0, 0.0);
    let a3 = Vector3::new(1.4214491371745352, 0.8206740419870793, -4.715205165785214);

    let lattice = Matrix3::from_columns(&[a1, a2, a3]);

    let elements = ["Co", "Li", "O", "O"];
    let species: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
    let coords = vec![
        Vector3::new(0.5, 0.5, 0.5),
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.2396, 0.2396, 0.2812),
        Vector3::new(0.7604, 0.7604, 0.7188),
    ];

    let mut s = Structure::new(lattice, species, coords, false);
    println!["{:#?}", s];
    let reducer = Reducer { dtol: 0.1 };

    let new_s = reducer.find_primitive_cell(&s);

    //println!["{:#?}", new_s];

    //let reducer = Reducer { dtol: 0.1 };

    let d_s = reducer.delaunay_reduce(&new_s);

    //println!["{:#?}", d_s];

    let sa = SymmetryAnalyzer { dtol: 0.1 };
    let system = sa.get_crystal_system(&s);
    println!("{:?}", system);
}
