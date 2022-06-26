use crystsymm::reduce::Reducer;
use crystsymm::structure::Structure;
use nalgebra::{Matrix3, Vector3};

fn main() {
    let a1 = Vector3::new(0.0, 9.501458, 0.0);
    let a2 = Vector3::new(3.748244, 0.0, 0.0);
    let a3 = Vector3::new(0.0, -1.5562529999999997, -6.257947999999999);

    let lattice = Matrix3::from_columns(&[a1, a2, a3]);

    let elements = ["Au", "Au", "Au", "Au", "Se", "Se", "Se", "Se"];
    let species: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
    let coords = vec![
        Vector3::new(1.8741, -0.7781, -3.1290),
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0000, 3.9726, -3.1290),
        Vector3::new(1.8741, 4.7507, 0.0000),
        Vector3::new(1.8741, 4.9239, -4.5307),
        Vector3::new(0.0000, 7.7720, -1.7273),
        Vector3::new(0.0000, 0.1732, -4.5307),
        Vector3::new(1.8741, 3.0213, -1.7273),
    ];

    let mut s = Structure::new(lattice, species, coords, true);
    println!["{:#?}", s];

    s.set_origin(Vector3::new(0.5, 0.0, 0.5));
    //println!["{:#?}", s];

    s.normalize_coords(0.1);
    //println!["{:#?}", s];

    let reducer = Reducer { dtol: 0.1 };

    let new_s = reducer.find_primitive_cell(&s);

    println!["{:#?}", new_s];
}
