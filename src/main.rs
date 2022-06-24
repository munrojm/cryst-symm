use crystsymm::Structure;
use nalgebra::{Matrix3, Vector3};

fn main() {
    let a1 = Vector3::new(-3.748244, 0.0, 0.0);
    let a2 = Vector3::new(1.874122, -4.750729, 0.0);
    let a3 = Vector3::new(0.0, 1.5562529999999997, 6.257947999999999);

    let lattice = Matrix3::from_columns(&[a1, a2, a3]);

    let elements = ["Au", "Au", "Se", "Se"];
    let species: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
    let coords = vec![
        Vector3::new(-1.8741, 0.7781, 3.1290),
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(-1.8741, -3.0213, 1.7273),
        Vector3::new(0.0, -0.1732, 4.5307),
    ];

    let s = Structure::new(lattice, species, coords);
    println!["{:#?}", s];
}
