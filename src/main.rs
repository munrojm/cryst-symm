use crystsymm::Structure;

fn main() {
    let lattice = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    let elements = ["Si", "Si", "Si", "Si", "O", "O"];
    let species: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
    let coords = vec![vec![0.0, 0.0, 0.5]];

    let s = Structure::new(lattice, species, coords);
    println!["{}", s]
}
