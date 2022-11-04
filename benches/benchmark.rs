use criterion::{black_box, criterion_group, criterion_main, Criterion, SamplingMode};
use crystsymm::{analyzer::SymmetryAnalyzer, structure::Structure};
use nalgebra::{Matrix3, Vector3};

fn create_structure(scaling_factor: f64) -> Structure {
    let a = Vector3::new(-3.748244, 0.0, 0.0);
    let b = Vector3::new(1.874122, -4.750729, 0.0);
    let c = Vector3::new(0.0, 1.5562529, 6.25794799);

    let lattice = Matrix3::from_columns(&[a, b, c]);

    let elements = ["Au", "Au", "Se", "Se"];
    let species: Vec<String> = elements.iter().map(|s| s.to_string()).collect();
    let coords = vec![
        Vector3::new(-1.8741, 0.7781, 3.1290),
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(-1.8741, -3.0213, 1.7273),
        Vector3::new(0.0, -0.1732, 4.5307),
    ];

    let mut s = Structure::new(lattice, species, coords, true);

    let tmat = Matrix3::new(
        scaling_factor,
        0.0,
        0.0,
        0.0,
        scaling_factor,
        0.0,
        0.0,
        0.0,
        scaling_factor,
    );

    s.apply_transformation(&tmat, &0.05);

    s
}

fn criterion_benchmark(c: &mut Criterion) {
    let structure = create_structure(2.0);
    let sa = SymmetryAnalyzer {
        dtol: 0.05,
        atol: 5.0,
    };

    let mut group = c.benchmark_group("analyzer");
    group.sampling_mode(SamplingMode::Flat);

    group.bench_function("spg_determination", |b| {
        b.iter(|| sa.get_space_group_operations(black_box(&structure)))
    });

    group.bench_function("std_conv_determination", |b| {
        b.iter(|| sa.get_standard_conventional_structure(black_box(&structure)))
    });

    group.bench_function("std_prim_determination", |b| {
        b.iter(|| sa.get_standard_primitive_structure(black_box(&structure)))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
