//! Criterion benchmarks for SDF parsing
//!
//! Run with: cargo bench --bench sdf_parse_benchmark

use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};
use sdfrust::{SdfIterator, parse_sdf_string, parse_sdf_string_multi};
use std::io::BufReader;

/// Simple methane molecule (5 atoms, 4 bonds)
const METHANE_SDF: &str = r#"methane


  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289    0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6289   -0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6289    0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289   -0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
"#;

/// Medium molecule (~21 atoms) - aspirin-like structure
const ASPIRIN_SDF: &str = r#"aspirin_like


 21 21  0  0  0  0  0  0  0  0999 V2000
    3.7321   -0.0600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.3301    1.4400    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5981    1.4400    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660   -1.5600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5981   -0.5600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4641   -0.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5981   -1.5600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3301   -0.5600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4641   -2.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3301   -1.5600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4641    0.9400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660   -0.5600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000   -0.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0611   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.8671   -0.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.4641   -2.6800    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.8671   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3100    0.4769    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4631    0.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6900   -0.5969    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.3301    2.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  5  1  0  0  0  0
  1 12  1  0  0  0  0
  2 11  1  0  0  0  0
  2 21  1  0  0  0  0
  3 11  2  0  0  0  0
  4 12  2  0  0  0  0
  5  6  1  0  0  0  0
  5  7  2  0  0  0  0
  6  8  2  0  0  0  0
  6 11  1  0  0  0  0
  7  9  1  0  0  0  0
  7 14  1  0  0  0  0
  8 10  1  0  0  0  0
  8 15  1  0  0  0  0
  9 10  2  0  0  0  0
  9 16  1  0  0  0  0
 10 17  1  0  0  0  0
 12 13  1  0  0  0  0
 13 18  1  0  0  0  0
 13 19  1  0  0  0  0
 13 20  1  0  0  0  0
M  END
$$$$
"#;

/// Generate a synthetic linear chain molecule with the specified number of atoms
fn generate_synthetic_chain(num_atoms: usize) -> String {
    let num_bonds = if num_atoms > 0 { num_atoms - 1 } else { 0 };
    let mut sdf = String::with_capacity(num_atoms * 80 + num_bonds * 30 + 200);

    // Header (line 1: name, line 2: program info, line 3: comment, line 4: counts)
    sdf.push_str(&format!("synthetic_chain_{}\n\n\n", num_atoms));
    sdf.push_str(&format!(
        "{:3}{:3}  0  0  0  0  0  0  0  0999 V2000\n",
        num_atoms, num_bonds
    ));

    // Atoms - linear chain of carbons
    for i in 0..num_atoms {
        let x = i as f64 * 1.5;
        sdf.push_str(&format!(
            "{:10.4}{:10.4}{:10.4} C   0  0  0  0  0  0  0  0  0  0  0  0\n",
            x, 0.0, 0.0
        ));
    }

    // Bonds - connect consecutive atoms
    for i in 0..num_bonds {
        sdf.push_str(&format!("{:3}{:3}  1  0  0  0  0\n", i + 1, i + 2));
    }

    sdf.push_str("M  END\n$$$$\n");
    sdf
}

/// Generate a multi-molecule SDF string
fn generate_multi_sdf(base_mol: &str, count: usize) -> String {
    base_mol.repeat(count)
}

/// Benchmark parsing single molecules of various sizes
fn bench_parse_single(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse_single");

    // Small molecule (5 atoms)
    group.throughput(Throughput::Elements(1));
    group.bench_function("methane_5_atoms", |b| {
        b.iter(|| parse_sdf_string(black_box(METHANE_SDF)).unwrap())
    });

    // Medium molecule (21 atoms)
    group.bench_function("aspirin_21_atoms", |b| {
        b.iter(|| parse_sdf_string(black_box(ASPIRIN_SDF)).unwrap())
    });

    // Synthetic molecules of various sizes (max 999 for V2000 format)
    let sizes = [10, 50, 100, 500];
    for size in sizes {
        let sdf = generate_synthetic_chain(size);
        group.bench_with_input(BenchmarkId::new("synthetic_chain", size), &sdf, |b, sdf| {
            b.iter(|| parse_sdf_string(black_box(sdf)).unwrap())
        });
    }

    group.finish();
}

/// Benchmark parsing multi-molecule files
fn bench_parse_multi(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse_multi");

    let counts = [10, 100, 1000, 10000];
    for count in counts {
        let multi_sdf = generate_multi_sdf(METHANE_SDF, count);
        group.throughput(Throughput::Elements(count as u64));
        group.bench_with_input(
            BenchmarkId::new("methane_molecules", count),
            &multi_sdf,
            |b, sdf| b.iter(|| parse_sdf_string_multi(black_box(sdf)).unwrap()),
        );
    }

    group.finish();
}

/// Benchmark streaming iterator parsing
fn bench_parse_iterator(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse_iterator");

    let counts = [100, 1000, 10000];
    for count in counts {
        let multi_sdf = generate_multi_sdf(METHANE_SDF, count);
        group.throughput(Throughput::Elements(count as u64));
        group.bench_with_input(
            BenchmarkId::new("iterator_methane", count),
            &multi_sdf,
            |b, sdf| {
                b.iter(|| {
                    let reader = BufReader::new(sdf.as_bytes());
                    let iter = SdfIterator::new(reader);
                    let mut count = 0;
                    for mol in iter {
                        black_box(mol.unwrap());
                        count += 1;
                    }
                    count
                })
            },
        );
    }

    group.finish();
}

/// Benchmark parsing from file (using test_data files)
fn bench_parse_real_files(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse_real_files");

    // Read test files at benchmark setup time
    let test_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/test_data");

    let files = [
        ("aspirin", "aspirin.sdf"),
        ("caffeine", "caffeine_pubchem.sdf"),
        ("glucose", "glucose.sdf"),
        ("acetaminophen", "acetaminophen.sdf"),
    ];

    for (name, filename) in files {
        let filepath = test_dir.join(filename);
        if filepath.exists() {
            let content = std::fs::read_to_string(&filepath).unwrap();
            group.bench_with_input(BenchmarkId::new("file", name), &content, |b, content| {
                b.iter(|| parse_sdf_string(black_box(content)).unwrap())
            });
        }
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_parse_single,
    bench_parse_multi,
    bench_parse_iterator,
    bench_parse_real_files,
);
criterion_main!(benches);
