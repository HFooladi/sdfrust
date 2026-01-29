//! Criterion benchmarks for MOL2 parsing
//!
//! Run with: cargo bench --bench mol2_parse_benchmark

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use sdfrust::{parse_mol2_string, parse_mol2_string_multi, Mol2Iterator};
use std::io::BufReader;

/// Simple methane molecule (5 atoms, 4 bonds) in MOL2 format
const METHANE_MOL2: &str = r#"@<TRIPOS>MOLECULE
methane
 5 4 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3       1 MOL       0.0000
      2 H1          0.6289    0.6289    0.6289 H         1 MOL       0.0000
      3 H2         -0.6289   -0.6289    0.6289 H         1 MOL       0.0000
      4 H3         -0.6289    0.6289   -0.6289 H         1 MOL       0.0000
      5 H4          0.6289   -0.6289   -0.6289 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
     3     1     4 1
     4     1     5 1
"#;

/// Benzene molecule (12 atoms, 12 bonds) in MOL2 format
const BENZENE_MOL2: &str = r#"@<TRIPOS>MOLECULE
benzene
 12 12 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          1.2124    0.7000    0.0000 C.ar      1 MOL       0.0000
      2 C2          1.2124   -0.7000    0.0000 C.ar      1 MOL       0.0000
      3 C3          0.0000   -1.4000    0.0000 C.ar      1 MOL       0.0000
      4 C4         -1.2124   -0.7000    0.0000 C.ar      1 MOL       0.0000
      5 C5         -1.2124    0.7000    0.0000 C.ar      1 MOL       0.0000
      6 C6          0.0000    1.4000    0.0000 C.ar      1 MOL       0.0000
      7 H1          2.1560    1.2432    0.0000 H         1 MOL       0.0000
      8 H2          2.1560   -1.2432    0.0000 H         1 MOL       0.0000
      9 H3          0.0000   -2.4864    0.0000 H         1 MOL       0.0000
     10 H4         -2.1560   -1.2432    0.0000 H         1 MOL       0.0000
     11 H5         -2.1560    1.2432    0.0000 H         1 MOL       0.0000
     12 H6          0.0000    2.4864    0.0000 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 ar
     2     2     3 ar
     3     3     4 ar
     4     4     5 ar
     5     5     6 ar
     6     6     1 ar
     7     1     7 1
     8     2     8 1
     9     3     9 1
    10     4    10 1
    11     5    11 1
    12     6    12 1
"#;

/// Generate a synthetic linear chain molecule in MOL2 format
fn generate_synthetic_chain_mol2(num_atoms: usize) -> String {
    let num_bonds = if num_atoms > 0 { num_atoms - 1 } else { 0 };
    let mut mol2 = String::with_capacity(num_atoms * 100 + num_bonds * 40 + 200);

    // MOLECULE section
    mol2.push_str("@<TRIPOS>MOLECULE\n");
    mol2.push_str(&format!("synthetic_chain_{}\n", num_atoms));
    mol2.push_str(&format!(" {} {} 0 0 0\n", num_atoms, num_bonds));
    mol2.push_str("SMALL\n");
    mol2.push_str("NO_CHARGES\n\n");

    // ATOM section
    mol2.push_str("@<TRIPOS>ATOM\n");
    for i in 0..num_atoms {
        let x = i as f64 * 1.5;
        mol2.push_str(&format!(
            "{:7} C{:<10}{:10.4}{:10.4}{:10.4} C.3       1 MOL       0.0000\n",
            i + 1,
            i + 1,
            x,
            0.0,
            0.0
        ));
    }

    // BOND section
    mol2.push_str("@<TRIPOS>BOND\n");
    for i in 0..num_bonds {
        mol2.push_str(&format!("{:6}{:6}{:6} 1\n", i + 1, i + 1, i + 2));
    }

    mol2
}

/// Generate a multi-molecule MOL2 string
fn generate_multi_mol2(base_mol: &str, count: usize) -> String {
    base_mol.repeat(count)
}

/// Benchmark parsing single MOL2 molecules of various sizes
fn bench_mol2_parse_single(c: &mut Criterion) {
    let mut group = c.benchmark_group("mol2_parse_single");

    // Small molecule (5 atoms)
    group.throughput(Throughput::Elements(1));
    group.bench_function("methane_5_atoms", |b| {
        b.iter(|| parse_mol2_string(black_box(METHANE_MOL2)).unwrap())
    });

    // Benzene (12 atoms with aromatics)
    group.bench_function("benzene_12_atoms", |b| {
        b.iter(|| parse_mol2_string(black_box(BENZENE_MOL2)).unwrap())
    });

    // Synthetic molecules of various sizes
    let sizes = [10, 50, 100, 500];
    for size in sizes {
        let mol2 = generate_synthetic_chain_mol2(size);
        group.bench_with_input(
            BenchmarkId::new("synthetic_chain", size),
            &mol2,
            |b, mol2| {
                b.iter(|| parse_mol2_string(black_box(mol2)).unwrap())
            },
        );
    }

    group.finish();
}

/// Benchmark parsing multi-molecule MOL2 files
fn bench_mol2_parse_multi(c: &mut Criterion) {
    let mut group = c.benchmark_group("mol2_parse_multi");

    let counts = [10, 100, 1000, 10000];
    for count in counts {
        let multi_mol2 = generate_multi_mol2(METHANE_MOL2, count);
        group.throughput(Throughput::Elements(count as u64));
        group.bench_with_input(
            BenchmarkId::new("methane_molecules", count),
            &multi_mol2,
            |b, mol2| {
                b.iter(|| parse_mol2_string_multi(black_box(mol2)).unwrap())
            },
        );
    }

    group.finish();
}

/// Benchmark streaming iterator parsing for MOL2
fn bench_mol2_parse_iterator(c: &mut Criterion) {
    let mut group = c.benchmark_group("mol2_parse_iterator");

    let counts = [100, 1000, 10000];
    for count in counts {
        let multi_mol2 = generate_multi_mol2(METHANE_MOL2, count);
        group.throughput(Throughput::Elements(count as u64));
        group.bench_with_input(
            BenchmarkId::new("iterator_methane", count),
            &multi_mol2,
            |b, mol2| {
                b.iter(|| {
                    let reader = BufReader::new(mol2.as_bytes());
                    let iter = Mol2Iterator::new(reader);
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

/// Benchmark parsing from real MOL2 test files
fn bench_mol2_parse_real_files(c: &mut Criterion) {
    let mut group = c.benchmark_group("mol2_parse_real_files");

    let test_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/test_data");

    let files = [("methane", "methane.mol2"), ("benzene", "benzene.mol2")];

    for (name, filename) in files {
        let filepath = test_dir.join(filename);
        if filepath.exists() {
            let content = std::fs::read_to_string(&filepath).unwrap();
            group.bench_with_input(BenchmarkId::new("file", name), &content, |b, content| {
                b.iter(|| parse_mol2_string(black_box(content)).unwrap())
            });
        }
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_mol2_parse_single,
    bench_mol2_parse_multi,
    bench_mol2_parse_iterator,
    bench_mol2_parse_real_files,
);
criterion_main!(benches);
