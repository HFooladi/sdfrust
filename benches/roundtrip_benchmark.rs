//! Criterion benchmarks for parse + write round-trip operations
//!
//! Run with: cargo bench --bench roundtrip_benchmark

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use sdfrust::{parse_sdf_string, parse_sdf_string_multi, write_sdf_string, Atom, Bond, BondOrder, Molecule};

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

/// Create a linear chain molecule directly (without parsing)
fn create_chain_molecule(num_atoms: usize) -> Molecule {
    let name = format!("synthetic_chain_{}", num_atoms);
    let mut mol = Molecule::new(&name);
    for i in 0..num_atoms {
        mol.atoms.push(Atom::new(i, "C", i as f64 * 1.5, 0.0, 0.0));
    }
    for i in 0..num_atoms.saturating_sub(1) {
        mol.bonds.push(Bond::new(i, i + 1, BondOrder::Single));
    }
    mol
}

/// Generate a synthetic linear chain molecule in SDF format
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

/// Benchmark parse -> write round-trip for single molecules
fn bench_roundtrip_single(c: &mut Criterion) {
    let mut group = c.benchmark_group("roundtrip_single");

    // Small molecule (5 atoms)
    group.throughput(Throughput::Elements(1));
    group.bench_function("methane_5_atoms", |b| {
        b.iter(|| {
            let mol = parse_sdf_string(black_box(METHANE_SDF)).unwrap();
            write_sdf_string(black_box(&mol)).unwrap()
        })
    });

    // Medium molecule (21 atoms)
    group.bench_function("aspirin_21_atoms", |b| {
        b.iter(|| {
            let mol = parse_sdf_string(black_box(ASPIRIN_SDF)).unwrap();
            write_sdf_string(black_box(&mol)).unwrap()
        })
    });

    // Various synthetic sizes (max 999 for V2000 format)
    let sizes = [10, 50, 100, 500];
    for size in sizes {
        let sdf = generate_synthetic_chain(size);
        group.bench_with_input(
            BenchmarkId::new("synthetic_chain", size),
            &sdf,
            |b, sdf| {
                b.iter(|| {
                    let mol = parse_sdf_string(black_box(sdf)).unwrap();
                    write_sdf_string(black_box(&mol)).unwrap()
                })
            },
        );
    }

    group.finish();
}

/// Benchmark write-only (pre-built molecules)
fn bench_write_only(c: &mut Criterion) {
    let mut group = c.benchmark_group("write_only");

    let sizes = [5, 21, 50, 100, 500];
    for size in sizes {
        let mol = create_chain_molecule(size);
        group.bench_with_input(BenchmarkId::new("write_molecule", size), &mol, |b, mol| {
            b.iter(|| write_sdf_string(black_box(mol)).unwrap())
        });
    }

    group.finish();
}

/// Benchmark roundtrip with multi-molecule files
fn bench_roundtrip_multi(c: &mut Criterion) {
    let mut group = c.benchmark_group("roundtrip_multi");

    let counts = [10, 100, 1000];
    for count in counts {
        let multi_sdf = METHANE_SDF.repeat(count);
        group.throughput(Throughput::Elements(count as u64));
        group.bench_with_input(
            BenchmarkId::new("roundtrip_batch", count),
            &multi_sdf,
            |b, sdf| {
                b.iter(|| {
                    let mols = parse_sdf_string_multi(black_box(sdf)).unwrap();
                    let mut output = String::new();
                    for mol in &mols {
                        output.push_str(&write_sdf_string(mol).unwrap());
                    }
                    output
                })
            },
        );
    }

    group.finish();
}

/// Benchmark molecule operations (not just parsing)
fn bench_molecule_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("molecule_operations");

    // Pre-build molecules of various sizes
    let sizes = [10, 50, 100, 500];
    for size in sizes {
        let mol = create_chain_molecule(size);

        // Formula calculation
        group.bench_with_input(BenchmarkId::new("formula", size), &mol, |b, mol| {
            b.iter(|| mol.formula())
        });

        // Centroid calculation
        group.bench_with_input(BenchmarkId::new("centroid", size), &mol, |b, mol| {
            b.iter(|| mol.centroid())
        });

        // Element counts
        group.bench_with_input(BenchmarkId::new("element_counts", size), &mol, |b, mol| {
            b.iter(|| mol.element_counts())
        });

        // Neighbor lookup (middle atom)
        let middle = size / 2;
        group.bench_with_input(BenchmarkId::new("neighbors", size), &mol, |b, mol| {
            b.iter(|| mol.neighbors(black_box(middle)))
        });
    }

    group.finish();
}

/// Benchmark property operations
fn bench_property_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("property_operations");

    // Create molecule with many properties
    let mut mol_with_props = Molecule::new("test");
    mol_with_props.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    for i in 0..100 {
        mol_with_props.set_property(&format!("PROP_{}", i), &format!("value_{}", i));
    }

    // Property lookup
    group.bench_function("property_get", |b| {
        b.iter(|| mol_with_props.get_property(black_box("PROP_50")))
    });

    // Property set (same key - update)
    group.bench_function("property_set_update", |b| {
        let mut mol = mol_with_props.clone();
        b.iter(|| mol.set_property(black_box("PROP_50"), black_box("new_value")))
    });

    // Property set (new key)
    group.bench_function("property_set_new", |b| {
        let mut mol = mol_with_props.clone();
        let mut counter = 0;
        b.iter(|| {
            mol.set_property(&format!("NEW_{}", counter), "value");
            counter += 1;
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_roundtrip_single,
    bench_write_only,
    bench_roundtrip_multi,
    bench_molecule_operations,
    bench_property_operations,
);
criterion_main!(benches);
