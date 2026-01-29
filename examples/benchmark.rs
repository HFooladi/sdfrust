//! Performance benchmark for sdfrust
//!
//! This example measures parsing and writing performance.
//!
//! Run with: cargo run --release --example benchmark

use sdfrust::{Atom, Bond, BondOrder, Molecule, parse_sdf_string, write_sdf_string};
use std::time::Instant;

fn main() {
    println!("=== sdfrust Performance Benchmark ===\n");

    // Benchmark 1: Parse simple molecules
    benchmark_parse_simple();

    // Benchmark 2: Parse and write round-trip
    benchmark_roundtrip();

    // Benchmark 3: Large molecule operations
    benchmark_large_molecule();

    // Benchmark 4: Multi-molecule file
    benchmark_multi_molecule();

    // Benchmark 5: Property access
    benchmark_properties();

    println!("\n=== Benchmark Complete ===");
}

fn benchmark_parse_simple() {
    println!("--- Benchmark: Parse Simple Molecules ---");

    let methane = r#"methane


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

    let iterations = 10000;
    let start = Instant::now();

    for _ in 0..iterations {
        let _ = parse_sdf_string(methane).unwrap();
    }

    let elapsed = start.elapsed();
    let per_mol = elapsed.as_nanos() as f64 / iterations as f64;

    println!("  Parsed {} molecules in {:?}", iterations, elapsed);
    println!("  Average: {:.2} ns/molecule", per_mol);
    println!("  Rate: {:.0} molecules/sec", 1e9 / per_mol);
}

fn benchmark_roundtrip() {
    println!("\n--- Benchmark: Parse + Write Round-trip ---");

    let aspirin_like = r#"test_molecule


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

    let iterations = 5000;
    let start = Instant::now();

    for _ in 0..iterations {
        let mol = parse_sdf_string(aspirin_like).unwrap();
        let _ = write_sdf_string(&mol).unwrap();
    }

    let elapsed = start.elapsed();
    let per_mol = elapsed.as_nanos() as f64 / iterations as f64;

    println!("  Round-tripped {} molecules in {:?}", iterations, elapsed);
    println!("  Average: {:.2} ns/molecule", per_mol);
    println!("  Rate: {:.0} round-trips/sec", 1e9 / per_mol);
}

fn benchmark_large_molecule() {
    println!("\n--- Benchmark: Large Molecule Operations ---");

    // Create a large molecule (100 atoms, 99 bonds - linear chain)
    let mut mol = Molecule::new("large_chain");
    for i in 0..100 {
        mol.atoms.push(Atom::new(i, "C", i as f64 * 1.5, 0.0, 0.0));
    }
    for i in 0..99 {
        mol.bonds.push(Bond::new(i, i + 1, BondOrder::Single));
    }

    let iterations = 10000;

    // Benchmark formula calculation
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = mol.formula();
    }
    let elapsed = start.elapsed();
    println!(
        "  Formula calculation: {:?} for {} iterations ({:.0} ns/call)",
        elapsed,
        iterations,
        elapsed.as_nanos() as f64 / iterations as f64
    );

    // Benchmark centroid calculation
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = mol.centroid();
    }
    let elapsed = start.elapsed();
    println!(
        "  Centroid calculation: {:?} for {} iterations ({:.0} ns/call)",
        elapsed,
        iterations,
        elapsed.as_nanos() as f64 / iterations as f64
    );

    // Benchmark neighbor lookup
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = mol.neighbors(50); // Middle atom
    }
    let elapsed = start.elapsed();
    println!(
        "  Neighbor lookup: {:?} for {} iterations ({:.0} ns/call)",
        elapsed,
        iterations,
        elapsed.as_nanos() as f64 / iterations as f64
    );

    // Benchmark element counts
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = mol.element_counts();
    }
    let elapsed = start.elapsed();
    println!(
        "  Element counts: {:?} for {} iterations ({:.0} ns/call)",
        elapsed,
        iterations,
        elapsed.as_nanos() as f64 / iterations as f64
    );
}

fn benchmark_multi_molecule() {
    println!("\n--- Benchmark: Multi-molecule File ---");

    // Generate a multi-molecule SDF string
    let single_mol = r#"mol


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

    let num_mols = 1000;
    let multi_sdf: String = single_mol.repeat(num_mols);

    let iterations = 100;
    let start = Instant::now();

    for _ in 0..iterations {
        let mols = sdfrust::parse_sdf_string_multi(&multi_sdf).unwrap();
        assert_eq!(mols.len(), num_mols);
    }

    let elapsed = start.elapsed();
    let total_mols = iterations * num_mols;
    let per_mol = elapsed.as_nanos() as f64 / total_mols as f64;

    println!(
        "  Parsed {} molecules ({}x{}) in {:?}",
        total_mols, iterations, num_mols, elapsed
    );
    println!("  Average: {:.2} ns/molecule", per_mol);
    println!("  Rate: {:.0} molecules/sec", 1e9 / per_mol);
}

fn benchmark_properties() {
    println!("\n--- Benchmark: Property Access ---");

    // Create molecule with many properties
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));

    for i in 0..100 {
        mol.set_property(&format!("PROP_{}", i), &format!("value_{}", i));
    }

    let iterations = 100000;

    // Benchmark property lookup
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = mol.get_property("PROP_50");
    }
    let elapsed = start.elapsed();
    println!(
        "  Property lookup: {:?} for {} iterations ({:.0} ns/call)",
        elapsed,
        iterations,
        elapsed.as_nanos() as f64 / iterations as f64
    );

    // Benchmark property set
    let start = Instant::now();
    for i in 0..iterations {
        mol.set_property("NEW_PROP", &format!("value_{}", i));
    }
    let elapsed = start.elapsed();
    println!(
        "  Property set: {:?} for {} iterations ({:.0} ns/call)",
        elapsed,
        iterations,
        elapsed.as_nanos() as f64 / iterations as f64
    );
}
