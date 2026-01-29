//! Example: Parsing SDF files and analyzing molecules
//!
//! This example demonstrates:
//! - Parsing SDF files from disk
//! - Iterating over multi-molecule files
//! - Accessing molecular properties
//! - Computing descriptors
//! - Round-trip parsing and writing
//!
//! Run with: cargo run --example parse_molecules

use sdfrust::{BondOrder, parse_sdf_file, parse_sdf_string, write_sdf_string};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== sdfrust Parsing Example ===\n");

    // Parse a single molecule from file
    println!("--- Parsing Aspirin from PubChem ---");
    let aspirin = parse_sdf_file("tests/test_data/aspirin.sdf")?;

    println!("Name: {}", aspirin.name);
    println!("Atoms: {}", aspirin.atom_count());
    println!("Bonds: {}", aspirin.bond_count());
    println!("Formula: {}", aspirin.formula());

    // Access PubChem properties
    if let Some(cid) = aspirin.get_property("PUBCHEM_COMPOUND_CID") {
        println!("PubChem CID: {}", cid);
    }
    if let Some(mw) = aspirin.get_property("PUBCHEM_MOLECULAR_WEIGHT") {
        println!("Molecular Weight: {} g/mol", mw);
    }
    if let Some(smiles) = aspirin.get_property("PUBCHEM_SMILES") {
        println!("SMILES: {}", smiles);
    }

    // Analyze bond types
    let single = aspirin
        .bonds()
        .filter(|b| b.order == BondOrder::Single)
        .count();
    let double = aspirin
        .bonds()
        .filter(|b| b.order == BondOrder::Double)
        .count();
    let aromatic = aspirin
        .bonds()
        .filter(|b| b.order == BondOrder::Aromatic)
        .count();
    println!("\nBond distribution:");
    println!("  Single: {}", single);
    println!("  Double: {}", double);
    println!("  Aromatic: {}", aromatic);

    // Heavy atom analysis
    let heavy_atoms: Vec<_> = aspirin.atoms().filter(|a| a.element != "H").collect();
    println!("\nHeavy atoms ({}):", heavy_atoms.len());
    for atom in &heavy_atoms {
        let bonds = aspirin.bonds_for_atom(atom.index);
        println!(
            "  {:>2} {}: {} bonds",
            atom.index,
            atom.element,
            bonds.len()
        );
    }

    println!("\n--- Geometry Analysis ---");

    // Compute centroid
    if let Some((cx, cy, cz)) = aspirin.centroid() {
        println!("Centroid: ({:.3}, {:.3}, {:.3})", cx, cy, cz);
    }

    // Find atom distances
    let c1 = &aspirin.atoms[4]; // First carbon
    let c2 = &aspirin.atoms[5]; // Second carbon
    let distance = c1.distance_to(c2);
    println!("C-C bond length: {:.3} Angstroms", distance);

    println!("\n--- Parsing Caffeine ---");
    let caffeine = parse_sdf_file("tests/test_data/caffeine_pubchem.sdf")?;

    println!("Name: {}", caffeine.name);
    println!("Formula: {}", caffeine.formula());

    // Element breakdown
    let counts = caffeine.element_counts();
    print!("Elements: ");
    let mut elements: Vec<_> = counts.iter().collect();
    elements.sort_by_key(|(e, _)| *e);
    for (element, count) in elements {
        print!("{}:{} ", element, count);
    }
    println!();

    // Check for nitrogen (caffeine has 4)
    let nitrogens = caffeine.atoms_by_element("N");
    println!(
        "Nitrogen atoms: {} (caffeine has nitrogen heterocycles)",
        nitrogens.len()
    );

    println!("\n--- Multi-molecule Iteration ---");

    // Create a combined SDF string for demonstration
    let combined_sdf = format!(
        "{}{}",
        std::fs::read_to_string("tests/test_data/aspirin.sdf")?,
        std::fs::read_to_string("tests/test_data/caffeine_pubchem.sdf")?
    );

    // Parse from string (multi-molecule)
    let cursor = std::io::Cursor::new(combined_sdf.as_bytes());
    let reader = std::io::BufReader::new(cursor);

    println!("Iterating over molecules in combined SDF:");
    for (i, result) in sdfrust::SdfIterator::new(reader).enumerate() {
        let mol = result?;
        println!(
            "  [{}] {} - {} atoms, {} bonds, formula: {}",
            i + 1,
            mol.name,
            mol.atom_count(),
            mol.bond_count(),
            mol.formula()
        );
    }

    println!("\n--- Iterating Over Directory ---");

    // Process all test files
    let test_files = [
        "tests/test_data/aspirin.sdf",
        "tests/test_data/caffeine_pubchem.sdf",
        "tests/test_data/glucose.sdf",
        "tests/test_data/galactose.sdf",
        "tests/test_data/acetaminophen.sdf",
        "tests/test_data/methionine.sdf",
    ];

    println!(
        "{:<12} {:>6} {:>6} {:>15}",
        "Name", "Atoms", "Bonds", "Formula"
    );
    println!("{}", "-".repeat(45));

    for file in &test_files {
        if let Ok(mol) = parse_sdf_file(file) {
            println!(
                "{:<12} {:>6} {:>6} {:>15}",
                &mol.name[..mol.name.len().min(12)],
                mol.atom_count(),
                mol.bond_count(),
                mol.formula()
            );
        }
    }

    println!("\n--- Round-trip Test ---");

    // Demonstrate round-trip: parse -> write -> parse
    let original = parse_sdf_file("tests/test_data/aspirin.sdf")?;
    let sdf_string = write_sdf_string(&original)?;
    let reparsed = parse_sdf_string(&sdf_string)?;

    println!("Original atoms: {}", original.atom_count());
    println!("Reparsed atoms: {}", reparsed.atom_count());
    println!(
        "Atoms match: {}",
        original.atom_count() == reparsed.atom_count()
    );

    // Compare coordinates
    let mut coord_match = true;
    for i in 0..original.atom_count() {
        let o = &original.atoms[i];
        let r = &reparsed.atoms[i];
        if (o.x - r.x).abs() > 0.001 || (o.y - r.y).abs() > 0.001 || (o.z - r.z).abs() > 0.001 {
            coord_match = false;
            break;
        }
    }
    println!("Coordinates preserved: {}", coord_match);

    println!("\n=== Example Complete ===");

    Ok(())
}
