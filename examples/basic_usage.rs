//! Basic usage example for sdfrust
//!
//! This example demonstrates:
//! - Creating molecules programmatically
//! - Adding atoms and bonds
//! - Setting properties
//! - Computing molecular formula
//! - Writing to SDF format
//!
//! Run with: cargo run --example basic_usage

use sdfrust::{Atom, Bond, BondOrder, Molecule, write_sdf_string};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== sdfrust Basic Usage Example ===\n");

    // Create a new molecule (water: H2O)
    let mut water = Molecule::new("water");

    // Add atoms with 3D coordinates
    water.atoms.push(Atom::new(0, "O", 0.0000, 0.0000, 0.1173));
    water.atoms.push(Atom::new(1, "H", 0.7572, 0.0000, -0.4692));
    water.atoms.push(Atom::new(2, "H", -0.7572, 0.0000, -0.4692));

    // Add bonds (O-H single bonds)
    water.bonds.push(Bond::new(0, 1, BondOrder::Single));
    water.bonds.push(Bond::new(0, 2, BondOrder::Single));

    // Add some properties
    water.set_property("MOLECULAR_WEIGHT", "18.015");
    water.set_property("COMMON_NAME", "Water");
    water.set_property("CAS_NUMBER", "7732-18-5");

    // Display molecule info
    println!("Molecule: {}", water.name);
    println!("Atoms: {}", water.atom_count());
    println!("Bonds: {}", water.bond_count());
    println!("Formula: {}", water.formula());
    println!("Total charge: {}", water.total_charge());
    println!();

    // Show atom details
    println!("Atom details:");
    for atom in water.atoms() {
        println!(
            "  {} ({:>2}): ({:7.4}, {:7.4}, {:7.4})",
            atom.index, atom.element, atom.x, atom.y, atom.z
        );
    }
    println!();

    // Show bond details
    println!("Bond details:");
    for bond in water.bonds() {
        println!(
            "  {} -- {} : {:?}",
            bond.atom1, bond.atom2, bond.order
        );
    }
    println!();

    // Show properties
    println!("Properties:");
    for (key, value) in &water.properties {
        println!("  {}: {}", key, value);
    }
    println!();

    // Compute centroid
    if let Some((cx, cy, cz)) = water.centroid() {
        println!("Centroid: ({:.4}, {:.4}, {:.4})", cx, cy, cz);
    }

    // Get neighbors of oxygen
    let oxygen_neighbors = water.neighbors(0);
    println!("Oxygen (atom 0) is bonded to atoms: {:?}", oxygen_neighbors);
    println!();

    // Center the molecule at origin
    water.center();
    println!("After centering:");
    if let Some((cx, cy, cz)) = water.centroid() {
        println!("  New centroid: ({:.4}, {:.4}, {:.4})", cx, cy, cz);
    }
    println!();

    // Write to SDF format
    let sdf_output = write_sdf_string(&water)?;
    println!("SDF Output:");
    println!("{}", sdf_output);

    // Create a more complex molecule (ethanol: C2H5OH)
    println!("\n=== Creating Ethanol ===\n");

    let mut ethanol = Molecule::new("ethanol");

    // Add heavy atoms first, then hydrogens
    ethanol.atoms.push(Atom::new(0, "C", -0.001, 1.086, 0.008));    // C1
    ethanol.atoms.push(Atom::new(1, "C", 0.002, -0.422, 0.002));    // C2
    ethanol.atoms.push(Atom::new(2, "O", 1.210, -0.907, -0.003));   // O
    ethanol.atoms.push(Atom::new(3, "H", 1.020, 1.465, 0.001));     // H on C1
    ethanol.atoms.push(Atom::new(4, "H", -0.536, 1.449, -0.873));   // H on C1
    ethanol.atoms.push(Atom::new(5, "H", -0.523, 1.437, 0.902));    // H on C1
    ethanol.atoms.push(Atom::new(6, "H", -0.524, -0.785, 0.891));   // H on C2
    ethanol.atoms.push(Atom::new(7, "H", -0.510, -0.793, -0.884));  // H on C2
    ethanol.atoms.push(Atom::new(8, "H", 1.192, -1.870, -0.003));   // H on O

    // Add bonds
    ethanol.bonds.push(Bond::new(0, 1, BondOrder::Single));  // C-C
    ethanol.bonds.push(Bond::new(1, 2, BondOrder::Single));  // C-O
    ethanol.bonds.push(Bond::new(0, 3, BondOrder::Single));  // C-H
    ethanol.bonds.push(Bond::new(0, 4, BondOrder::Single));  // C-H
    ethanol.bonds.push(Bond::new(0, 5, BondOrder::Single));  // C-H
    ethanol.bonds.push(Bond::new(1, 6, BondOrder::Single));  // C-H
    ethanol.bonds.push(Bond::new(1, 7, BondOrder::Single));  // C-H
    ethanol.bonds.push(Bond::new(2, 8, BondOrder::Single));  // O-H

    ethanol.set_property("MOLECULAR_WEIGHT", "46.07");
    ethanol.set_property("BOILING_POINT", "78.37 C");
    ethanol.set_property("SMILES", "CCO");

    println!("Molecule: {}", ethanol.name);
    println!("Formula: {}", ethanol.formula());
    println!("Atoms: {} ({} heavy atoms)",
        ethanol.atom_count(),
        ethanol.atoms().filter(|a| a.element != "H").count()
    );
    println!("Bonds: {}", ethanol.bond_count());

    // Element counts
    let counts = ethanol.element_counts();
    println!("\nElement counts:");
    for (element, count) in &counts {
        println!("  {}: {}", element, count);
    }

    // Check connectivity
    println!("\nConnectivity:");
    for i in 0..3 {  // Just heavy atoms
        let neighbors = ethanol.neighbors(i);
        let element = &ethanol.atoms[i].element;
        println!("  {} (atom {}): {} bonds to {:?}",
            element, i, neighbors.len(), neighbors);
    }

    println!("\n=== Example Complete ===");

    Ok(())
}
