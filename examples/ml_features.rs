//! ML-ready feature computation example for sdfrust
//!
//! Demonstrates:
//! - OGB-compatible GNN atom and bond featurization
//! - ECFP/Morgan fingerprints and Tanimoto similarity
//! - Gasteiger partial charges
//! - Ring perception (SSSR), aromaticity, hybridization, conjugation
//! - Valence and hydrogen count
//!
//! Run with: cargo run --example ml_features

use sdfrust::descriptors::{
    all_aromatic_atoms, all_aromatic_bonds, all_conjugated_bonds, all_hybridizations,
    all_total_hydrogen_counts, gasteiger_charges, sssr,
};
use sdfrust::featurize::ogb;
use sdfrust::fingerprints::ecfp;
use sdfrust::graph::AdjacencyList;
use sdfrust::{Atom, Bond, BondOrder, Molecule};

fn make_benzene() -> Molecule {
    let mut mol = Molecule::new("benzene");
    let coords: [(f64, f64); 6] = [
        (1.21, 0.70),
        (1.21, -0.70),
        (0.00, -1.40),
        (-1.21, -0.70),
        (-1.21, 0.70),
        (0.00, 1.40),
    ];
    for (i, &(x, y)) in coords.iter().enumerate() {
        mol.atoms.push(Atom::new(i, "C", x, y, 0.0));
    }
    for i in 0..6 {
        mol.bonds
            .push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
    }
    mol
}

fn make_ethanol() -> Molecule {
    let mut mol = Molecule::new("ethanol");
    mol.atoms.push(Atom::new(0, "C", -0.001, 1.086, 0.008));
    mol.atoms.push(Atom::new(1, "C", 0.002, -0.422, 0.002));
    mol.atoms.push(Atom::new(2, "O", 1.210, -0.907, -0.003));
    mol.atoms.push(Atom::new(3, "H", 1.020, 1.465, 0.001));
    mol.atoms.push(Atom::new(4, "H", -0.536, 1.449, -0.873));
    mol.atoms.push(Atom::new(5, "H", -0.523, 1.437, 0.902));
    mol.atoms.push(Atom::new(6, "H", -0.524, -0.785, 0.891));
    mol.atoms.push(Atom::new(7, "H", -0.510, -0.793, -0.884));
    mol.atoms.push(Atom::new(8, "H", 1.192, -1.870, -0.003));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
    mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
    mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
    mol.bonds.push(Bond::new(0, 4, BondOrder::Single));
    mol.bonds.push(Bond::new(0, 5, BondOrder::Single));
    mol.bonds.push(Bond::new(1, 6, BondOrder::Single));
    mol.bonds.push(Bond::new(1, 7, BondOrder::Single));
    mol.bonds.push(Bond::new(2, 8, BondOrder::Single));
    mol
}

fn make_aspirin() -> Molecule {
    // Simplified aspirin (just the ring + functional groups, no H)
    let mut mol = Molecule::new("aspirin");
    // Benzene ring: C0-C5
    for i in 0..6 {
        mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
    }
    for i in 0..6 {
        let order = if i % 2 == 0 {
            BondOrder::Double
        } else {
            BondOrder::Single
        };
        mol.bonds.push(Bond::new(i, (i + 1) % 6, order));
    }
    // Carboxylic acid: C6(=O7)-O8
    mol.atoms.push(Atom::new(6, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(7, "O", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(8, "O", 0.0, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 6, BondOrder::Single));
    mol.bonds.push(Bond::new(6, 7, BondOrder::Double));
    mol.bonds.push(Bond::new(6, 8, BondOrder::Single));
    // Ester: O9-C10(=O11)-C12
    mol.atoms.push(Atom::new(9, "O", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(10, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(11, "O", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(12, "C", 0.0, 0.0, 0.0));
    mol.bonds.push(Bond::new(3, 9, BondOrder::Single));
    mol.bonds.push(Bond::new(9, 10, BondOrder::Single));
    mol.bonds.push(Bond::new(10, 11, BondOrder::Double));
    mol.bonds.push(Bond::new(10, 12, BondOrder::Single));
    mol
}

fn main() {
    println!("=== sdfrust ML Features Example ===\n");

    // --- 1. Ring Perception & Aromaticity ---
    println!("--- Ring Perception & Aromaticity ---\n");
    let benzene = make_benzene();
    let rings = sssr(&benzene);
    println!("Benzene SSSR: {} ring(s)", rings.len());
    for (i, ring) in rings.iter().enumerate() {
        println!(
            "  Ring {}: atoms {:?} (size {})",
            i,
            ring.atoms,
            ring.atoms.len()
        );
    }

    let aromatic_atoms = all_aromatic_atoms(&benzene);
    println!("Aromatic atoms: {:?}", aromatic_atoms);

    let aromatic_bonds = all_aromatic_bonds(&benzene);
    println!("Aromatic bonds: {:?}", aromatic_bonds);
    println!();

    // --- 2. Hybridization ---
    println!("--- Hybridization ---\n");
    let ethanol = make_ethanol();
    let hybs = all_hybridizations(&ethanol);
    for (i, hyb) in hybs.iter().enumerate() {
        let elem = &ethanol.atoms[i].element;
        println!("  Atom {} ({}): {:?}", i, elem, hyb);
    }
    println!();

    // --- 3. Conjugation ---
    println!("--- Conjugation ---\n");
    let aspirin = make_aspirin();
    let conjugated = all_conjugated_bonds(&aspirin);
    for (i, &conj) in conjugated.iter().enumerate() {
        let bond = &aspirin.bonds[i];
        println!(
            "  Bond {} ({}-{}): conjugated={}",
            i, bond.atom1, bond.atom2, conj
        );
    }
    println!();

    // --- 4. Hydrogen Counts ---
    println!("--- Hydrogen Counts ---\n");
    let h_counts = all_total_hydrogen_counts(&ethanol);
    for (i, &count) in h_counts.iter().enumerate() {
        let elem = &ethanol.atoms[i].element;
        println!("  Atom {} ({}): {} total H", i, elem, count);
    }
    println!();

    // --- 5. OGB GNN Features ---
    println!("--- OGB GNN Features ---\n");
    let atom_feats = ogb::ogb_atom_features(&benzene);
    println!("Benzene atom features [N={}, F=9]:", atom_feats.num_atoms);
    for (i, feat) in atom_feats.features.iter().enumerate() {
        println!("  Atom {}: {:?}", i, feat);
    }

    let bond_feats = ogb::ogb_bond_features(&benzene);
    println!("\nBenzene bond features [E={}, F=3]:", bond_feats.num_bonds);
    for (i, feat) in bond_feats.features.iter().enumerate() {
        println!("  Bond {}: {:?}", i, feat);
    }

    let graph = ogb::ogb_graph_features(&benzene);
    println!(
        "\nFull graph: {} atoms, {} directed edges",
        graph.atom_features.num_atoms,
        graph.edge_src.len()
    );
    println!();

    // --- 6. ECFP Fingerprints ---
    println!("--- ECFP Fingerprints ---\n");
    let fp_benzene = ecfp::ecfp(&benzene, 2, 2048);
    let fp_ethanol = ecfp::ecfp(&ethanol, 2, 2048);
    println!(
        "Benzene ECFP4: {} bits set / {}",
        fp_benzene.num_on_bits(),
        fp_benzene.n_bits
    );
    println!(
        "Ethanol ECFP4: {} bits set / {}",
        fp_ethanol.num_on_bits(),
        fp_ethanol.n_bits
    );
    println!("Benzene density: {:.4}", fp_benzene.density());
    println!(
        "Tanimoto(benzene, ethanol): {:.4}",
        fp_benzene.tanimoto(&fp_ethanol)
    );
    println!(
        "Tanimoto(benzene, benzene): {:.4}",
        fp_benzene.tanimoto(&fp_benzene)
    );

    // Count fingerprint
    let counts = ecfp::ecfp_counts(&benzene, 2);
    println!(
        "\nBenzene ECFP4 count fingerprint: {} unique features",
        counts.counts.len()
    );
    println!();

    // --- 7. Gasteiger Charges ---
    println!("--- Gasteiger Charges ---\n");
    let charges = gasteiger_charges(&ethanol);
    for (i, &charge) in charges.iter().enumerate() {
        let elem = &ethanol.atoms[i].element;
        println!("  Atom {} ({}): charge = {:+.4}", i, elem, charge);
    }
    println!();

    // --- 8. Graph Adjacency ---
    println!("--- Graph Adjacency ---\n");
    let adj = AdjacencyList::from_molecule(&ethanol);
    for i in 0..3 {
        let elem = &ethanol.atoms[i].element;
        let neighbors = adj.neighbor_atoms(i);
        println!(
            "  Atom {} ({}): degree={}, heavy_degree={}, neighbors={:?}",
            i,
            elem,
            adj.degree(i),
            adj.heavy_degree(i),
            neighbors
        );
    }

    println!("\n=== Example Complete ===");
}
