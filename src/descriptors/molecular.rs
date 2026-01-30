//! Molecular property calculations.
//!
//! This module provides functions for calculating molecular properties
//! such as molecular weight, exact mass, and atom/bond counts.

use std::collections::HashMap;

use crate::bond::BondOrder;
use crate::molecule::Molecule;

use super::elements;

/// Calculate the molecular weight (sum of atomic weights).
///
/// Uses standard atomic weights (IUPAC 2021) for each element.
/// Returns `None` if any atom has an unknown element.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, descriptors};
///
/// let mut mol = Molecule::new("water");
/// mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "H", -0.3, 0.95, 0.0));
///
/// let mw = descriptors::molecular_weight(&mol).unwrap();
/// assert!((mw - 18.015).abs() < 0.01);  // H2O ≈ 18.015 Da
/// ```
pub fn molecular_weight(mol: &Molecule) -> Option<f64> {
    let mut total = 0.0;
    for atom in &mol.atoms {
        let weight = elements::atomic_weight(&atom.element)?;
        total += weight;
    }
    Some(total)
}

/// Calculate the exact (monoisotopic) mass.
///
/// Uses the mass of the most abundant isotope for each element.
/// Returns `None` if any atom has an unknown element.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, descriptors};
///
/// let mut mol = Molecule::new("methane");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "H", -1.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(3, "H", 0.0, 1.0, 0.0));
/// mol.atoms.push(Atom::new(4, "H", 0.0, -1.0, 0.0));
///
/// let mass = descriptors::exact_mass(&mol).unwrap();
/// // CH4: 12.0 + 4*1.00783 ≈ 16.031
/// assert!((mass - 16.031).abs() < 0.01);
/// ```
pub fn exact_mass(mol: &Molecule) -> Option<f64> {
    let mut total = 0.0;
    for atom in &mol.atoms {
        let mass = elements::monoisotopic_mass(&atom.element)?;
        total += mass;
    }
    Some(total)
}

/// Count non-hydrogen atoms (heavy atoms).
///
/// Heavy atoms are all atoms except hydrogen (H), deuterium (D), and tritium (T).
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, descriptors};
///
/// let mut mol = Molecule::new("water");
/// mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "H", -0.3, 0.95, 0.0));
///
/// assert_eq!(descriptors::heavy_atom_count(&mol), 1);  // Only O
/// ```
pub fn heavy_atom_count(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .filter(|a| !is_hydrogen(&a.element))
        .count()
}

/// Count bonds by bond order.
///
/// Returns a HashMap mapping each BondOrder to its count.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder, descriptors};
///
/// let mut mol = Molecule::new("ethene");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "C", 1.3, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "H", -0.5, 0.9, 0.0));
/// mol.atoms.push(Atom::new(3, "H", -0.5, -0.9, 0.0));
/// mol.atoms.push(Atom::new(4, "H", 1.8, 0.9, 0.0));
/// mol.atoms.push(Atom::new(5, "H", 1.8, -0.9, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Double));
/// mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
/// mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
/// mol.bonds.push(Bond::new(1, 4, BondOrder::Single));
/// mol.bonds.push(Bond::new(1, 5, BondOrder::Single));
///
/// let counts = descriptors::bond_type_counts(&mol);
/// assert_eq!(counts.get(&BondOrder::Single), Some(&4));
/// assert_eq!(counts.get(&BondOrder::Double), Some(&1));
/// ```
pub fn bond_type_counts(mol: &Molecule) -> HashMap<BondOrder, usize> {
    let mut counts = HashMap::new();
    for bond in &mol.bonds {
        *counts.entry(bond.order).or_insert(0) += 1;
    }
    counts
}

/// Check if an element symbol represents hydrogen.
fn is_hydrogen(element: &str) -> bool {
    let elem = element.trim().to_uppercase();
    elem == "H" || elem == "D" || elem == "T"
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::Bond;

    fn make_water() -> Molecule {
        let mut mol = Molecule::new("water");
        mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 0.96, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "H", -0.24, 0.93, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol
    }

    fn make_methane() -> Molecule {
        let mut mol = Molecule::new("methane");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 0.63, 0.63, 0.63));
        mol.atoms.push(Atom::new(2, "H", -0.63, -0.63, 0.63));
        mol.atoms.push(Atom::new(3, "H", -0.63, 0.63, -0.63));
        mol.atoms.push(Atom::new(4, "H", 0.63, -0.63, -0.63));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));
        mol
    }

    #[test]
    fn test_molecular_weight_water() {
        let mol = make_water();
        let mw = molecular_weight(&mol).unwrap();
        // H2O: 2*1.008 + 15.999 = 18.015
        assert!((mw - 18.015).abs() < 0.001);
    }

    #[test]
    fn test_molecular_weight_methane() {
        let mol = make_methane();
        let mw = molecular_weight(&mol).unwrap();
        // CH4: 12.011 + 4*1.008 = 16.043
        assert!((mw - 16.043).abs() < 0.001);
    }

    #[test]
    fn test_molecular_weight_unknown_element() {
        let mut mol = Molecule::new("unknown");
        mol.atoms.push(Atom::new(0, "Xx", 0.0, 0.0, 0.0));
        assert!(molecular_weight(&mol).is_none());
    }

    #[test]
    fn test_molecular_weight_empty() {
        let mol = Molecule::new("empty");
        let mw = molecular_weight(&mol).unwrap();
        assert!((mw - 0.0).abs() < 0.001);
    }

    #[test]
    fn test_exact_mass_water() {
        let mol = make_water();
        let mass = exact_mass(&mol).unwrap();
        // H2O: 2*1.00783 + 15.99491 = 18.01056
        assert!((mass - 18.01056).abs() < 0.001);
    }

    #[test]
    fn test_exact_mass_methane() {
        let mol = make_methane();
        let mass = exact_mass(&mol).unwrap();
        // CH4: 12.0 + 4*1.00783 = 16.03130
        assert!((mass - 16.03130).abs() < 0.001);
    }

    #[test]
    fn test_heavy_atom_count_water() {
        let mol = make_water();
        assert_eq!(heavy_atom_count(&mol), 1);
    }

    #[test]
    fn test_heavy_atom_count_methane() {
        let mol = make_methane();
        assert_eq!(heavy_atom_count(&mol), 1);
    }

    #[test]
    fn test_heavy_atom_count_empty() {
        let mol = Molecule::new("empty");
        assert_eq!(heavy_atom_count(&mol), 0);
    }

    #[test]
    fn test_heavy_atom_count_all_heavy() {
        let mut mol = Molecule::new("co2");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "O", -1.2, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "O", 1.2, 0.0, 0.0));
        assert_eq!(heavy_atom_count(&mol), 3);
    }

    #[test]
    fn test_bond_type_counts_water() {
        let mol = make_water();
        let counts = bond_type_counts(&mol);
        assert_eq!(counts.get(&BondOrder::Single), Some(&2));
        assert_eq!(counts.get(&BondOrder::Double), None);
    }

    #[test]
    fn test_bond_type_counts_mixed() {
        let mut mol = Molecule::new("formaldehyde");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "O", 1.2, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "H", -0.5, 0.9, 0.0));
        mol.atoms.push(Atom::new(3, "H", -0.5, -0.9, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));

        let counts = bond_type_counts(&mol);
        assert_eq!(counts.get(&BondOrder::Single), Some(&2));
        assert_eq!(counts.get(&BondOrder::Double), Some(&1));
    }

    #[test]
    fn test_bond_type_counts_empty() {
        let mol = Molecule::new("empty");
        let counts = bond_type_counts(&mol);
        assert!(counts.is_empty());
    }

    #[test]
    fn test_is_hydrogen() {
        assert!(is_hydrogen("H"));
        assert!(is_hydrogen("h"));
        assert!(is_hydrogen("D"));
        assert!(is_hydrogen("T"));
        assert!(!is_hydrogen("C"));
        assert!(!is_hydrogen("He"));
    }
}
