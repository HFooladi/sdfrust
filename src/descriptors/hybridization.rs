//! Atom hybridization inference.
//!
//! Determines SP/SP2/SP3/SP3D/SP3D2 hybridization from bond orders
//! and neighbor counts. Prerequisite for aromaticity perception and
//! Gasteiger charges. Provides 1 of 9 OGB atom features.
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder};
//! use sdfrust::descriptors::hybridization::{Hybridization, atom_hybridization};
//!
//! // Ethylene: C=C → SP2
//! let mut mol = Molecule::new("ethylene");
//! mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "C", 1.3, 0.0, 0.0));
//! mol.bonds.push(Bond::new(0, 1, BondOrder::Double));
//!
//! assert_eq!(atom_hybridization(&mol, 0), Hybridization::SP2);
//! ```

use crate::bond::BondOrder;
use crate::graph::{AdjacencyList, is_hydrogen};
use crate::molecule::Molecule;

/// Atom hybridization state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Hybridization {
    /// s hybridization (not common in organic chemistry).
    S,
    /// sp hybridization (linear, e.g., acetylene C, nitrile C).
    SP,
    /// sp2 hybridization (trigonal planar, e.g., ethylene C, aromatic C).
    SP2,
    /// sp3 hybridization (tetrahedral, e.g., methane C).
    SP3,
    /// sp3d hybridization (trigonal bipyramidal, e.g., PCl5).
    SP3D,
    /// sp3d2 hybridization (octahedral, e.g., SF6).
    SP3D2,
    /// Unknown or cannot determine.
    Other,
}

impl Hybridization {
    /// Returns the OGB-compatible integer index (0-5 for S through SP3D2, 6 for Other).
    pub fn to_ogb_index(&self) -> u8 {
        match self {
            Hybridization::S => 0,
            Hybridization::SP => 1,
            Hybridization::SP2 => 2,
            Hybridization::SP3 => 3,
            Hybridization::SP3D => 4,
            Hybridization::SP3D2 => 5,
            Hybridization::Other => 6,
        }
    }
}

/// Determine the hybridization of atom `idx`.
///
/// Inference rules:
/// - Hydrogen atoms → S (single s orbital)
/// - Any triple bond → SP
/// - Any double bond or aromatic bond → SP2
/// - Otherwise → SP3
/// - Hypervalent: 5 bonds → SP3D, 6 bonds → SP3D2
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::hybridization::{Hybridization, atom_hybridization};
///
/// // Acetylene: C≡C → SP
/// let mut mol = Molecule::new("acetylene");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "C", 1.2, 0.0, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Triple));
///
/// assert_eq!(atom_hybridization(&mol, 0), Hybridization::SP);
/// ```
pub fn atom_hybridization(mol: &Molecule, idx: usize) -> Hybridization {
    if idx >= mol.atoms.len() {
        return Hybridization::Other;
    }

    // Hydrogen atoms use a single s orbital
    if is_hydrogen(&mol.atoms[idx].element) {
        return Hybridization::S;
    }

    let adj = AdjacencyList::from_molecule(mol);
    let neighbors = adj.neighbors(idx);
    let num_neighbors = neighbors.len();

    // Collect bond orders for this atom
    let mut has_triple = false;
    let mut has_double = false;
    let mut has_aromatic = false;

    for &(_neighbor, bond_idx) in neighbors {
        if let Some(bond) = mol.bonds.get(bond_idx) {
            match bond.order {
                BondOrder::Triple => has_triple = true,
                BondOrder::Double => has_double = true,
                BondOrder::Aromatic | BondOrder::SingleOrAromatic | BondOrder::DoubleOrAromatic => {
                    has_aromatic = true
                }
                _ => {}
            }
        }
    }

    // Hypervalent atoms
    if num_neighbors >= 6 {
        return Hybridization::SP3D2;
    }
    if num_neighbors >= 5 {
        return Hybridization::SP3D;
    }

    // Standard hybridization from bond orders
    if has_triple {
        Hybridization::SP
    } else if has_double || has_aromatic {
        Hybridization::SP2
    } else if num_neighbors > 0 {
        Hybridization::SP3
    } else {
        // Isolated atom
        Hybridization::S
    }
}

/// Compute hybridization for all atoms in the molecule.
///
/// Uses a two-pass approach:
/// 1. Initial hybridization from bond orders
/// 2. Upgrade O/N/S atoms from SP3 to SP2 if bonded to any SP2 neighbor
///    (lone pair conjugation, matching RDKit behavior)
///
/// Returns a vector of length `mol.atom_count()`.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::hybridization::{Hybridization, all_hybridizations};
///
/// let mut mol = Molecule::new("ethane");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
///
/// let hybs = all_hybridizations(&mol);
/// assert_eq!(hybs[0], Hybridization::SP3);
/// assert_eq!(hybs[1], Hybridization::SP3);
/// ```
pub fn all_hybridizations(mol: &Molecule) -> Vec<Hybridization> {
    let n = mol.atom_count();
    let mut hybs: Vec<Hybridization> = (0..n).map(|i| atom_hybridization(mol, i)).collect();

    // Second pass: upgrade lone-pair atoms (O, N, S) from SP3 to SP2
    // if they are bonded to any SP2 neighbor (lone pair conjugation).
    // This matches RDKit's behavior for esters, carboxylic acids, amides, etc.
    let adj = AdjacencyList::from_molecule(mol);
    for i in 0..n {
        if hybs[i] != Hybridization::SP3 {
            continue;
        }
        if !has_lone_pair(&mol.atoms[i].element) {
            continue;
        }
        // Check if any neighbor is SP2 or SP
        for &(neighbor, _bond_idx) in adj.neighbors(i) {
            if hybs[neighbor] == Hybridization::SP2 || hybs[neighbor] == Hybridization::SP {
                hybs[i] = Hybridization::SP2;
                break;
            }
        }
    }

    hybs
}

/// Check if an element has lone pairs that can participate in conjugation.
fn has_lone_pair(element: &str) -> bool {
    let elem = element.trim();
    let upper: String = elem
        .chars()
        .next()
        .map(|c| c.to_uppercase().collect::<String>())
        .unwrap_or_default()
        + &elem.chars().skip(1).collect::<String>().to_lowercase();

    matches!(upper.as_str(), "O" | "N" | "S" | "Se" | "P")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};

    #[test]
    fn test_sp3_methane() {
        let mut mol = Molecule::new("methane");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "H", -1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(3, "H", 0.0, 1.0, 0.0));
        mol.atoms.push(Atom::new(4, "H", 0.0, -1.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));

        assert_eq!(atom_hybridization(&mol, 0), Hybridization::SP3);
        // Hydrogen atoms should be S hybridization
        assert_eq!(atom_hybridization(&mol, 1), Hybridization::S);
    }

    #[test]
    fn test_sp2_ethylene() {
        let mut mol = Molecule::new("ethylene");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.3, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));

        assert_eq!(atom_hybridization(&mol, 0), Hybridization::SP2);
        assert_eq!(atom_hybridization(&mol, 1), Hybridization::SP2);
    }

    #[test]
    fn test_sp_acetylene() {
        let mut mol = Molecule::new("acetylene");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.2, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Triple));

        assert_eq!(atom_hybridization(&mol, 0), Hybridization::SP);
    }

    #[test]
    fn test_sp2_aromatic() {
        let mut mol = Molecule::new("benzene");
        for i in 0..6 {
            mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
        }
        for i in 0..6 {
            mol.bonds
                .push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
        }

        for i in 0..6 {
            assert_eq!(atom_hybridization(&mol, i), Hybridization::SP2);
        }
    }

    #[test]
    fn test_sp3d_phosphorus() {
        // PCl5: P with 5 bonds
        let mut mol = Molecule::new("PCl5");
        mol.atoms.push(Atom::new(0, "P", 0.0, 0.0, 0.0));
        for i in 1..6 {
            mol.atoms.push(Atom::new(i, "Cl", i as f64, 0.0, 0.0));
            mol.bonds.push(Bond::new(0, i, BondOrder::Single));
        }

        assert_eq!(atom_hybridization(&mol, 0), Hybridization::SP3D);
    }

    #[test]
    fn test_sp3d2_sulfur() {
        // SF6: S with 6 bonds
        let mut mol = Molecule::new("SF6");
        mol.atoms.push(Atom::new(0, "S", 0.0, 0.0, 0.0));
        for i in 1..7 {
            mol.atoms.push(Atom::new(i, "F", i as f64, 0.0, 0.0));
            mol.bonds.push(Bond::new(0, i, BondOrder::Single));
        }

        assert_eq!(atom_hybridization(&mol, 0), Hybridization::SP3D2);
    }

    #[test]
    fn test_isolated_atom() {
        let mut mol = Molecule::new("He");
        mol.atoms.push(Atom::new(0, "He", 0.0, 0.0, 0.0));

        assert_eq!(atom_hybridization(&mol, 0), Hybridization::S);
    }

    #[test]
    fn test_out_of_bounds() {
        let mol = Molecule::new("empty");
        assert_eq!(atom_hybridization(&mol, 0), Hybridization::Other);
    }

    #[test]
    fn test_ogb_index() {
        assert_eq!(Hybridization::S.to_ogb_index(), 0);
        assert_eq!(Hybridization::SP.to_ogb_index(), 1);
        assert_eq!(Hybridization::SP2.to_ogb_index(), 2);
        assert_eq!(Hybridization::SP3.to_ogb_index(), 3);
        assert_eq!(Hybridization::SP3D.to_ogb_index(), 4);
        assert_eq!(Hybridization::SP3D2.to_ogb_index(), 5);
        assert_eq!(Hybridization::Other.to_ogb_index(), 6);
    }

    #[test]
    fn test_all_hybridizations() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.2, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Triple));

        let hybs = all_hybridizations(&mol);
        assert_eq!(hybs.len(), 2);
        assert_eq!(hybs[0], Hybridization::SP);
        assert_eq!(hybs[1], Hybridization::SP);
    }

    #[test]
    fn test_lone_pair_upgrade_ester() {
        // Ester C(=O)-O: the single-bonded O should be upgraded to SP2
        let mut mol = Molecule::new("ester");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0)); // carbonyl C
        mol.atoms.push(Atom::new(1, "O", 1.2, 0.0, 0.0)); // carbonyl O
        mol.atoms.push(Atom::new(2, "O", -0.6, 1.0, 0.0)); // ester O
        mol.atoms.push(Atom::new(3, "C", -1.5, 1.5, 0.0)); // methyl C
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(2, 3, BondOrder::Single));

        let hybs = all_hybridizations(&mol);
        assert_eq!(hybs[0], Hybridization::SP2); // C with double bond
        assert_eq!(hybs[1], Hybridization::SP2); // O with double bond
        assert_eq!(hybs[2], Hybridization::SP2); // O upgraded (lone pair conjugation)
        assert_eq!(hybs[3], Hybridization::SP3); // Methyl C stays SP3
    }

    #[test]
    fn test_water_stays_sp3() {
        // Water: O with two H neighbors — should remain SP3
        let mut mol = Molecule::new("water");
        mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 0.96, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "H", -0.24, 0.93, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));

        let hybs = all_hybridizations(&mol);
        assert_eq!(hybs[0], Hybridization::SP3); // No SP2 neighbor to trigger upgrade
    }
}
