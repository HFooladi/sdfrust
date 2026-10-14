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

use crate::descriptors::conjugation::{AtomEnvironment, conjugated_bonds};
use crate::descriptors::elements::get_element;
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
    /// Returns the index into OGB's hybridization list
    /// `[SP, SP2, SP3, SP3D, SP3D2, misc]`; `S` and `Other` map to misc (5).
    pub fn to_ogb_index(&self) -> u8 {
        match self {
            Hybridization::SP => 0,
            Hybridization::SP2 => 1,
            Hybridization::SP3 => 2,
            Hybridization::SP3D => 3,
            Hybridization::SP3D2 => 4,
            Hybridization::S | Hybridization::Other => 5,
        }
    }
}

/// Number of σ bonds plus lone pairs (RDKit's `numBondsPlusLonePairs`).
fn bonds_plus_lone_pairs(mol: &Molecule, env: &AtomEnvironment, idx: usize) -> i32 {
    let atom = &mol.atoms[idx];
    let degree = env.total_degree(idx);
    let Some(data) = get_element(&atom.element) else {
        return degree;
    };
    if data.atomic_number <= 1 {
        return degree;
    }
    if data.atomic_number >= 89 {
        return degree; // actinides and beyond: just use the degree
    }
    let outer = data.outer_electrons as i32;
    let valence = env.total_valence(idx);
    let charge = atom.formal_charge as i32;
    let free_electrons = outer - (valence + charge);
    if valence + outer - charge < 8 {
        // Below an octet: unpaired electrons occupy orbitals too
        let radicals = atom.radical_electrons() as i32;
        degree + (free_electrons - radicals) / 2 + radicals
    } else {
        degree + free_electrons / 2
    }
}

/// Hybridization of every atom from precomputed environments (RDKit's
/// `setHybridization`).
fn hybridizations(
    mol: &Molecule,
    env: &AtomEnvironment,
    conjugated: &[bool],
) -> Vec<Hybridization> {
    (0..mol.atom_count())
        .map(|i| match bonds_plus_lone_pairs(mol, env, i) {
            i32::MIN..=1 => Hybridization::S,
            2 => Hybridization::SP,
            3 => Hybridization::SP2,
            4 => {
                // An sp3-like atom next to a π system (amide N, ester O,
                // pyrrole N) is sp2, unless it has 4 explicit neighbours.
                let has_conjugated_bond = env.adj.neighbors(i).iter().any(|&(_, b)| conjugated[b]);
                if env.adj.degree(i) < 4 && has_conjugated_bond {
                    Hybridization::SP2
                } else {
                    Hybridization::SP3
                }
            }
            5 => Hybridization::SP3D,
            6 => Hybridization::SP3D2,
            _ => Hybridization::Other,
        })
        .collect()
}

/// Determine the hybridization of atom `idx`.
///
/// Uses the number of σ bonds plus lone pairs, as RDKit does: 2 → SP,
/// 3 → SP2, 4 → SP3 (SP2 if the atom has a conjugated bond and fewer than four
/// explicit neighbours), 5 → SP3D, 6 → SP3D2; hydrogen and isolated ions are S.
/// For many atoms, prefer [`all_hybridizations`], which perceives the molecule once.
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
    all_hybridizations(mol)[idx]
}

/// Compute hybridization for all atoms in the molecule.
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
    let env = AtomEnvironment::new(mol);
    let conjugated = conjugated_bonds(mol, &env);
    hybridizations(mol, &env, &conjugated)
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
        assert_eq!(Hybridization::SP.to_ogb_index(), 0);
        assert_eq!(Hybridization::SP2.to_ogb_index(), 1);
        assert_eq!(Hybridization::SP3.to_ogb_index(), 2);
        assert_eq!(Hybridization::SP3D.to_ogb_index(), 3);
        assert_eq!(Hybridization::SP3D2.to_ogb_index(), 4);
        assert_eq!(Hybridization::S.to_ogb_index(), 5);
        assert_eq!(Hybridization::Other.to_ogb_index(), 5);
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
