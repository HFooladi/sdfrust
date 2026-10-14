//! Double-bond (E/Z) stereo perception from coordinates.
//!
//! For each non-aromatic double bond whose two ends both carry distinguishable
//! substituents, the higher-ranked substituent on each side (by CIP-style
//! ranking) is located from the 2D/3D coordinates: on the same side gives Z,
//! on opposite sides gives E. Bonds marked "either" in the file are
//! [`DoubleBondStereo::Any`]. Like RDKit, double bonds in rings smaller than
//! eight atoms get no stereo.

use crate::bond::{BondOrder, BondStereo};
use crate::descriptors::aromaticity::all_aromatic_bonds;
use crate::descriptors::chirality::{CipPriority, compute_cip_priority};
use crate::descriptors::rings::sssr;
use crate::descriptors::valence::all_implicit_hydrogen_counts;
use crate::graph::AdjacencyList;
use crate::molecule::Molecule;

/// Configuration of a double bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DoubleBondStereo {
    /// Not a stereo double bond (or not a double bond).
    None,
    /// Higher-priority substituents on the same side (cis-like).
    Z,
    /// Higher-priority substituents on opposite sides (trans-like).
    E,
    /// Stereo bond whose configuration is explicitly unknown ("either").
    Any,
}

/// Perceive E/Z configuration for every bond (`None` for non-stereo bonds).
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::stereo::{DoubleBondStereo, double_bond_stereo};
///
/// // trans-2-butene drawn in 2D
/// let mut mol = Molecule::new("e-2-butene");
/// mol.atoms.push(Atom::new(0, "C", -1.3, 0.75, 0.0));
/// mol.atoms.push(Atom::new(1, "C", -0.65, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "C", 0.65, 0.0, 0.0));
/// mol.atoms.push(Atom::new(3, "C", 1.3, -0.75, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
/// mol.bonds.push(Bond::new(1, 2, BondOrder::Double));
/// mol.bonds.push(Bond::new(2, 3, BondOrder::Single));
///
/// assert_eq!(double_bond_stereo(&mol)[1], DoubleBondStereo::E);
/// ```
pub fn double_bond_stereo(mol: &Molecule) -> Vec<DoubleBondStereo> {
    let m = mol.bond_count();
    let mut result = vec![DoubleBondStereo::None; m];
    if !mol.bonds.iter().any(|b| b.order == BondOrder::Double) {
        return result;
    }

    let adj = AdjacencyList::from_molecule(mol);
    let aromatic = all_aromatic_bonds(mol);
    let implicit_h = all_implicit_hydrogen_counts(mol);
    let mut small_ring_bond = vec![false; m];
    for ring in sssr(mol).iter().filter(|r| r.size() < 8) {
        for &b in &ring.bonds {
            small_ring_bond[b] = true;
        }
    }

    for (b, bond) in mol.bonds.iter().enumerate() {
        if bond.order != BondOrder::Double || aromatic[b] || small_ring_bond[b] {
            continue;
        }
        let (a1, a2) = (bond.atom1, bond.atom2);
        // Cumulated double bonds (allenes, azides) are not stereo bonds
        let cumulated = |atom: usize| {
            adj.neighbors(atom).iter().any(|&(_, other)| {
                other != b
                    && matches!(
                        mol.bonds[other].order,
                        BondOrder::Double | BondOrder::Triple
                    )
            })
        };
        if cumulated(a1) || cumulated(a2) {
            continue;
        }
        // An explicit "either" flag in the file is kept as-is (as RDKit does)
        if bond.stereo == BondStereo::Either {
            result[b] = DoubleBondStereo::Any;
            continue;
        }
        let (Some(s1), Some(s2)) = (
            ranked_substituent(mol, &adj, &implicit_h, a1, a2),
            ranked_substituent(mol, &adj, &implicit_h, a2, a1),
        ) else {
            continue;
        };
        result[b] = match side(mol, s1, a1, a2, s2) {
            Some(true) => DoubleBondStereo::Z,
            Some(false) => DoubleBondStereo::E,
            None => DoubleBondStereo::None,
        };
    }
    result
}

/// The single highest-ranked substituent of `atom` (other than `partner`),
/// or `None` if the end cannot be stereogenic (no substituent, more than two,
/// or two equal ones).
fn ranked_substituent(
    mol: &Molecule,
    adj: &AdjacencyList,
    implicit_h: &[u8],
    atom: usize,
    partner: usize,
) -> Option<usize> {
    let others: Vec<usize> = adj
        .neighbor_atoms(atom)
        .into_iter()
        .filter(|&j| j != partner)
        .collect();
    let total = others.len() + implicit_h[atom] as usize;
    match (others.len(), total) {
        (1, 1) | (1, 2) => Some(others[0]),
        (2, 2) => {
            let p: Vec<CipPriority> = others
                .iter()
                .map(|&j| compute_cip_priority(mol, adj, atom, j))
                .collect();
            match p[0].cmp(&p[1]) {
                std::cmp::Ordering::Greater => Some(others[0]),
                std::cmp::Ordering::Less => Some(others[1]),
                std::cmp::Ordering::Equal => None,
            }
        }
        _ => None,
    }
}

/// Whether `x` (on `a`) and `y` (on `b`) are on the same side of the a=b
/// bond (`|dihedral x-a-b-y| < 90°`). `None` if the geometry is degenerate.
fn side(mol: &Molecule, x: usize, a: usize, b: usize, y: usize) -> Option<bool> {
    let p = |i: usize| {
        let at = &mol.atoms[i];
        [at.x, at.y, at.z]
    };
    let sub = |u: [f64; 3], v: [f64; 3]| [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
    let cross = |u: [f64; 3], v: [f64; 3]| {
        [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ]
    };
    let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    let axis = sub(p(b), p(a));
    let n1 = cross(axis, sub(p(x), p(a)));
    let n2 = cross(axis, sub(p(y), p(b)));
    let d = dot(n1, n2);
    let scale = dot(n1, n1).sqrt() * dot(n2, n2).sqrt();
    if scale < 1e-8 || (d / scale).abs() < 1e-3 {
        return None;
    }
    Some(d > 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::Bond;

    fn butene(y_end: f64) -> Molecule {
        let mut mol = Molecule::new("2-butene");
        mol.atoms.push(Atom::new(0, "C", -1.3, 0.75, 0.0));
        mol.atoms.push(Atom::new(1, "C", -0.65, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 0.65, 0.0, 0.0));
        mol.atoms.push(Atom::new(3, "C", 1.3, y_end, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Double));
        mol.bonds.push(Bond::new(2, 3, BondOrder::Single));
        mol
    }

    #[test]
    fn test_e_and_z() {
        assert_eq!(double_bond_stereo(&butene(-0.75))[1], DoubleBondStereo::E);
        assert_eq!(double_bond_stereo(&butene(0.75))[1], DoubleBondStereo::Z);
    }

    #[test]
    fn test_either_is_any() {
        let mut mol = butene(0.75);
        mol.bonds[1].stereo = BondStereo::Either;
        assert_eq!(double_bond_stereo(&mol)[1], DoubleBondStereo::Any);
    }

    #[test]
    fn test_terminal_alkene_has_no_stereo() {
        // propene: CH2= end has two hydrogens
        let mut mol = butene(0.75);
        mol.atoms.pop();
        mol.bonds.pop();
        assert_eq!(double_bond_stereo(&mol)[1], DoubleBondStereo::None);
    }

    #[test]
    fn test_identical_substituents_have_no_stereo() {
        // 2-methylbut-2-ene: one end carries two methyls
        let mut mol = butene(0.75);
        mol.atoms.push(Atom::new(4, "C", 1.3, -0.75, 0.0));
        mol.bonds.push(Bond::new(2, 4, BondOrder::Single));
        assert_eq!(double_bond_stereo(&mol)[1], DoubleBondStereo::None);
    }
}
