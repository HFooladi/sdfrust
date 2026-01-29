use std::collections::HashMap;

use crate::atom::Atom;
use crate::bond::{Bond, BondOrder};
use crate::collection::Collection;
use crate::sgroup::SGroup;
use crate::stereogroup::StereoGroup;

/// SDF format version.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum SdfFormat {
    /// V2000 format (fixed-width columns, max 999 atoms/bonds).
    #[default]
    V2000,
    /// V3000 format (variable-width, unlimited atoms/bonds).
    V3000,
}

impl SdfFormat {
    /// Returns the format string for the counts line.
    pub fn to_str(&self) -> &'static str {
        match self {
            SdfFormat::V2000 => "V2000",
            SdfFormat::V3000 => "V3000",
        }
    }

    /// Creates an SdfFormat from a string.
    pub fn parse(s: &str) -> Option<Self> {
        if s.contains("V3000") {
            Some(SdfFormat::V3000)
        } else if s.contains("V2000") {
            Some(SdfFormat::V2000)
        } else {
            None
        }
    }
}

/// Represents a molecule with atoms, bonds, and properties.
#[derive(Debug, Clone, PartialEq)]
pub struct Molecule {
    /// Molecule name (from the first line of the molfile).
    pub name: String,

    /// Program/timestamp line (second line of molfile).
    pub program_line: Option<String>,

    /// Comment line (third line of molfile).
    pub comment: Option<String>,

    /// List of atoms in the molecule.
    pub atoms: Vec<Atom>,

    /// List of bonds in the molecule.
    pub bonds: Vec<Bond>,

    /// Properties from the SDF data block (key-value pairs).
    pub properties: HashMap<String, String>,

    /// Format version (V2000 or V3000).
    pub format_version: SdfFormat,

    /// Enhanced stereochemistry groups (V3000).
    pub stereogroups: Vec<StereoGroup>,

    /// SGroups for superatoms, polymers, etc. (V3000).
    pub sgroups: Vec<SGroup>,

    /// Collections for atom lists, R-groups, etc. (V3000).
    pub collections: Vec<Collection>,
}

impl Molecule {
    /// Creates a new empty molecule with the given name.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            program_line: None,
            comment: None,
            atoms: Vec::new(),
            bonds: Vec::new(),
            properties: HashMap::new(),
            format_version: SdfFormat::V2000,
            stereogroups: Vec::new(),
            sgroups: Vec::new(),
            collections: Vec::new(),
        }
    }

    /// Returns true if this molecule requires V3000 format.
    ///
    /// Returns true if:
    /// - The molecule has more than 999 atoms or bonds
    /// - The molecule has V3000-only features (stereogroups, sgroups, collections)
    /// - Any atom has V3000-specific fields set
    /// - Any bond has extended bond types (coordination, hydrogen)
    pub fn needs_v3000(&self) -> bool {
        // V2000 limit is 999 atoms/bonds
        if self.atoms.len() > 999 || self.bonds.len() > 999 {
            return true;
        }

        // V3000-specific structures
        if !self.stereogroups.is_empty() || !self.sgroups.is_empty() || !self.collections.is_empty()
        {
            return true;
        }

        // V3000-only atom fields
        for atom in &self.atoms {
            if atom.v3000_id.is_some()
                || atom.atom_atom_mapping.is_some()
                || atom.rgroup_label.is_some()
                || atom.radical.is_some()
            {
                return true;
            }
        }

        // V3000-only bond types
        for bond in &self.bonds {
            if matches!(
                bond.order,
                crate::bond::BondOrder::Coordination | crate::bond::BondOrder::Hydrogen
            ) {
                return true;
            }
            if bond.v3000_id.is_some() || bond.reacting_center.is_some() {
                return true;
            }
        }

        false
    }

    /// Returns the number of atoms in the molecule.
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Returns the number of bonds in the molecule.
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Returns true if the molecule has no atoms.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Returns an iterator over all atoms.
    pub fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.atoms.iter()
    }

    /// Returns an iterator over all bonds.
    pub fn bonds(&self) -> impl Iterator<Item = &Bond> {
        self.bonds.iter()
    }

    /// Returns the atom at the given index, if it exists.
    pub fn get_atom(&self, index: usize) -> Option<&Atom> {
        self.atoms.get(index)
    }

    /// Returns all bonds connected to the given atom index.
    pub fn bonds_for_atom(&self, atom_index: usize) -> Vec<&Bond> {
        self.bonds
            .iter()
            .filter(|b| b.contains_atom(atom_index))
            .collect()
    }

    /// Returns the indices of atoms connected to the given atom.
    pub fn neighbors(&self, atom_index: usize) -> Vec<usize> {
        self.bonds
            .iter()
            .filter_map(|b| b.other_atom(atom_index))
            .collect()
    }

    /// Returns the molecular formula as a string (e.g., "C6H12O6").
    pub fn formula(&self) -> String {
        let mut counts: HashMap<&str, usize> = HashMap::new();
        for atom in &self.atoms {
            *counts.entry(atom.element.as_str()).or_insert(0) += 1;
        }

        // Standard order: C, H, then alphabetical
        let mut formula = String::new();

        if let Some(&c) = counts.get("C") {
            formula.push('C');
            if c > 1 {
                formula.push_str(&c.to_string());
            }
            counts.remove("C");
        }

        if let Some(&h) = counts.get("H") {
            formula.push('H');
            if h > 1 {
                formula.push_str(&h.to_string());
            }
            counts.remove("H");
        }

        let mut remaining: Vec<_> = counts.into_iter().collect();
        remaining.sort_by_key(|(elem, _)| *elem);

        for (elem, count) in remaining {
            formula.push_str(elem);
            if count > 1 {
                formula.push_str(&count.to_string());
            }
        }

        formula
    }

    /// Returns the total formal charge of the molecule.
    pub fn total_charge(&self) -> i32 {
        self.atoms.iter().map(|a| a.formal_charge as i32).sum()
    }

    /// Returns the centroid (geometric center) of the molecule.
    pub fn centroid(&self) -> Option<(f64, f64, f64)> {
        if self.atoms.is_empty() {
            return None;
        }

        let n = self.atoms.len() as f64;
        let sum_x: f64 = self.atoms.iter().map(|a| a.x).sum();
        let sum_y: f64 = self.atoms.iter().map(|a| a.y).sum();
        let sum_z: f64 = self.atoms.iter().map(|a| a.z).sum();

        Some((sum_x / n, sum_y / n, sum_z / n))
    }

    /// Translates the molecule by the given vector.
    pub fn translate(&mut self, dx: f64, dy: f64, dz: f64) {
        for atom in &mut self.atoms {
            atom.x += dx;
            atom.y += dy;
            atom.z += dz;
        }
    }

    /// Centers the molecule at the origin.
    pub fn center(&mut self) {
        if let Some((cx, cy, cz)) = self.centroid() {
            self.translate(-cx, -cy, -cz);
        }
    }

    /// Returns a property value by key, if it exists.
    pub fn get_property(&self, key: &str) -> Option<&str> {
        self.properties.get(key).map(|s| s.as_str())
    }

    /// Sets a property value.
    pub fn set_property(&mut self, key: &str, value: &str) {
        self.properties.insert(key.to_string(), value.to_string());
    }

    /// Returns true if the molecule contains any aromatic bonds.
    pub fn has_aromatic_bonds(&self) -> bool {
        self.bonds.iter().any(|b| b.is_aromatic())
    }

    /// Returns true if the molecule contains any charged atoms.
    pub fn has_charges(&self) -> bool {
        self.atoms.iter().any(|a| a.is_charged())
    }

    /// Returns the count of each element in the molecule.
    pub fn element_counts(&self) -> HashMap<String, usize> {
        let mut counts = HashMap::new();
        for atom in &self.atoms {
            *counts.entry(atom.element.clone()).or_insert(0) += 1;
        }
        counts
    }

    /// Calculates the sum of bond orders (useful for validation).
    pub fn total_bond_order(&self) -> f64 {
        self.bonds.iter().map(|b| b.order.order()).sum()
    }

    /// Returns atoms that match the given element.
    pub fn atoms_by_element(&self, element: &str) -> Vec<&Atom> {
        self.atoms.iter().filter(|a| a.element == element).collect()
    }

    /// Returns bonds with the given order.
    pub fn bonds_by_order(&self, order: BondOrder) -> Vec<&Bond> {
        self.bonds.iter().filter(|b| b.order == order).collect()
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new("")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};

    #[test]
    fn test_needs_v3000_basic() {
        let mol = Molecule::new("test");
        assert!(!mol.needs_v3000());
    }

    #[test]
    fn test_needs_v3000_stereogroups() {
        let mut mol = Molecule::new("test");
        mol.stereogroups.push(StereoGroup::default());
        assert!(mol.needs_v3000());
    }

    #[test]
    fn test_needs_v3000_coordination_bond() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "N", 1.0, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Coordination));
        assert!(mol.needs_v3000());
    }
}
