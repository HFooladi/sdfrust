use std::collections::HashMap;

use crate::atom::Atom;
use crate::bond::{Bond, BondOrder};

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
        }
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
        self.atoms
            .iter()
            .filter(|a| a.element == element)
            .collect()
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
