//! Python wrapper for Molecule struct.

use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use std::collections::HashMap;

use sdfrust::{Atom, Bond, BondOrder, Molecule, SdfFormat};

use crate::atom::PyAtom;
use crate::bond::{PyBond, PyBondOrder};

/// SDF format version.
#[pyclass(name = "SdfFormat")]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct PySdfFormat {
    pub(crate) inner: SdfFormat,
}

impl From<SdfFormat> for PySdfFormat {
    fn from(format: SdfFormat) -> Self {
        PySdfFormat { inner: format }
    }
}

impl From<PySdfFormat> for SdfFormat {
    fn from(format: PySdfFormat) -> Self {
        format.inner
    }
}

#[pymethods]
impl PySdfFormat {
    /// V2000 format (fixed-width columns, max 999 atoms/bonds).
    #[staticmethod]
    pub fn v2000() -> Self {
        PySdfFormat {
            inner: SdfFormat::V2000,
        }
    }

    /// V3000 format (variable-width, unlimited atoms/bonds).
    #[staticmethod]
    pub fn v3000() -> Self {
        PySdfFormat {
            inner: SdfFormat::V3000,
        }
    }

    /// Returns the format string ("V2000" or "V3000").
    pub fn to_str(&self) -> &'static str {
        self.inner.to_str()
    }

    fn __repr__(&self) -> String {
        match self.inner {
            SdfFormat::V2000 => "SdfFormat.V2000".to_string(),
            SdfFormat::V3000 => "SdfFormat.V3000".to_string(),
        }
    }

    fn __str__(&self) -> String {
        self.inner.to_str().to_string()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }
}

/// A molecule with atoms, bonds, and properties.
///
/// The Molecule class is the main container for chemical structure data.
/// It holds atoms, bonds, and key-value properties parsed from SDF or MOL2 files.
#[pyclass(name = "Molecule")]
#[derive(Clone)]
pub struct PyMolecule {
    pub(crate) inner: Molecule,
}

impl From<Molecule> for PyMolecule {
    fn from(mol: Molecule) -> Self {
        PyMolecule { inner: mol }
    }
}

impl From<PyMolecule> for Molecule {
    fn from(mol: PyMolecule) -> Self {
        mol.inner
    }
}

#[pymethods]
impl PyMolecule {
    /// Create a new empty molecule with the given name.
    ///
    /// Args:
    ///     name: The molecule name.
    ///
    /// Returns:
    ///     A new Molecule instance.
    #[new]
    #[pyo3(signature = (name=""))]
    pub fn new(name: &str) -> Self {
        PyMolecule {
            inner: Molecule::new(name),
        }
    }

    // ============================================================
    // Basic Properties
    // ============================================================

    /// The molecule name (from the first line of the molfile).
    #[getter]
    pub fn name(&self) -> &str {
        &self.inner.name
    }

    /// Set the molecule name.
    #[setter]
    pub fn set_name(&mut self, name: &str) {
        self.inner.name = name.to_string();
    }

    /// The number of atoms in the molecule.
    #[getter]
    pub fn num_atoms(&self) -> usize {
        self.inner.atom_count()
    }

    /// The number of bonds in the molecule.
    #[getter]
    pub fn num_bonds(&self) -> usize {
        self.inner.bond_count()
    }

    /// List of all atoms in the molecule.
    #[getter]
    pub fn atoms(&self) -> Vec<PyAtom> {
        self.inner.atoms.iter().map(PyAtom::from).collect()
    }

    /// List of all bonds in the molecule.
    #[getter]
    pub fn bonds(&self) -> Vec<PyBond> {
        self.inner.bonds.iter().map(PyBond::from).collect()
    }

    /// Dictionary of properties from the SDF data block.
    #[getter]
    pub fn properties(&self) -> HashMap<String, String> {
        self.inner.properties.clone()
    }

    /// The SDF format version (V2000 or V3000).
    #[getter]
    pub fn format_version(&self) -> PySdfFormat {
        PySdfFormat::from(self.inner.format_version)
    }

    /// Set the SDF format version.
    #[setter]
    pub fn set_format_version(&mut self, format: PySdfFormat) {
        self.inner.format_version = format.inner;
    }

    /// Program/timestamp line (second line of molfile).
    #[getter]
    pub fn program_line(&self) -> Option<String> {
        self.inner.program_line.clone()
    }

    /// Comment line (third line of molfile).
    #[getter]
    pub fn comment(&self) -> Option<String> {
        self.inner.comment.clone()
    }

    // ============================================================
    // Atom Access
    // ============================================================

    /// Get the atom at the given index.
    ///
    /// Args:
    ///     index: The atom index (0-based).
    ///
    /// Returns:
    ///     The Atom at the given index.
    ///
    /// Raises:
    ///     IndexError: If the index is out of bounds.
    pub fn get_atom(&self, index: usize) -> PyResult<PyAtom> {
        self.inner
            .get_atom(index)
            .map(PyAtom::from)
            .ok_or_else(|| PyIndexError::new_err(format!("Atom index {} out of range", index)))
    }

    /// Add an atom to the molecule.
    ///
    /// Args:
    ///     atom: The Atom to add.
    pub fn add_atom(&mut self, atom: PyAtom) {
        self.inner.atoms.push(atom.inner);
    }

    /// Add a bond to the molecule.
    ///
    /// Args:
    ///     bond: The Bond to add.
    pub fn add_bond(&mut self, bond: PyBond) {
        self.inner.bonds.push(bond.inner);
    }

    /// Get the indices of atoms connected to the given atom.
    ///
    /// Args:
    ///     atom_index: The atom index to find neighbors for.
    ///
    /// Returns:
    ///     List of neighbor atom indices.
    pub fn neighbors(&self, atom_index: usize) -> Vec<usize> {
        self.inner.neighbors(atom_index)
    }

    /// Get all bonds connected to the given atom.
    ///
    /// Args:
    ///     atom_index: The atom index.
    ///
    /// Returns:
    ///     List of bonds involving this atom.
    pub fn bonds_for_atom(&self, atom_index: usize) -> Vec<PyBond> {
        self.inner
            .bonds_for_atom(atom_index)
            .into_iter()
            .map(PyBond::from)
            .collect()
    }

    // ============================================================
    // Property Access
    // ============================================================

    /// Get a property value by key.
    ///
    /// Args:
    ///     key: The property key.
    ///
    /// Returns:
    ///     The property value, or None if not found.
    pub fn get_property(&self, key: &str) -> Option<String> {
        self.inner.get_property(key).map(|s| s.to_string())
    }

    /// Set a property value.
    ///
    /// Args:
    ///     key: The property key.
    ///     value: The property value.
    pub fn set_property(&mut self, key: &str, value: &str) {
        self.inner.set_property(key, value);
    }

    // ============================================================
    // Molecular Formula & Charges
    // ============================================================

    /// Returns the molecular formula as a string (e.g., "C6H12O6").
    ///
    /// Elements are ordered: C, H, then alphabetically.
    pub fn formula(&self) -> String {
        self.inner.formula()
    }

    /// Returns the total formal charge of the molecule.
    pub fn total_charge(&self) -> i32 {
        self.inner.total_charge()
    }

    /// Returns True if the molecule contains any charged atoms.
    pub fn has_charges(&self) -> bool {
        self.inner.has_charges()
    }

    /// Returns True if the molecule contains any aromatic bonds.
    pub fn has_aromatic_bonds(&self) -> bool {
        self.inner.has_aromatic_bonds()
    }

    // ============================================================
    // Geometry
    // ============================================================

    /// Returns the centroid (geometric center) of the molecule.
    ///
    /// Returns:
    ///     A tuple (x, y, z) of the centroid coordinates, or None if the molecule is empty.
    pub fn centroid(&self) -> Option<(f64, f64, f64)> {
        self.inner.centroid()
    }

    /// Translate the molecule by the given vector.
    ///
    /// Args:
    ///     dx: Translation in x direction.
    ///     dy: Translation in y direction.
    ///     dz: Translation in z direction.
    pub fn translate(&mut self, dx: f64, dy: f64, dz: f64) {
        self.inner.translate(dx, dy, dz);
    }

    /// Center the molecule at the origin (move centroid to (0, 0, 0)).
    pub fn center(&mut self) {
        self.inner.center();
    }

    // ============================================================
    // Element & Bond Filtering
    // ============================================================

    /// Returns a count of each element in the molecule.
    ///
    /// Returns:
    ///     A dictionary mapping element symbols to counts.
    pub fn element_counts(&self) -> HashMap<String, usize> {
        self.inner.element_counts()
    }

    /// Returns atoms that match the given element.
    ///
    /// Args:
    ///     element: The element symbol to filter by (e.g., "C", "N").
    ///
    /// Returns:
    ///     List of atoms with the given element.
    pub fn atoms_by_element(&self, element: &str) -> Vec<PyAtom> {
        self.inner
            .atoms_by_element(element)
            .into_iter()
            .map(PyAtom::from)
            .collect()
    }

    /// Returns bonds with the given order.
    ///
    /// Args:
    ///     order: The bond order to filter by.
    ///
    /// Returns:
    ///     List of bonds with the given order.
    pub fn bonds_by_order(&self, order: PyBondOrder) -> Vec<PyBond> {
        self.inner
            .bonds_by_order(order.inner)
            .into_iter()
            .map(PyBond::from)
            .collect()
    }

    // ============================================================
    // Descriptors
    // ============================================================

    /// Calculate the molecular weight (sum of atomic weights).
    ///
    /// Uses standard atomic weights (IUPAC 2021) for each element.
    ///
    /// Returns:
    ///     The molecular weight in g/mol, or None if any atom has an unknown element.
    pub fn molecular_weight(&self) -> Option<f64> {
        self.inner.molecular_weight()
    }

    /// Calculate the exact (monoisotopic) mass.
    ///
    /// Uses the mass of the most abundant isotope for each element.
    ///
    /// Returns:
    ///     The exact mass in g/mol, or None if any atom has an unknown element.
    pub fn exact_mass(&self) -> Option<f64> {
        self.inner.exact_mass()
    }

    /// Count non-hydrogen atoms (heavy atoms).
    ///
    /// Heavy atoms are all atoms except hydrogen (H), deuterium (D), and tritium (T).
    pub fn heavy_atom_count(&self) -> usize {
        self.inner.heavy_atom_count()
    }

    /// Count the number of rings in the molecule.
    ///
    /// Uses the Euler characteristic formula: rings = bonds - atoms + components.
    pub fn ring_count(&self) -> usize {
        self.inner.ring_count()
    }

    /// Count rotatable bonds.
    ///
    /// A bond is rotatable if it is a single bond, not in a ring,
    /// not terminal, and doesn't involve hydrogen atoms.
    pub fn rotatable_bond_count(&self) -> usize {
        self.inner.rotatable_bond_count()
    }

    /// Check if an atom at the given index is in a ring.
    ///
    /// Args:
    ///     index: The atom index.
    ///
    /// Returns:
    ///     True if the atom is in a ring, False otherwise.
    pub fn is_atom_in_ring(&self, index: usize) -> bool {
        self.inner.is_atom_in_ring(index)
    }

    /// Check if a bond at the given index is in a ring.
    ///
    /// Args:
    ///     index: The bond index.
    ///
    /// Returns:
    ///     True if the bond is in a ring, False otherwise.
    pub fn is_bond_in_ring(&self, index: usize) -> bool {
        self.inner.is_bond_in_ring(index)
    }

    /// Count bonds by bond order.
    ///
    /// Returns:
    ///     A dictionary mapping BondOrder to count.
    pub fn bond_type_counts(&self) -> HashMap<String, usize> {
        let counts = self.inner.bond_type_counts();
        counts
            .into_iter()
            .map(|(order, count)| {
                let name = match order {
                    BondOrder::Single => "single",
                    BondOrder::Double => "double",
                    BondOrder::Triple => "triple",
                    BondOrder::Aromatic => "aromatic",
                    BondOrder::SingleOrDouble => "single_or_double",
                    BondOrder::SingleOrAromatic => "single_or_aromatic",
                    BondOrder::DoubleOrAromatic => "double_or_aromatic",
                    BondOrder::Any => "any",
                    BondOrder::Coordination => "coordination",
                    BondOrder::Hydrogen => "hydrogen",
                };
                (name.to_string(), count)
            })
            .collect()
    }

    // ============================================================
    // Format Detection
    // ============================================================

    /// Returns True if this molecule requires V3000 format.
    ///
    /// Returns True if:
    /// - The molecule has more than 999 atoms or bonds
    /// - The molecule has V3000-only features (stereogroups, sgroups, collections)
    /// - Any atom has V3000-specific fields set
    /// - Any bond has extended bond types (coordination, hydrogen)
    pub fn needs_v3000(&self) -> bool {
        self.inner.needs_v3000()
    }

    // ============================================================
    // Utility
    // ============================================================

    /// Returns True if the molecule has no atoms.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Returns the sum of bond orders (useful for validation).
    pub fn total_bond_order(&self) -> f64 {
        self.inner.total_bond_order()
    }

    // ============================================================
    // NumPy support (when feature is enabled)
    // ============================================================

    #[cfg(feature = "numpy")]
    /// Get atom coordinates as a NumPy array.
    ///
    /// Returns:
    ///     A NumPy array of shape (N, 3) where N is the number of atoms.
    ///     Each row contains [x, y, z] coordinates in Angstroms.
    pub fn get_coords_array<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, numpy::PyArray2<f64>>> {
        use numpy::PyArray2;

        // Convert to Vec<Vec<f64>> for from_vec2
        let coords: Vec<Vec<f64>> = self
            .inner
            .atoms
            .iter()
            .map(|a| vec![a.x, a.y, a.z])
            .collect();

        PyArray2::from_vec2(py, &coords)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[cfg(feature = "numpy")]
    /// Set atom coordinates from a NumPy array.
    ///
    /// Args:
    ///     coords: A NumPy array of shape (N, 3) where N is the number of atoms.
    ///
    /// Raises:
    ///     ValueError: If the array shape doesn't match (num_atoms, 3).
    pub fn set_coords_array(&mut self, coords: &Bound<'_, numpy::PyArray2<f64>>) -> PyResult<()> {
        use numpy::{PyArrayMethods, PyUntypedArrayMethods};

        let shape = coords.shape();
        let n_atoms = self.inner.atoms.len();

        if shape.len() != 2 || shape[0] != n_atoms || shape[1] != 3 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Expected array of shape ({}, 3), got {:?}",
                n_atoms, shape
            )));
        }

        let readonly = coords.readonly();
        let data = readonly.as_slice()?;

        for (i, atom) in self.inner.atoms.iter_mut().enumerate() {
            atom.x = data[i * 3];
            atom.y = data[i * 3 + 1];
            atom.z = data[i * 3 + 2];
        }

        Ok(())
    }

    #[cfg(feature = "numpy")]
    /// Get atomic numbers as a NumPy array.
    ///
    /// Returns:
    ///     A NumPy array of shape (N,) containing atomic numbers.
    ///     Unknown elements are assigned atomic number 0.
    pub fn get_atomic_numbers<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<u8>> {
        let atomic_nums: Vec<u8> = self
            .inner
            .atoms
            .iter()
            .map(|a| element_to_atomic_number(&a.element))
            .collect();

        numpy::PyArray1::from_vec(py, atomic_nums)
    }

    // ============================================================
    // Geometry operations (when feature is enabled)
    // ============================================================

    #[cfg(feature = "geometry")]
    /// Rotate the molecule around an axis by a given angle.
    ///
    /// The rotation is performed around the origin. Center the molecule first
    /// if rotation around the centroid is desired.
    ///
    /// Args:
    ///     axis: The rotation axis as [x, y, z] (will be normalized).
    ///     angle: The rotation angle in radians.
    ///
    /// Example:
    ///     >>> import math
    ///     >>> mol.rotate([0, 0, 1], math.pi / 2)  # 90° around Z
    pub fn rotate(&mut self, axis: [f64; 3], angle: f64) {
        self.inner.rotate(axis, angle);
    }

    #[cfg(feature = "geometry")]
    /// Apply a 3x3 rotation matrix to the molecule.
    ///
    /// Args:
    ///     matrix: A 3x3 rotation matrix as a list of lists.
    ///
    /// Example:
    ///     >>> identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    ///     >>> mol.apply_rotation_matrix(identity)
    pub fn apply_rotation_matrix(&mut self, matrix: [[f64; 3]; 3]) {
        self.inner.apply_rotation_matrix(&matrix);
    }

    #[cfg(feature = "geometry")]
    /// Apply a rotation matrix and translation to the molecule.
    ///
    /// First applies the rotation, then the translation.
    ///
    /// Args:
    ///     rotation: A 3x3 rotation matrix as a list of lists.
    ///     translation: A translation vector [dx, dy, dz].
    ///
    /// Example:
    ///     >>> identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    ///     >>> mol.apply_transform(identity, [1, 2, 3])
    pub fn apply_transform(&mut self, rotation: [[f64; 3]; 3], translation: [f64; 3]) {
        self.inner.apply_transform(&rotation, translation);
    }

    #[cfg(feature = "geometry")]
    /// Compute the pairwise distance matrix for all atoms.
    ///
    /// Returns an NxN matrix where entry [i][j] is the Euclidean distance
    /// between atom i and atom j in Angstroms.
    ///
    /// Returns:
    ///     A list of lists containing pairwise distances.
    ///
    /// Example:
    ///     >>> matrix = mol.distance_matrix()
    ///     >>> print(matrix[0][1])  # Distance between atoms 0 and 1
    pub fn distance_matrix(&self) -> Vec<Vec<f64>> {
        self.inner.distance_matrix()
    }

    #[cfg(feature = "geometry")]
    /// Calculate RMSD to another molecule.
    ///
    /// Computes the root mean square deviation of atomic positions.
    /// The molecules must have the same number of atoms.
    /// No alignment is performed - atoms are compared directly by index.
    ///
    /// Args:
    ///     other: The other molecule to compare to.
    ///
    /// Returns:
    ///     The RMSD value in Angstroms.
    ///
    /// Raises:
    ///     ValueError: If the molecules have different atom counts.
    ///
    /// Example:
    ///     >>> rmsd = mol1.rmsd_to(mol2)
    pub fn rmsd_to(&self, other: &PyMolecule) -> PyResult<f64> {
        self.inner
            .rmsd_to(&other.inner)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "Molecule(name='{}', atoms={}, bonds={})",
            self.inner.name,
            self.inner.atom_count(),
            self.inner.bond_count()
        )
    }

    fn __str__(&self) -> String {
        format!(
            "{} ({} atoms, {} bonds)",
            self.inner.name,
            self.inner.atom_count(),
            self.inner.bond_count()
        )
    }

    fn __len__(&self) -> usize {
        self.inner.atom_count()
    }

    fn __bool__(&self) -> bool {
        !self.inner.is_empty()
    }
}

// Helper function to create a molecule from atoms and bonds
impl PyMolecule {
    /// Create a PyMolecule from a name and atoms/bonds.
    #[allow(dead_code)]
    pub fn from_parts(
        name: &str,
        atoms: Vec<Atom>,
        bonds: Vec<Bond>,
        properties: HashMap<String, String>,
    ) -> Self {
        let mut mol = Molecule::new(name);
        mol.atoms = atoms;
        mol.bonds = bonds;
        mol.properties = properties;
        PyMolecule { inner: mol }
    }
}

/// Convert element symbol to atomic number.
/// Returns 0 for unknown elements.
#[cfg(feature = "numpy")]
fn element_to_atomic_number(element: &str) -> u8 {
    match element.trim() {
        "H" => 1,
        "He" => 2,
        "Li" => 3,
        "Be" => 4,
        "B" => 5,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "Ne" => 10,
        "Na" => 11,
        "Mg" => 12,
        "Al" => 13,
        "Si" => 14,
        "P" => 15,
        "S" => 16,
        "Cl" => 17,
        "Ar" => 18,
        "K" => 19,
        "Ca" => 20,
        "Sc" => 21,
        "Ti" => 22,
        "V" => 23,
        "Cr" => 24,
        "Mn" => 25,
        "Fe" => 26,
        "Co" => 27,
        "Ni" => 28,
        "Cu" => 29,
        "Zn" => 30,
        "Ga" => 31,
        "Ge" => 32,
        "As" => 33,
        "Se" => 34,
        "Br" => 35,
        "Kr" => 36,
        "Rb" => 37,
        "Sr" => 38,
        "Y" => 39,
        "Zr" => 40,
        "Nb" => 41,
        "Mo" => 42,
        "Tc" => 43,
        "Ru" => 44,
        "Rh" => 45,
        "Pd" => 46,
        "Ag" => 47,
        "Cd" => 48,
        "In" => 49,
        "Sn" => 50,
        "Sb" => 51,
        "Te" => 52,
        "I" => 53,
        "Xe" => 54,
        "Cs" => 55,
        "Ba" => 56,
        "La" => 57,
        "Ce" => 58,
        "Pr" => 59,
        "Nd" => 60,
        "Pm" => 61,
        "Sm" => 62,
        "Eu" => 63,
        "Gd" => 64,
        "Tb" => 65,
        "Dy" => 66,
        "Ho" => 67,
        "Er" => 68,
        "Tm" => 69,
        "Yb" => 70,
        "Lu" => 71,
        "Hf" => 72,
        "Ta" => 73,
        "W" => 74,
        "Re" => 75,
        "Os" => 76,
        "Ir" => 77,
        "Pt" => 78,
        "Au" => 79,
        "Hg" => 80,
        "Tl" => 81,
        "Pb" => 82,
        "Bi" => 83,
        "Po" => 84,
        "At" => 85,
        "Rn" => 86,
        "Fr" => 87,
        "Ra" => 88,
        "Ac" => 89,
        "Th" => 90,
        "Pa" => 91,
        "U" => 92,
        "Np" => 93,
        "Pu" => 94,
        "Am" => 95,
        "Cm" => 96,
        "Bk" => 97,
        "Cf" => 98,
        "Es" => 99,
        "Fm" => 100,
        "Md" => 101,
        "No" => 102,
        "Lr" => 103,
        "Rf" => 104,
        "Db" => 105,
        "Sg" => 106,
        "Bh" => 107,
        "Hs" => 108,
        "Mt" => 109,
        "Ds" => 110,
        "Rg" => 111,
        "Cn" => 112,
        "Nh" => 113,
        "Fl" => 114,
        "Mc" => 115,
        "Lv" => 116,
        "Ts" => 117,
        "Og" => 118,
        "D" => 1,
        "T" => 1,
        _ => 0,
    }
}
