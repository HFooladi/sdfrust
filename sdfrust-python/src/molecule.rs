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

    // ============================================================
    // Bond Inference
    // ============================================================

    /// Infer single bonds from 3D coordinates and covalent radii.
    ///
    /// Two atoms are bonded if their distance is within the sum of their
    /// covalent radii plus a tolerance. All inferred bonds are single bonds.
    /// Existing bonds are cleared before inference.
    ///
    /// This is useful for XYZ files which contain no bond information.
    ///
    /// Args:
    ///     tolerance: Optional tolerance in Angstroms (default: 0.45).
    ///
    /// Raises:
    ///     ValueError: If any atom has an unknown element.
    ///
    /// Example:
    ///     >>> mol = sdfrust.parse_xyz_file("water.xyz")
    ///     >>> mol.infer_bonds()
    ///     >>> print(mol.num_bonds)  # 2
    #[pyo3(signature = (tolerance=None))]
    pub fn infer_bonds(&mut self, tolerance: Option<f64>) -> PyResult<()> {
        self.inner
            .infer_bonds(tolerance)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

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

    // ============================================================
    // Chemical Perception & ML Features
    // ============================================================

    /// Compute the Smallest Set of Smallest Rings (SSSR).
    ///
    /// Returns a list of rings, where each ring is a list of atom indices.
    ///
    /// Returns:
    ///     List of lists of atom indices, one per ring.
    ///
    /// Example:
    ///     >>> mol = sdfrust.parse_sdf_file("benzene.sdf")
    ///     >>> rings = mol.sssr()
    ///     >>> print(len(rings))  # 1
    pub fn sssr(&self) -> Vec<Vec<usize>> {
        sdfrust::descriptors::sssr(&self.inner)
            .into_iter()
            .map(|r| r.atoms)
            .collect()
    }

    /// Get ring sizes that an atom participates in.
    ///
    /// Args:
    ///     atom_index: The atom index.
    ///
    /// Returns:
    ///     List of ring sizes (e.g., [5, 6] if in both a 5- and 6-membered ring).
    pub fn ring_sizes(&self, atom_index: usize) -> Vec<usize> {
        sdfrust::descriptors::ring_sizes(&self.inner, atom_index)
    }

    /// Get the smallest ring size for an atom.
    ///
    /// Args:
    ///     atom_index: The atom index.
    ///
    /// Returns:
    ///     Smallest ring size, or None if the atom is not in any ring.
    pub fn smallest_ring_size(&self, atom_index: usize) -> Option<usize> {
        sdfrust::descriptors::smallest_ring_size(&self.inner, atom_index)
    }

    /// Get the hybridization of an atom.
    ///
    /// Args:
    ///     atom_index: The atom index.
    ///
    /// Returns:
    ///     Hybridization string: "S", "SP", "SP2", "SP3", "SP3D", "SP3D2", or "Other".
    pub fn atom_hybridization(&self, atom_index: usize) -> String {
        format!(
            "{:?}",
            sdfrust::descriptors::atom_hybridization(&self.inner, atom_index)
        )
    }

    /// Get hybridizations for all atoms.
    ///
    /// Returns:
    ///     List of hybridization strings.
    pub fn all_hybridizations(&self) -> Vec<String> {
        sdfrust::descriptors::all_hybridizations(&self.inner)
            .into_iter()
            .map(|h| format!("{:?}", h))
            .collect()
    }

    /// Check if an atom is aromatic (Hückel 4n+2 rule or file annotation).
    ///
    /// Args:
    ///     atom_index: The atom index.
    ///
    /// Returns:
    ///     True if the atom is aromatic.
    pub fn is_aromatic_atom(&self, atom_index: usize) -> bool {
        sdfrust::descriptors::is_aromatic_atom(&self.inner, atom_index)
    }

    /// Check if a bond is aromatic.
    ///
    /// Args:
    ///     bond_index: The bond index.
    ///
    /// Returns:
    ///     True if the bond is aromatic.
    pub fn is_aromatic_bond(&self, bond_index: usize) -> bool {
        sdfrust::descriptors::is_aromatic_bond(&self.inner, bond_index)
    }

    /// Get aromaticity for all atoms.
    ///
    /// Returns:
    ///     List of booleans, one per atom.
    pub fn all_aromatic_atoms(&self) -> Vec<bool> {
        sdfrust::descriptors::all_aromatic_atoms(&self.inner)
    }

    /// Get aromaticity for all bonds.
    ///
    /// Returns:
    ///     List of booleans, one per bond.
    pub fn all_aromatic_bonds(&self) -> Vec<bool> {
        sdfrust::descriptors::all_aromatic_bonds(&self.inner)
    }

    /// Check if a bond is conjugated.
    ///
    /// Args:
    ///     bond_index: The bond index.
    ///
    /// Returns:
    ///     True if the bond is conjugated.
    pub fn is_conjugated_bond(&self, bond_index: usize) -> bool {
        sdfrust::descriptors::is_conjugated_bond(&self.inner, bond_index)
    }

    /// Get conjugation status for all bonds.
    ///
    /// Returns:
    ///     List of booleans, one per bond.
    pub fn all_conjugated_bonds(&self) -> Vec<bool> {
        sdfrust::descriptors::all_conjugated_bonds(&self.inner)
    }

    /// Get the degree (number of bonds) for an atom.
    ///
    /// Args:
    ///     atom_index: The atom index.
    ///
    /// Returns:
    ///     Number of bonds connected to this atom.
    pub fn atom_degree(&self, atom_index: usize) -> usize {
        sdfrust::descriptors::atom_degree(&self.inner, atom_index)
    }

    /// Get the total hydrogen count for an atom (implicit + explicit).
    ///
    /// Args:
    ///     atom_index: The atom index.
    ///
    /// Returns:
    ///     Total number of hydrogens on this atom.
    pub fn total_hydrogen_count(&self, atom_index: usize) -> u8 {
        sdfrust::descriptors::total_hydrogen_count(&self.inner, atom_index)
    }

    /// Get implicit hydrogen count for an atom.
    ///
    /// Args:
    ///     atom_index: The atom index.
    ///
    /// Returns:
    ///     Number of implicit hydrogens.
    pub fn implicit_hydrogen_count(&self, atom_index: usize) -> u8 {
        sdfrust::descriptors::implicit_hydrogen_count(&self.inner, atom_index)
    }

    /// Compute Gasteiger partial charges for all atoms.
    ///
    /// Uses the PEOE (Partial Equalization of Orbital Electronegativity) algorithm
    /// with 6 iterations and 0.5 damping factor.
    ///
    /// Returns:
    ///     List of partial charges (float), one per atom.
    ///
    /// Example:
    ///     >>> mol = sdfrust.parse_sdf_file("aspirin.sdf")
    ///     >>> charges = mol.gasteiger_charges()
    ///     >>> print(f"O charge: {charges[0]:.3f}")
    pub fn gasteiger_charges(&self) -> Vec<f64> {
        sdfrust::descriptors::gasteiger_charges(&self.inner)
    }

    /// Compute Gasteiger charges with custom parameters.
    ///
    /// Args:
    ///     max_iter: Maximum number of iterations (default: 6).
    ///     damping: Damping factor (default: 0.5).
    ///
    /// Returns:
    ///     List of partial charges (float), one per atom.
    #[pyo3(signature = (max_iter=6, damping=0.5))]
    pub fn gasteiger_charges_with_params(&self, max_iter: usize, damping: f64) -> Vec<f64> {
        sdfrust::descriptors::gasteiger_charges_with_params(&self.inner, max_iter, damping)
    }

    /// Compute OGB-compatible atom features.
    ///
    /// Returns a 2D list of shape [N, 9] where N is the number of atoms.
    /// Features: atomic_number, chirality, degree, formal_charge+5,
    /// num_hs, radical, hybridization, is_aromatic, is_in_ring.
    ///
    /// Returns:
    ///     List of lists of integers (9 features per atom).
    ///
    /// Example:
    ///     >>> feats = mol.ogb_atom_features()
    ///     >>> print(feats[0])  # [6, 0, 3, 5, 0, 0, 2, 1, 1] for aromatic C
    pub fn ogb_atom_features(&self) -> Vec<Vec<i32>> {
        sdfrust::featurize::ogb::ogb_atom_features(&self.inner).features
    }

    /// Compute OGB-compatible bond features.
    ///
    /// Returns a 2D list of shape [E, 3] where E is the number of bonds.
    /// Features: bond_type (0-3), stereo (0-3), is_conjugated (0/1).
    ///
    /// Returns:
    ///     List of lists of integers (3 features per bond).
    pub fn ogb_bond_features(&self) -> Vec<Vec<i32>> {
        sdfrust::featurize::ogb::ogb_bond_features(&self.inner).features
    }

    /// Compute complete OGB graph features with directed edge index.
    ///
    /// Returns a dictionary with:
    ///   - "atom_features": [N, 9] atom feature matrix
    ///   - "bond_features": [2E, 3] directed bond feature matrix
    ///   - "edge_src": source atoms for directed edges
    ///   - "edge_dst": destination atoms for directed edges
    ///
    /// Returns:
    ///     Dictionary of graph features matching PyTorch Geometric format.
    pub fn ogb_graph_features(&self, py: Python<'_>) -> PyResult<PyObject> {
        let graph = sdfrust::featurize::ogb::ogb_graph_features(&self.inner);
        let dict = pyo3::types::PyDict::new(py);
        dict.set_item("atom_features", &graph.atom_features.features)?;
        dict.set_item("bond_features", &graph.bond_features.features)?;
        dict.set_item("edge_src", &graph.edge_src)?;
        dict.set_item("edge_dst", &graph.edge_dst)?;
        dict.set_item("num_atoms", graph.atom_features.num_atoms)?;
        dict.set_item("num_bonds", graph.bond_features.num_bonds)?;
        Ok(dict.into())
    }

    /// Compute ECFP/Morgan fingerprint as a bit vector.
    ///
    /// Args:
    ///     radius: Fingerprint radius (2 for ECFP4, 3 for ECFP6). Default: 2.
    ///     n_bits: Length of the bit vector. Default: 2048.
    ///
    /// Returns:
    ///     List of booleans representing the fingerprint bit vector.
    ///
    /// Example:
    ///     >>> fp = mol.ecfp(radius=2, n_bits=2048)
    ///     >>> print(sum(fp))  # Number of set bits
    #[pyo3(signature = (radius=2, n_bits=2048))]
    pub fn ecfp(&self, radius: usize, n_bits: usize) -> Vec<bool> {
        sdfrust::fingerprints::ecfp::ecfp(&self.inner, radius, n_bits).bits
    }

    /// Compute ECFP/Morgan fingerprint as on-bit indices.
    ///
    /// Args:
    ///     radius: Fingerprint radius. Default: 2.
    ///     n_bits: Length of the bit vector. Default: 2048.
    ///
    /// Returns:
    ///     List of indices where bits are set.
    #[pyo3(signature = (radius=2, n_bits=2048))]
    pub fn ecfp_on_bits(&self, radius: usize, n_bits: usize) -> Vec<usize> {
        sdfrust::fingerprints::ecfp::ecfp(&self.inner, radius, n_bits).on_bits()
    }

    /// Compute Tanimoto similarity between two molecules using ECFP.
    ///
    /// Args:
    ///     other: The other molecule.
    ///     radius: Fingerprint radius. Default: 2.
    ///     n_bits: Length of the bit vector. Default: 2048.
    ///
    /// Returns:
    ///     Tanimoto similarity coefficient in [0, 1].
    #[pyo3(signature = (other, radius=2, n_bits=2048))]
    pub fn tanimoto_similarity(&self, other: &PyMolecule, radius: usize, n_bits: usize) -> f64 {
        let fp1 = sdfrust::fingerprints::ecfp::ecfp(&self.inner, radius, n_bits);
        let fp2 = sdfrust::fingerprints::ecfp::ecfp(&other.inner, radius, n_bits);
        fp1.tanimoto(&fp2)
    }

    /// Compute ECFP count fingerprint (hash → count map).
    ///
    /// Args:
    ///     radius: Fingerprint radius. Default: 2.
    ///
    /// Returns:
    ///     Dictionary mapping feature hash (int) to count (int).
    #[pyo3(signature = (radius=2))]
    pub fn ecfp_counts(&self, radius: usize) -> HashMap<u32, u32> {
        sdfrust::fingerprints::ecfp::ecfp_counts(&self.inner, radius).counts
    }

    // ============================================================
    // NumPy ML features (when numpy feature is enabled)
    // ============================================================

    #[cfg(feature = "numpy")]
    /// Get Gasteiger charges as a NumPy array.
    ///
    /// Returns:
    ///     NumPy array of shape (N,) with partial charges.
    pub fn get_gasteiger_charges_array<'py>(
        &self,
        py: Python<'py>,
    ) -> Bound<'py, numpy::PyArray1<f64>> {
        let charges = sdfrust::descriptors::gasteiger_charges(&self.inner);
        numpy::PyArray1::from_vec(py, charges)
    }

    #[cfg(feature = "numpy")]
    /// Get OGB atom features as a NumPy array.
    ///
    /// Returns:
    ///     NumPy array of shape (N, 9) with integer features.
    pub fn get_ogb_atom_features_array<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, numpy::PyArray2<i32>>> {
        let feats = sdfrust::featurize::ogb::ogb_atom_features(&self.inner);
        numpy::PyArray2::from_vec2(py, &feats.features)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[cfg(feature = "numpy")]
    /// Get OGB bond features as a NumPy array.
    ///
    /// Returns:
    ///     NumPy array of shape (E, 3) with integer features.
    pub fn get_ogb_bond_features_array<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, numpy::PyArray2<i32>>> {
        let feats = sdfrust::featurize::ogb::ogb_bond_features(&self.inner);
        numpy::PyArray2::from_vec2(py, &feats.features)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[cfg(feature = "numpy")]
    /// Get ECFP fingerprint as a NumPy boolean array.
    ///
    /// Args:
    ///     radius: Fingerprint radius. Default: 2.
    ///     n_bits: Length of the bit vector. Default: 2048.
    ///
    /// Returns:
    ///     NumPy boolean array of shape (n_bits,).
    #[pyo3(signature = (radius=2, n_bits=2048))]
    pub fn get_ecfp_array<'py>(
        &self,
        py: Python<'py>,
        radius: usize,
        n_bits: usize,
    ) -> Bound<'py, numpy::PyArray1<bool>> {
        let fp = sdfrust::fingerprints::ecfp::ecfp(&self.inner, radius, n_bits);
        numpy::PyArray1::from_vec(py, fp.bits)
    }

    // ============================================================
    // Geometry-based ML features (when geometry feature is enabled)
    // ============================================================

    #[cfg(feature = "geometry")]
    /// Compute cutoff-based neighbor list for 3D GNNs.
    ///
    /// Returns a dictionary with:
    ///   - "edge_src": source atom indices
    ///   - "edge_dst": destination atom indices
    ///   - "distances": pairwise distances
    ///   - "num_edges": number of directed edges
    ///
    /// Args:
    ///     cutoff: Distance cutoff in Angstroms.
    ///
    /// Returns:
    ///     Dictionary with edge index and distances.
    pub fn neighbor_list(&self, cutoff: f64, py: Python<'_>) -> PyResult<PyObject> {
        let nl = sdfrust::geometry::neighbor_list(&self.inner, cutoff);
        let dict = pyo3::types::PyDict::new(py);
        dict.set_item("edge_src", &nl.edge_src)?;
        dict.set_item("edge_dst", &nl.edge_dst)?;
        dict.set_item("distances", &nl.distances)?;
        dict.set_item("num_edges", nl.num_edges())?;
        Ok(dict.into())
    }

    #[cfg(feature = "geometry")]
    /// Compute a bond angle (in radians) for atoms i-j-k.
    ///
    /// Args:
    ///     i: First atom index.
    ///     j: Central atom index.
    ///     k: Third atom index.
    ///
    /// Returns:
    ///     Angle in radians, or None if undefined.
    pub fn bond_angle(&self, i: usize, j: usize, k: usize) -> Option<f64> {
        sdfrust::geometry::bond_angle(&self.inner, i, j, k)
    }

    #[cfg(feature = "geometry")]
    /// Compute a dihedral angle (in radians) for atoms i-j-k-l.
    ///
    /// Args:
    ///     i: First atom index.
    ///     j: Second atom index.
    ///     k: Third atom index.
    ///     l: Fourth atom index.
    ///
    /// Returns:
    ///     Dihedral angle in radians (-pi to pi), or None if undefined.
    pub fn dihedral_angle(&self, i: usize, j: usize, k: usize, l: usize) -> Option<f64> {
        sdfrust::geometry::dihedral_angle(&self.inner, i, j, k, l)
    }

    #[cfg(feature = "geometry")]
    /// Compute all bond angles in the molecule.
    ///
    /// Returns a dictionary with:
    ///   - "triplets": list of [i, j, k] triplets (j is central atom)
    ///   - "angles": list of angles in radians
    ///
    /// Returns:
    ///     Dictionary with triplet indices and angle values.
    pub fn all_bond_angles(&self, py: Python<'_>) -> PyResult<PyObject> {
        let result = sdfrust::geometry::all_bond_angles(&self.inner);
        let dict = pyo3::types::PyDict::new(py);
        let triplets: Vec<Vec<usize>> = result.triplets.iter().map(|t| t.to_vec()).collect();
        dict.set_item("triplets", triplets)?;
        dict.set_item("angles", &result.angles)?;
        Ok(dict.into())
    }

    #[cfg(feature = "geometry")]
    /// Compute all dihedral angles in the molecule.
    ///
    /// Returns a dictionary with:
    ///   - "quadruplets": list of [i, j, k, l] quadruplets
    ///   - "angles": list of dihedral angles in radians (-pi to pi)
    ///
    /// Returns:
    ///     Dictionary with quadruplet indices and angle values.
    pub fn all_dihedral_angles(&self, py: Python<'_>) -> PyResult<PyObject> {
        let result = sdfrust::geometry::all_dihedral_angles(&self.inner);
        let dict = pyo3::types::PyDict::new(py);
        let quads: Vec<Vec<usize>> = result.quadruplets.iter().map(|q| q.to_vec()).collect();
        dict.set_item("quadruplets", quads)?;
        dict.set_item("angles", &result.angles)?;
        Ok(dict.into())
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
