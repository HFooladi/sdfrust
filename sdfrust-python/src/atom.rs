//! Python wrapper for Atom struct.

use pyo3::prelude::*;
use sdfrust::Atom;

/// A Python-accessible wrapper for an atom in a molecule.
///
/// Atoms have 3D coordinates, an element symbol, and optional properties
/// like formal charge and mass difference.
#[pyclass(name = "Atom")]
#[derive(Clone)]
pub struct PyAtom {
    pub(crate) inner: Atom,
}

impl From<Atom> for PyAtom {
    fn from(atom: Atom) -> Self {
        PyAtom { inner: atom }
    }
}

impl From<&Atom> for PyAtom {
    fn from(atom: &Atom) -> Self {
        PyAtom {
            inner: atom.clone(),
        }
    }
}

impl From<PyAtom> for Atom {
    fn from(atom: PyAtom) -> Self {
        atom.inner
    }
}

#[pymethods]
impl PyAtom {
    /// Create a new atom with the given properties.
    ///
    /// Args:
    ///     index: The atom index (0-based).
    ///     element: The element symbol (e.g., "C", "N", "O").
    ///     x: The x coordinate in Angstroms.
    ///     y: The y coordinate in Angstroms.
    ///     z: The z coordinate in Angstroms.
    ///
    /// Returns:
    ///     A new Atom instance.
    #[new]
    #[pyo3(signature = (index, element, x, y, z))]
    pub fn new(index: usize, element: &str, x: f64, y: f64, z: f64) -> Self {
        PyAtom {
            inner: Atom::new(index, element, x, y, z),
        }
    }

    /// The atom index (0-based).
    #[getter]
    pub fn index(&self) -> usize {
        self.inner.index
    }

    /// The element symbol (e.g., "C", "N", "O").
    #[getter]
    pub fn element(&self) -> &str {
        &self.inner.element
    }

    /// The x coordinate in Angstroms.
    #[getter]
    pub fn x(&self) -> f64 {
        self.inner.x
    }

    /// The y coordinate in Angstroms.
    #[getter]
    pub fn y(&self) -> f64 {
        self.inner.y
    }

    /// The z coordinate in Angstroms.
    #[getter]
    pub fn z(&self) -> f64 {
        self.inner.z
    }

    /// The formal charge (-15 to +15, 0 = uncharged).
    #[getter]
    pub fn formal_charge(&self) -> i8 {
        self.inner.formal_charge
    }

    /// Set the formal charge.
    #[setter]
    pub fn set_formal_charge(&mut self, charge: i8) {
        self.inner.formal_charge = charge;
    }

    /// The mass difference from monoisotopic mass (-3 to +4).
    #[getter]
    pub fn mass_difference(&self) -> i8 {
        self.inner.mass_difference
    }

    /// Set the mass difference.
    #[setter]
    pub fn set_mass_difference(&mut self, diff: i8) {
        self.inner.mass_difference = diff;
    }

    /// Stereo parity (0 = not stereo, 1 = odd, 2 = even, 3 = either/unknown).
    #[getter]
    pub fn stereo_parity(&self) -> Option<u8> {
        self.inner.stereo_parity
    }

    /// Hydrogen count (0 = use default, 1 = H0, 2 = H1, etc.).
    #[getter]
    pub fn hydrogen_count(&self) -> Option<u8> {
        self.inner.hydrogen_count
    }

    /// Valence (0 = use default, 15 = zero valence).
    #[getter]
    pub fn valence(&self) -> Option<u8> {
        self.inner.valence
    }

    /// Returns the 3D coordinates as a tuple (x, y, z).
    pub fn coords(&self) -> (f64, f64, f64) {
        self.inner.coords()
    }

    /// Returns the distance to another atom in Angstroms.
    ///
    /// Args:
    ///     other: The other atom to calculate distance to.
    ///
    /// Returns:
    ///     The Euclidean distance in Angstroms.
    pub fn distance_to(&self, other: &PyAtom) -> f64 {
        self.inner.distance_to(&other.inner)
    }

    /// Returns True if this atom has a non-zero formal charge.
    pub fn is_charged(&self) -> bool {
        self.inner.is_charged()
    }

    fn __repr__(&self) -> String {
        format!(
            "Atom(index={}, element='{}', x={:.4}, y={:.4}, z={:.4}, charge={})",
            self.inner.index,
            self.inner.element,
            self.inner.x,
            self.inner.y,
            self.inner.z,
            self.inner.formal_charge
        )
    }

    fn __str__(&self) -> String {
        format!(
            "{}[{}] at ({:.4}, {:.4}, {:.4})",
            self.inner.element, self.inner.index, self.inner.x, self.inner.y, self.inner.z
        )
    }

    fn __eq__(&self, other: &PyAtom) -> bool {
        self.inner == other.inner
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.inner.index.hash(&mut hasher);
        self.inner.element.hash(&mut hasher);
        // Hash coordinates with fixed precision to avoid floating point issues
        ((self.inner.x * 10000.0) as i64).hash(&mut hasher);
        ((self.inner.y * 10000.0) as i64).hash(&mut hasher);
        ((self.inner.z * 10000.0) as i64).hash(&mut hasher);
        hasher.finish()
    }
}
