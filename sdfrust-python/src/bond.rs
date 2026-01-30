//! Python wrappers for Bond, BondOrder, and BondStereo.

use pyo3::prelude::*;
use sdfrust::{Bond, BondOrder, BondStereo};

/// Bond order types as defined in the SDF specification.
///
/// Available orders:
/// - SINGLE (1): Single bond
/// - DOUBLE (2): Double bond
/// - TRIPLE (3): Triple bond
/// - AROMATIC (4): Aromatic bond
/// - SINGLE_OR_DOUBLE (5): Single or double query bond
/// - SINGLE_OR_AROMATIC (6): Single or aromatic query bond
/// - DOUBLE_OR_AROMATIC (7): Double or aromatic query bond
/// - ANY (8): Any bond type
/// - COORDINATION (9): Dative/coordinate bond (V3000 only)
/// - HYDROGEN (10): Hydrogen bond (V3000 only)
#[pyclass(name = "BondOrder")]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct PyBondOrder {
    pub(crate) inner: BondOrder,
}

impl From<BondOrder> for PyBondOrder {
    fn from(order: BondOrder) -> Self {
        PyBondOrder { inner: order }
    }
}

impl From<PyBondOrder> for BondOrder {
    fn from(order: PyBondOrder) -> Self {
        order.inner
    }
}

#[pymethods]
impl PyBondOrder {
    /// Create a single bond order.
    #[staticmethod]
    pub fn single() -> Self {
        PyBondOrder {
            inner: BondOrder::Single,
        }
    }

    /// Create a double bond order.
    #[staticmethod]
    pub fn double() -> Self {
        PyBondOrder {
            inner: BondOrder::Double,
        }
    }

    /// Create a triple bond order.
    #[staticmethod]
    pub fn triple() -> Self {
        PyBondOrder {
            inner: BondOrder::Triple,
        }
    }

    /// Create an aromatic bond order.
    #[staticmethod]
    pub fn aromatic() -> Self {
        PyBondOrder {
            inner: BondOrder::Aromatic,
        }
    }

    /// Create a single-or-double query bond order.
    #[staticmethod]
    pub fn single_or_double() -> Self {
        PyBondOrder {
            inner: BondOrder::SingleOrDouble,
        }
    }

    /// Create a single-or-aromatic query bond order.
    #[staticmethod]
    pub fn single_or_aromatic() -> Self {
        PyBondOrder {
            inner: BondOrder::SingleOrAromatic,
        }
    }

    /// Create a double-or-aromatic query bond order.
    #[staticmethod]
    pub fn double_or_aromatic() -> Self {
        PyBondOrder {
            inner: BondOrder::DoubleOrAromatic,
        }
    }

    /// Create an any-type query bond order.
    #[staticmethod]
    pub fn any() -> Self {
        PyBondOrder {
            inner: BondOrder::Any,
        }
    }

    /// Create a coordination (dative) bond order (V3000 only).
    #[staticmethod]
    pub fn coordination() -> Self {
        PyBondOrder {
            inner: BondOrder::Coordination,
        }
    }

    /// Create a hydrogen bond order (V3000 only).
    #[staticmethod]
    pub fn hydrogen() -> Self {
        PyBondOrder {
            inner: BondOrder::Hydrogen,
        }
    }

    /// Create a BondOrder from an SDF bond type value (1-10).
    ///
    /// Args:
    ///     value: The SDF bond type value.
    ///
    /// Returns:
    ///     A BondOrder instance, or None if the value is invalid.
    #[staticmethod]
    pub fn from_sdf(value: u8) -> Option<Self> {
        BondOrder::from_sdf(value).map(|o| PyBondOrder { inner: o })
    }

    /// Convert to SDF bond type value (1-10).
    pub fn to_sdf(&self) -> u8 {
        self.inner.to_sdf()
    }

    /// Returns the nominal bond order as a float.
    ///
    /// For aromatic bonds, returns 1.5.
    pub fn order(&self) -> f64 {
        self.inner.order()
    }

    /// Returns True if this is an aromatic or aromatic-query bond.
    pub fn is_aromatic(&self) -> bool {
        matches!(
            self.inner,
            BondOrder::Aromatic | BondOrder::SingleOrAromatic | BondOrder::DoubleOrAromatic
        )
    }

    fn __repr__(&self) -> String {
        let name = match self.inner {
            BondOrder::Single => "SINGLE",
            BondOrder::Double => "DOUBLE",
            BondOrder::Triple => "TRIPLE",
            BondOrder::Aromatic => "AROMATIC",
            BondOrder::SingleOrDouble => "SINGLE_OR_DOUBLE",
            BondOrder::SingleOrAromatic => "SINGLE_OR_AROMATIC",
            BondOrder::DoubleOrAromatic => "DOUBLE_OR_AROMATIC",
            BondOrder::Any => "ANY",
            BondOrder::Coordination => "COORDINATION",
            BondOrder::Hydrogen => "HYDROGEN",
        };
        format!("BondOrder.{}", name)
    }

    fn __str__(&self) -> String {
        match self.inner {
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
        }
        .to_string()
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

    // Class attributes for easy access
    #[classattr]
    const SINGLE: u8 = 1;
    #[classattr]
    const DOUBLE: u8 = 2;
    #[classattr]
    const TRIPLE: u8 = 3;
    #[classattr]
    const AROMATIC: u8 = 4;
    #[classattr]
    const SINGLE_OR_DOUBLE: u8 = 5;
    #[classattr]
    const SINGLE_OR_AROMATIC: u8 = 6;
    #[classattr]
    const DOUBLE_OR_AROMATIC: u8 = 7;
    #[classattr]
    const ANY: u8 = 8;
    #[classattr]
    const COORDINATION: u8 = 9;
    #[classattr]
    const HYDROGEN: u8 = 10;
}

/// Bond stereochemistry types.
///
/// Available types:
/// - NONE (0): Not stereo (default)
/// - UP (1): Up (wedge)
/// - EITHER (4): Either (wiggly)
/// - DOWN (6): Down (dashed)
#[pyclass(name = "BondStereo")]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct PyBondStereo {
    pub(crate) inner: BondStereo,
}

impl From<BondStereo> for PyBondStereo {
    fn from(stereo: BondStereo) -> Self {
        PyBondStereo { inner: stereo }
    }
}

impl From<PyBondStereo> for BondStereo {
    fn from(stereo: PyBondStereo) -> Self {
        stereo.inner
    }
}

#[pymethods]
impl PyBondStereo {
    /// Create a non-stereo bond.
    #[staticmethod]
    pub fn none() -> Self {
        PyBondStereo {
            inner: BondStereo::None,
        }
    }

    /// Create an up (wedge) stereo bond.
    #[staticmethod]
    pub fn up() -> Self {
        PyBondStereo {
            inner: BondStereo::Up,
        }
    }

    /// Create an either (wiggly) stereo bond.
    #[staticmethod]
    pub fn either() -> Self {
        PyBondStereo {
            inner: BondStereo::Either,
        }
    }

    /// Create a down (dashed) stereo bond.
    #[staticmethod]
    pub fn down() -> Self {
        PyBondStereo {
            inner: BondStereo::Down,
        }
    }

    /// Create a BondStereo from an SDF stereo value.
    #[staticmethod]
    pub fn from_sdf(value: u8) -> Self {
        PyBondStereo {
            inner: BondStereo::from_sdf(value),
        }
    }

    /// Convert to SDF stereo value.
    pub fn to_sdf(&self) -> u8 {
        self.inner.to_sdf()
    }

    fn __repr__(&self) -> String {
        let name = match self.inner {
            BondStereo::None => "NONE",
            BondStereo::Up => "UP",
            BondStereo::Either => "EITHER",
            BondStereo::Down => "DOWN",
        };
        format!("BondStereo.{}", name)
    }

    fn __str__(&self) -> String {
        match self.inner {
            BondStereo::None => "none",
            BondStereo::Up => "up",
            BondStereo::Either => "either",
            BondStereo::Down => "down",
        }
        .to_string()
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

    // Class attributes
    #[classattr]
    const NONE: u8 = 0;
    #[classattr]
    const UP: u8 = 1;
    #[classattr]
    const EITHER: u8 = 4;
    #[classattr]
    const DOWN: u8 = 6;
}

/// A bond between two atoms in a molecule.
#[pyclass(name = "Bond")]
#[derive(Clone)]
pub struct PyBond {
    pub(crate) inner: Bond,
}

impl From<Bond> for PyBond {
    fn from(bond: Bond) -> Self {
        PyBond { inner: bond }
    }
}

impl From<&Bond> for PyBond {
    fn from(bond: &Bond) -> Self {
        PyBond {
            inner: bond.clone(),
        }
    }
}

impl From<PyBond> for Bond {
    fn from(bond: PyBond) -> Self {
        bond.inner
    }
}

#[pymethods]
impl PyBond {
    /// Create a new bond between two atoms.
    ///
    /// Args:
    ///     atom1: Index of the first atom (0-based).
    ///     atom2: Index of the second atom (0-based).
    ///     order: The bond order.
    ///
    /// Returns:
    ///     A new Bond instance.
    #[new]
    #[pyo3(signature = (atom1, atom2, order))]
    pub fn new(atom1: usize, atom2: usize, order: PyBondOrder) -> Self {
        PyBond {
            inner: Bond::new(atom1, atom2, order.inner),
        }
    }

    /// Create a new bond with stereochemistry.
    ///
    /// Args:
    ///     atom1: Index of the first atom (0-based).
    ///     atom2: Index of the second atom (0-based).
    ///     order: The bond order.
    ///     stereo: The bond stereochemistry.
    ///
    /// Returns:
    ///     A new Bond instance with stereochemistry.
    #[staticmethod]
    #[pyo3(signature = (atom1, atom2, order, stereo))]
    pub fn with_stereo(
        atom1: usize,
        atom2: usize,
        order: PyBondOrder,
        stereo: PyBondStereo,
    ) -> Self {
        PyBond {
            inner: Bond::with_stereo(atom1, atom2, order.inner, stereo.inner),
        }
    }

    /// Index of the first atom (0-based).
    #[getter]
    pub fn atom1(&self) -> usize {
        self.inner.atom1
    }

    /// Index of the second atom (0-based).
    #[getter]
    pub fn atom2(&self) -> usize {
        self.inner.atom2
    }

    /// The bond order.
    #[getter]
    pub fn order(&self) -> PyBondOrder {
        PyBondOrder {
            inner: self.inner.order,
        }
    }

    /// The bond stereochemistry.
    #[getter]
    pub fn stereo(&self) -> PyBondStereo {
        PyBondStereo {
            inner: self.inner.stereo,
        }
    }

    /// Bond topology (0 = either, 1 = ring, 2 = chain).
    #[getter]
    pub fn topology(&self) -> Option<u8> {
        self.inner.topology
    }

    /// Returns True if this bond involves the given atom index.
    pub fn contains_atom(&self, atom_index: usize) -> bool {
        self.inner.contains_atom(atom_index)
    }

    /// Returns the other atom in the bond given one atom index.
    ///
    /// Args:
    ///     atom_index: One of the atom indices in this bond.
    ///
    /// Returns:
    ///     The index of the other atom, or None if the given index is not part of this bond.
    pub fn other_atom(&self, atom_index: usize) -> Option<usize> {
        self.inner.other_atom(atom_index)
    }

    /// Returns True if this is an aromatic bond.
    pub fn is_aromatic(&self) -> bool {
        self.inner.is_aromatic()
    }

    fn __repr__(&self) -> String {
        format!(
            "Bond(atom1={}, atom2={}, order={}, stereo={})",
            self.inner.atom1,
            self.inner.atom2,
            PyBondOrder::from(self.inner.order).__str__(),
            PyBondStereo::from(self.inner.stereo).__str__()
        )
    }

    fn __str__(&self) -> String {
        format!(
            "{}-{} ({})",
            self.inner.atom1,
            self.inner.atom2,
            PyBondOrder::from(self.inner.order).__str__()
        )
    }

    fn __eq__(&self, other: &PyBond) -> bool {
        self.inner == other.inner
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.inner.atom1.hash(&mut hasher);
        self.inner.atom2.hash(&mut hasher);
        self.inner.order.hash(&mut hasher);
        self.inner.stereo.hash(&mut hasher);
        hasher.finish()
    }
}
