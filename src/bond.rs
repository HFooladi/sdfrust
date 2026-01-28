/// Bond order types as defined in the SDF specification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    Single = 1,
    Double = 2,
    Triple = 3,
    Aromatic = 4,
    SingleOrDouble = 5,
    SingleOrAromatic = 6,
    DoubleOrAromatic = 7,
    Any = 8,
}

impl BondOrder {
    /// Creates a BondOrder from an SDF bond type value.
    pub fn from_sdf(value: u8) -> Option<Self> {
        match value {
            1 => Some(BondOrder::Single),
            2 => Some(BondOrder::Double),
            3 => Some(BondOrder::Triple),
            4 => Some(BondOrder::Aromatic),
            5 => Some(BondOrder::SingleOrDouble),
            6 => Some(BondOrder::SingleOrAromatic),
            7 => Some(BondOrder::DoubleOrAromatic),
            8 => Some(BondOrder::Any),
            _ => None,
        }
    }

    /// Converts to SDF bond type value.
    pub fn to_sdf(&self) -> u8 {
        *self as u8
    }

    /// Returns the nominal bond order as a float (for aromatic bonds, returns 1.5).
    pub fn order(&self) -> f64 {
        match self {
            BondOrder::Single => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Aromatic => 1.5,
            BondOrder::SingleOrDouble => 1.5,
            BondOrder::SingleOrAromatic => 1.25,
            BondOrder::DoubleOrAromatic => 1.75,
            BondOrder::Any => 1.0,
        }
    }
}

/// Bond stereochemistry types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondStereo {
    /// Not stereo (default).
    None = 0,
    /// Up (wedge).
    Up = 1,
    /// Either (wiggly).
    Either = 4,
    /// Down (dashed).
    Down = 6,
}

impl BondStereo {
    /// Creates a BondStereo from an SDF stereo value.
    pub fn from_sdf(value: u8) -> Self {
        match value {
            1 => BondStereo::Up,
            4 => BondStereo::Either,
            6 => BondStereo::Down,
            _ => BondStereo::None,
        }
    }

    /// Converts to SDF stereo value.
    pub fn to_sdf(&self) -> u8 {
        *self as u8
    }
}

/// Represents a bond between two atoms.
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    /// Index of the first atom (0-based).
    pub atom1: usize,

    /// Index of the second atom (0-based).
    pub atom2: usize,

    /// Bond order.
    pub order: BondOrder,

    /// Bond stereochemistry.
    pub stereo: BondStereo,

    /// Bond topology (0 = either, 1 = ring, 2 = chain).
    pub topology: Option<u8>,
}

impl Bond {
    /// Creates a new bond between two atoms.
    pub fn new(atom1: usize, atom2: usize, order: BondOrder) -> Self {
        Self {
            atom1,
            atom2,
            order,
            stereo: BondStereo::None,
            topology: None,
        }
    }

    /// Creates a new bond with stereochemistry.
    pub fn with_stereo(atom1: usize, atom2: usize, order: BondOrder, stereo: BondStereo) -> Self {
        Self {
            atom1,
            atom2,
            order,
            stereo,
            topology: None,
        }
    }

    /// Returns true if this bond involves the given atom index.
    pub fn contains_atom(&self, atom_index: usize) -> bool {
        self.atom1 == atom_index || self.atom2 == atom_index
    }

    /// Returns the other atom in the bond given one atom index.
    /// Returns None if the given index is not part of this bond.
    pub fn other_atom(&self, atom_index: usize) -> Option<usize> {
        if self.atom1 == atom_index {
            Some(self.atom2)
        } else if self.atom2 == atom_index {
            Some(self.atom1)
        } else {
            None
        }
    }

    /// Returns true if this is an aromatic bond.
    pub fn is_aromatic(&self) -> bool {
        matches!(
            self.order,
            BondOrder::Aromatic | BondOrder::SingleOrAromatic | BondOrder::DoubleOrAromatic
        )
    }
}
