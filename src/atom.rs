/// Represents an atom in a molecule.
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    /// Atom index (0-based).
    pub index: usize,

    /// Element symbol (e.g., "C", "N", "O").
    pub element: String,

    /// X coordinate in Angstroms.
    pub x: f64,

    /// Y coordinate in Angstroms.
    pub y: f64,

    /// Z coordinate in Angstroms.
    pub z: f64,

    /// Formal charge (-15 to +15, 0 = uncharged).
    pub formal_charge: i8,

    /// Mass difference from monoisotopic mass (-3 to +4).
    pub mass_difference: i8,

    /// Stereo parity (0 = not stereo, 1 = odd, 2 = even, 3 = either/unknown).
    pub stereo_parity: Option<u8>,

    /// Hydrogen count (0 = use default, 1 = H0, 2 = H1, etc.).
    pub hydrogen_count: Option<u8>,

    /// Valence (0 = use default, 15 = zero valence).
    pub valence: Option<u8>,
}

impl Atom {
    /// Creates a new atom with the given element and coordinates.
    pub fn new(index: usize, element: &str, x: f64, y: f64, z: f64) -> Self {
        Self {
            index,
            element: element.to_string(),
            x,
            y,
            z,
            formal_charge: 0,
            mass_difference: 0,
            stereo_parity: None,
            hydrogen_count: None,
            valence: None,
        }
    }

    /// Returns the 3D coordinates as a tuple.
    pub fn coords(&self) -> (f64, f64, f64) {
        (self.x, self.y, self.z)
    }

    /// Returns the distance to another atom.
    pub fn distance_to(&self, other: &Atom) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Returns true if this atom has a non-zero formal charge.
    pub fn is_charged(&self) -> bool {
        self.formal_charge != 0
    }
}

impl Default for Atom {
    fn default() -> Self {
        Self {
            index: 0,
            element: String::new(),
            x: 0.0,
            y: 0.0,
            z: 0.0,
            formal_charge: 0,
            mass_difference: 0,
            stereo_parity: None,
            hydrogen_count: None,
            valence: None,
        }
    }
}
