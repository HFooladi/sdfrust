//! Enhanced stereochemistry types for SDF V3000 format.

/// Stereo group types as defined in V3000 format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StereoGroupType {
    /// Absolute stereochemistry - the depicted configuration is the actual configuration.
    Absolute,
    /// OR group - one of the depicted configurations exists, but which one is unknown.
    Or,
    /// AND group - a mixture of both configurations exists.
    And,
}

impl StereoGroupType {
    /// Creates a StereoGroupType from a V3000 string.
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_uppercase().as_str() {
            "ABS" | "ABSOLUTE" => Some(StereoGroupType::Absolute),
            "OR" => Some(StereoGroupType::Or),
            "AND" | "&" => Some(StereoGroupType::And),
            _ => None,
        }
    }

    /// Converts to V3000 format string.
    pub fn to_v3000_str(&self) -> &'static str {
        match self {
            StereoGroupType::Absolute => "ABS",
            StereoGroupType::Or => "OR",
            StereoGroupType::And => "AND",
        }
    }
}

/// Represents an enhanced stereochemistry group in V3000 format.
#[derive(Debug, Clone, PartialEq)]
pub struct StereoGroup {
    /// The type of stereo group.
    pub group_type: StereoGroupType,

    /// Group number (for distinguishing multiple OR or AND groups).
    pub group_number: u32,

    /// Indices of atoms in this stereo group (0-based).
    pub atoms: Vec<usize>,
}

impl StereoGroup {
    /// Creates a new stereo group.
    pub fn new(group_type: StereoGroupType, group_number: u32, atoms: Vec<usize>) -> Self {
        Self {
            group_type,
            group_number,
            atoms,
        }
    }

    /// Creates an absolute stereo group.
    pub fn absolute(atoms: Vec<usize>) -> Self {
        Self::new(StereoGroupType::Absolute, 0, atoms)
    }

    /// Creates an OR stereo group.
    pub fn or_group(group_number: u32, atoms: Vec<usize>) -> Self {
        Self::new(StereoGroupType::Or, group_number, atoms)
    }

    /// Creates an AND stereo group.
    pub fn and_group(group_number: u32, atoms: Vec<usize>) -> Self {
        Self::new(StereoGroupType::And, group_number, atoms)
    }

    /// Returns true if this group contains the given atom index.
    pub fn contains_atom(&self, atom_index: usize) -> bool {
        self.atoms.contains(&atom_index)
    }

    /// Returns the number of atoms in this group.
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }
}

impl Default for StereoGroup {
    fn default() -> Self {
        Self {
            group_type: StereoGroupType::Absolute,
            group_number: 0,
            atoms: Vec::new(),
        }
    }
}
