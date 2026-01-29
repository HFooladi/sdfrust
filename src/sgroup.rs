//! SGroup types for SDF V3000 format.
//!
//! SGroups (Sgroup or S-groups) are used to represent superatom abbreviations,
//! structural repeat units (polymers), data groups, and other structural features.

/// SGroup types as defined in V3000 format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SGroupType {
    /// Superatom abbreviation (e.g., Ph for phenyl, Me for methyl).
    Superatom,
    /// Multiple group for contracted representations.
    Multiple,
    /// Structural repeat unit (polymer).
    StructureRepeatUnit,
    /// Data group for attaching data to atoms.
    Data,
    /// Generic SGroup.
    Generic,
    /// Monomer unit.
    Monomer,
    /// Mer unit (for polymers).
    Mer,
    /// Copolymer unit.
    Copolymer,
    /// Component group.
    Component,
    /// Mixture group.
    Mixture,
    /// Formulation group.
    Formulation,
}

impl SGroupType {
    /// Creates an SGroupType from a V3000 string.
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_uppercase().as_str() {
            "SUP" => Some(SGroupType::Superatom),
            "MUL" => Some(SGroupType::Multiple),
            "SRU" => Some(SGroupType::StructureRepeatUnit),
            "DAT" => Some(SGroupType::Data),
            "GEN" => Some(SGroupType::Generic),
            "MON" => Some(SGroupType::Monomer),
            "MER" => Some(SGroupType::Mer),
            "COP" => Some(SGroupType::Copolymer),
            "COM" => Some(SGroupType::Component),
            "MIX" => Some(SGroupType::Mixture),
            "FOR" => Some(SGroupType::Formulation),
            _ => None,
        }
    }

    /// Converts to V3000 format string.
    pub fn to_v3000_str(&self) -> &'static str {
        match self {
            SGroupType::Superatom => "SUP",
            SGroupType::Multiple => "MUL",
            SGroupType::StructureRepeatUnit => "SRU",
            SGroupType::Data => "DAT",
            SGroupType::Generic => "GEN",
            SGroupType::Monomer => "MON",
            SGroupType::Mer => "MER",
            SGroupType::Copolymer => "COP",
            SGroupType::Component => "COM",
            SGroupType::Mixture => "MIX",
            SGroupType::Formulation => "FOR",
        }
    }
}

/// Represents an SGroup in V3000 format.
#[derive(Debug, Clone, PartialEq)]
pub struct SGroup {
    /// Unique identifier for this SGroup.
    pub id: u32,

    /// Type of SGroup.
    pub sgroup_type: SGroupType,

    /// Indices of atoms in this SGroup (0-based).
    pub atoms: Vec<usize>,

    /// Indices of bonds in this SGroup (0-based).
    pub bonds: Vec<usize>,

    /// Subscript label (e.g., "n" for polymers, abbreviation for superatoms).
    pub subscript: Option<String>,

    /// Superscript label.
    pub superscript: Option<String>,

    /// Parent SGroup ID (for nested SGroups).
    pub parent_id: Option<u32>,

    /// Connectivity type for SRU (hh = head-to-head, ht = head-to-tail, eu = either/unknown).
    pub connectivity: Option<String>,

    /// Label for the SGroup (display text).
    pub label: Option<String>,

    /// Data field name (for DAT type).
    pub field_name: Option<String>,

    /// Data field value (for DAT type).
    pub field_value: Option<String>,

    /// Crossing bonds for SRU.
    pub crossing_bonds: Vec<usize>,

    /// Bracket coordinates (for display).
    pub brackets: Vec<(f64, f64, f64, f64)>,
}

impl SGroup {
    /// Creates a new SGroup with the given ID and type.
    pub fn new(id: u32, sgroup_type: SGroupType) -> Self {
        Self {
            id,
            sgroup_type,
            atoms: Vec::new(),
            bonds: Vec::new(),
            subscript: None,
            superscript: None,
            parent_id: None,
            connectivity: None,
            label: None,
            field_name: None,
            field_value: None,
            crossing_bonds: Vec::new(),
            brackets: Vec::new(),
        }
    }

    /// Creates a new superatom SGroup.
    pub fn superatom(id: u32, label: &str, atoms: Vec<usize>) -> Self {
        let mut sg = Self::new(id, SGroupType::Superatom);
        sg.atoms = atoms;
        sg.label = Some(label.to_string());
        sg.subscript = Some(label.to_string());
        sg
    }

    /// Creates a new polymer SGroup.
    pub fn polymer(id: u32, subscript: &str, atoms: Vec<usize>) -> Self {
        let mut sg = Self::new(id, SGroupType::StructureRepeatUnit);
        sg.atoms = atoms;
        sg.subscript = Some(subscript.to_string());
        sg
    }

    /// Returns true if this SGroup contains the given atom index.
    pub fn contains_atom(&self, atom_index: usize) -> bool {
        self.atoms.contains(&atom_index)
    }

    /// Returns true if this SGroup contains the given bond index.
    pub fn contains_bond(&self, bond_index: usize) -> bool {
        self.bonds.contains(&bond_index)
    }

    /// Returns the number of atoms in this SGroup.
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Returns the number of bonds in this SGroup.
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }
}

impl Default for SGroup {
    fn default() -> Self {
        Self::new(0, SGroupType::Generic)
    }
}
