//! Collection types for SDF V3000 format.
//!
//! Collections are used to group atoms for R-groups, atom lists, and other purposes.

/// Collection types as defined in V3000 format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CollectionType {
    /// Atom list (e.g., "NOT \[C,N,O\]" or "\[C,N\]").
    AtomList,
    /// R-group definition.
    RGroup,
    /// Highlight group.
    Highlight,
    /// External connection point.
    ExternalConnectionPoint,
}

impl CollectionType {
    /// Creates a CollectionType from a V3000 string.
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_uppercase().as_str() {
            "ATOMLIST" | "ALS" => Some(CollectionType::AtomList),
            "RGROUP" | "RGP" => Some(CollectionType::RGroup),
            "HILITE" | "HIGHLIGHT" => Some(CollectionType::Highlight),
            "EXTERNALCONNECTIONPOINT" | "ECP" => Some(CollectionType::ExternalConnectionPoint),
            _ => None,
        }
    }

    /// Converts to V3000 format string.
    pub fn to_v3000_str(&self) -> &'static str {
        match self {
            CollectionType::AtomList => "ATOMLIST",
            CollectionType::RGroup => "RGROUP",
            CollectionType::Highlight => "HILITE",
            CollectionType::ExternalConnectionPoint => "EXTERNALCONNECTIONPOINT",
        }
    }
}

/// Represents a collection in V3000 format.
#[derive(Debug, Clone, PartialEq)]
pub struct Collection {
    /// Type of collection.
    pub collection_type: CollectionType,

    /// Unique identifier for this collection.
    pub id: u32,

    /// Indices of atoms in this collection (0-based).
    pub atoms: Vec<usize>,

    /// Element symbols for atom lists (e.g., ["C", "N", "O"]).
    pub elements: Option<Vec<String>>,

    /// Whether the atom list is a NOT list (exclude these elements).
    pub is_not_list: bool,

    /// R-group number (for RGroup type).
    pub rgroup_number: Option<u8>,

    /// Label for the collection.
    pub label: Option<String>,
}

impl Collection {
    /// Creates a new collection with the given type and ID.
    pub fn new(collection_type: CollectionType, id: u32) -> Self {
        Self {
            collection_type,
            id,
            atoms: Vec::new(),
            elements: None,
            is_not_list: false,
            rgroup_number: None,
            label: None,
        }
    }

    /// Creates an atom list collection.
    pub fn atom_list(id: u32, atoms: Vec<usize>, elements: Vec<String>, is_not: bool) -> Self {
        Self {
            collection_type: CollectionType::AtomList,
            id,
            atoms,
            elements: Some(elements),
            is_not_list: is_not,
            rgroup_number: None,
            label: None,
        }
    }

    /// Creates an R-group collection.
    pub fn rgroup(id: u32, rgroup_number: u8, atoms: Vec<usize>) -> Self {
        Self {
            collection_type: CollectionType::RGroup,
            id,
            atoms,
            elements: None,
            is_not_list: false,
            rgroup_number: Some(rgroup_number),
            label: None,
        }
    }

    /// Returns true if this collection contains the given atom index.
    pub fn contains_atom(&self, atom_index: usize) -> bool {
        self.atoms.contains(&atom_index)
    }

    /// Returns the number of atoms in this collection.
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }
}

impl Default for Collection {
    fn default() -> Self {
        Self::new(CollectionType::AtomList, 0)
    }
}
