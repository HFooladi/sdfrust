//! SDF V3000 format parser.
//!
//! V3000 is an extended molfile format that supports:
//! - Molecules with >999 atoms/bonds
//! - Variable-width space-separated fields
//! - Enhanced stereochemistry
//! - SGroups and collections

use std::collections::HashMap;
use std::io::BufRead;

use crate::atom::Atom;
use crate::bond::{Bond, BondOrder, BondStereo};
use crate::collection::{Collection, CollectionType};
use crate::error::{Result, SdfError};
use crate::molecule::{Molecule, SdfFormat};
use crate::sgroup::{SGroup, SGroupType};
use crate::stereogroup::{StereoGroup, StereoGroupType};

const V30_PREFIX: &str = "M  V30 ";

/// SDF V3000 format parser.
pub struct SdfV3000Parser<R> {
    reader: R,
    line_number: usize,
    current_line: String,
    peeked_line: Option<String>,
}

impl<R: BufRead> SdfV3000Parser<R> {
    /// Creates a new V3000 parser from a buffered reader.
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_number: 0,
            current_line: String::new(),
            peeked_line: None,
        }
    }

    /// Reads the next line from the input.
    fn read_line(&mut self) -> Result<bool> {
        if let Some(line) = self.peeked_line.take() {
            self.current_line = line;
            return Ok(true);
        }

        self.current_line.clear();
        let bytes_read = self.reader.read_line(&mut self.current_line)?;
        if bytes_read > 0 {
            self.line_number += 1;
            // Remove trailing newline
            if self.current_line.ends_with('\n') {
                self.current_line.pop();
                if self.current_line.ends_with('\r') {
                    self.current_line.pop();
                }
            }
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Reads a V3000 line, handling line continuations (lines ending with -).
    fn read_v30_line(&mut self) -> Result<Option<String>> {
        let mut result = String::new();

        loop {
            if !self.read_line()? {
                if result.is_empty() {
                    return Ok(None);
                }
                break;
            }

            // Check for V30 prefix
            if !self.current_line.starts_with(V30_PREFIX) {
                // Not a V30 line, put it back
                self.peeked_line = Some(self.current_line.clone());
                if result.is_empty() {
                    return Ok(None);
                }
                break;
            }

            // Extract content after prefix
            let content = &self.current_line[V30_PREFIX.len()..];

            // Check for line continuation
            if let Some(stripped) = content.strip_suffix('-') {
                result.push_str(stripped);
            } else {
                result.push_str(content);
                break;
            }
        }

        Ok(Some(result))
    }

    /// Parses a single molecule from the input.
    /// The header (name, program line, comment, counts line) must already be read.
    /// This method starts parsing from the CTAB block.
    pub fn parse_molecule_body(
        &mut self,
        name: String,
        program_line: Option<String>,
        comment: Option<String>,
    ) -> Result<Molecule> {
        // Parse BEGIN CTAB
        let line = self
            .read_v30_line()?
            .ok_or_else(|| SdfError::InvalidV3000Block("Expected BEGIN CTAB".to_string()))?;

        if !line.starts_with("BEGIN CTAB") {
            return Err(SdfError::InvalidV3000Block(format!(
                "Expected BEGIN CTAB, got: {}",
                line
            )));
        }

        // Parse COUNTS line
        let counts_line = self
            .read_v30_line()?
            .ok_or_else(|| SdfError::InvalidV3000Block("Expected COUNTS line".to_string()))?;

        let (atom_count, bond_count) = self.parse_counts(&counts_line)?;

        // Parse atom block
        let (atoms, id_to_index) = self.parse_atom_block(atom_count)?;

        // Parse bond block
        let bonds = self.parse_bond_block(bond_count, &id_to_index)?;

        // Parse optional sections (SGROUP, COLLECTION, etc.)
        let mut stereogroups = Vec::new();
        let mut sgroups = Vec::new();
        let mut collections = Vec::new();

        loop {
            let Some(line) = self.read_v30_line()? else {
                break;
            };

            if line.starts_with("END CTAB") {
                break;
            } else if line.starts_with("BEGIN SGROUP") {
                sgroups = self.parse_sgroup_block(&id_to_index)?;
            } else if line.starts_with("BEGIN COLLECTION") {
                collections = self.parse_collection_block(&id_to_index)?;
            } else if line.starts_with("BEGIN STEREO") {
                stereogroups = self.parse_stereo_block(&id_to_index)?;
            }
            // Skip other blocks we don't handle
        }

        // Read until M  END or $$$$
        self.read_until_end()?;

        Ok(Molecule {
            name,
            program_line,
            comment,
            atoms,
            bonds,
            properties: HashMap::new(),
            format_version: SdfFormat::V3000,
            stereogroups,
            sgroups,
            collections,
        })
    }

    /// Parses the COUNTS line.
    fn parse_counts(&self, line: &str) -> Result<(usize, usize)> {
        // Format: COUNTS na nb nsg n3d chiral [regno]
        if !line.starts_with("COUNTS ") {
            return Err(SdfError::InvalidV3000Block(format!(
                "Expected COUNTS, got: {}",
                line
            )));
        }

        let parts: Vec<&str> = line[7..].split_whitespace().collect();
        if parts.len() < 2 {
            return Err(SdfError::InvalidV3000Block(format!(
                "Invalid COUNTS line: {}",
                line
            )));
        }

        let atom_count: usize = parts[0].parse().map_err(|_| {
            SdfError::InvalidV3000Block(format!("Invalid atom count: {}", parts[0]))
        })?;

        let bond_count: usize = parts[1].parse().map_err(|_| {
            SdfError::InvalidV3000Block(format!("Invalid bond count: {}", parts[1]))
        })?;

        Ok((atom_count, bond_count))
    }

    /// Parses the atom block.
    fn parse_atom_block(&mut self, count: usize) -> Result<(Vec<Atom>, HashMap<u32, usize>)> {
        // If count is 0, the ATOM block may not be present
        if count == 0 {
            return Ok((Vec::new(), HashMap::new()));
        }

        // Read BEGIN ATOM
        let line = self
            .read_v30_line()?
            .ok_or_else(|| SdfError::InvalidV3000Block("Expected BEGIN ATOM".to_string()))?;

        if !line.starts_with("BEGIN ATOM") {
            return Err(SdfError::InvalidV3000Block(format!(
                "Expected BEGIN ATOM, got: {}",
                line
            )));
        }

        let mut atoms = Vec::with_capacity(count);
        let mut id_to_index: HashMap<u32, usize> = HashMap::with_capacity(count);

        for i in 0..count {
            let atom_line = self.read_v30_line()?.ok_or(SdfError::AtomCountMismatch {
                expected: count,
                found: i,
            })?;

            let (atom, v3000_id) = self.parse_atom_line(&atom_line, i)?;
            id_to_index.insert(v3000_id, i);
            atoms.push(atom);
        }

        // Read END ATOM
        let end_line = self
            .read_v30_line()?
            .ok_or_else(|| SdfError::InvalidV3000Block("Expected END ATOM".to_string()))?;

        if !end_line.starts_with("END ATOM") {
            return Err(SdfError::InvalidV3000Block(format!(
                "Expected END ATOM, got: {}",
                end_line
            )));
        }

        Ok((atoms, id_to_index))
    }

    /// Parses a single atom line.
    fn parse_atom_line(&self, line: &str, index: usize) -> Result<(Atom, u32)> {
        // Format: index type x y z aamap [key=value...]
        // Example: 1 C 0.0000 0.0000 0.0000 0 CHG=-1 RAD=2
        let parts: Vec<&str> = line.split_whitespace().collect();

        if parts.len() < 5 {
            return Err(SdfError::InvalidV3000AtomLine {
                line: self.line_number,
                message: format!("Not enough fields: {}", line),
            });
        }

        let v3000_id: u32 = parts[0]
            .parse()
            .map_err(|_| SdfError::InvalidV3000AtomLine {
                line: self.line_number,
                message: format!("Invalid atom ID: {}", parts[0]),
            })?;

        let element = parts[1].to_string();

        let x: f64 = parts[2]
            .parse()
            .map_err(|_| SdfError::InvalidV3000AtomLine {
                line: self.line_number,
                message: format!("Invalid x coordinate: {}", parts[2]),
            })?;

        let y: f64 = parts[3]
            .parse()
            .map_err(|_| SdfError::InvalidV3000AtomLine {
                line: self.line_number,
                message: format!("Invalid y coordinate: {}", parts[3]),
            })?;

        let z: f64 = parts[4]
            .parse()
            .map_err(|_| SdfError::InvalidV3000AtomLine {
                line: self.line_number,
                message: format!("Invalid z coordinate: {}", parts[4]),
            })?;

        // Parse aamap (atom-atom mapping, 0 = not mapped)
        let atom_atom_mapping: Option<u32> = if parts.len() > 5 {
            let val: u32 = parts[5].parse().unwrap_or(0);
            if val > 0 { Some(val) } else { None }
        } else {
            None
        };

        // Parse key=value pairs
        let mut formal_charge: i8 = 0;
        let mut mass_difference: i8 = 0;
        let mut radical: Option<u8> = None;
        let mut valence: Option<u8> = None;
        let mut hydrogen_count: Option<u8> = None;
        let mut stereo_parity: Option<u8> = None;
        let mut rgroup_label: Option<u8> = None;

        for part in parts.iter().skip(6) {
            if let Some((key, value)) = part.split_once('=') {
                match key.to_uppercase().as_str() {
                    "CHG" => {
                        formal_charge = value.parse().unwrap_or(0);
                    }
                    "RAD" => {
                        radical = value.parse().ok();
                    }
                    "MASS" => {
                        mass_difference = value.parse().unwrap_or(0);
                    }
                    "VAL" => {
                        valence = value.parse().ok();
                    }
                    "HCOUNT" => {
                        hydrogen_count = value.parse().ok();
                    }
                    "CFG" => {
                        stereo_parity = value.parse().ok();
                    }
                    "RGROUPS" => {
                        // Parse R-group attachment, format: (n rg1 rg2...)
                        if value.starts_with('(') && value.ends_with(')') {
                            let inner = &value[1..value.len() - 1];
                            let nums: Vec<&str> = inner.split_whitespace().collect();
                            if nums.len() >= 2 {
                                rgroup_label = nums[1].parse().ok();
                            }
                        }
                    }
                    _ => {} // Ignore unknown keys
                }
            }
        }

        let atom = Atom {
            index,
            element,
            x,
            y,
            z,
            formal_charge,
            mass_difference,
            stereo_parity,
            hydrogen_count,
            valence,
            v3000_id: Some(v3000_id),
            atom_atom_mapping,
            rgroup_label,
            radical,
        };

        Ok((atom, v3000_id))
    }

    /// Parses the bond block.
    fn parse_bond_block(
        &mut self,
        count: usize,
        id_to_index: &HashMap<u32, usize>,
    ) -> Result<Vec<Bond>> {
        // If count is 0, the BOND block may not be present
        if count == 0 {
            return Ok(Vec::new());
        }

        // Read BEGIN BOND
        let line = self
            .read_v30_line()?
            .ok_or_else(|| SdfError::InvalidV3000Block("Expected BEGIN BOND".to_string()))?;

        if !line.starts_with("BEGIN BOND") {
            return Err(SdfError::InvalidV3000Block(format!(
                "Expected BEGIN BOND, got: {}",
                line
            )));
        }

        let mut bonds = Vec::with_capacity(count);

        for i in 0..count {
            let bond_line = self.read_v30_line()?.ok_or(SdfError::BondCountMismatch {
                expected: count,
                found: i,
            })?;

            let bond = self.parse_bond_line(&bond_line, id_to_index)?;
            bonds.push(bond);
        }

        // Read END BOND
        let end_line = self
            .read_v30_line()?
            .ok_or_else(|| SdfError::InvalidV3000Block("Expected END BOND".to_string()))?;

        if !end_line.starts_with("END BOND") {
            return Err(SdfError::InvalidV3000Block(format!(
                "Expected END BOND, got: {}",
                end_line
            )));
        }

        Ok(bonds)
    }

    /// Parses a single bond line.
    fn parse_bond_line(&self, line: &str, id_to_index: &HashMap<u32, usize>) -> Result<Bond> {
        // Format: index type atom1 atom2 [key=value...]
        // Example: 1 1 1 2 CFG=1
        let parts: Vec<&str> = line.split_whitespace().collect();

        if parts.len() < 4 {
            return Err(SdfError::InvalidV3000BondLine {
                line: self.line_number,
                message: format!("Not enough fields: {}", line),
            });
        }

        let v3000_id: u32 = parts[0]
            .parse()
            .map_err(|_| SdfError::InvalidV3000BondLine {
                line: self.line_number,
                message: format!("Invalid bond ID: {}", parts[0]),
            })?;

        let bond_type: u8 = parts[1]
            .parse()
            .map_err(|_| SdfError::InvalidV3000BondLine {
                line: self.line_number,
                message: format!("Invalid bond type: {}", parts[1]),
            })?;

        let atom1_id: u32 = parts[2]
            .parse()
            .map_err(|_| SdfError::InvalidV3000BondLine {
                line: self.line_number,
                message: format!("Invalid atom1 ID: {}", parts[2]),
            })?;

        let atom2_id: u32 = parts[3]
            .parse()
            .map_err(|_| SdfError::InvalidV3000BondLine {
                line: self.line_number,
                message: format!("Invalid atom2 ID: {}", parts[3]),
            })?;

        // Convert V3000 atom IDs to 0-based indices
        let atom1 = *id_to_index
            .get(&atom1_id)
            .ok_or(SdfError::AtomIdNotFound { id: atom1_id })?;

        let atom2 = *id_to_index
            .get(&atom2_id)
            .ok_or(SdfError::AtomIdNotFound { id: atom2_id })?;

        let order = BondOrder::from_sdf(bond_type).ok_or(SdfError::InvalidBondOrder(bond_type))?;

        // Parse key=value pairs
        let mut stereo = BondStereo::None;
        let mut topology: Option<u8> = None;
        let mut reacting_center: Option<u8> = None;

        for part in parts.iter().skip(4) {
            if let Some((key, value)) = part.split_once('=') {
                match key.to_uppercase().as_str() {
                    "CFG" => {
                        let cfg: u8 = value.parse().unwrap_or(0);
                        stereo = match cfg {
                            1 => BondStereo::Up,
                            3 => BondStereo::Down,
                            2 => BondStereo::Either,
                            _ => BondStereo::None,
                        };
                    }
                    "TOPO" => {
                        topology = value.parse().ok();
                    }
                    "RXCTR" => {
                        reacting_center = value.parse().ok();
                    }
                    _ => {} // Ignore unknown keys
                }
            }
        }

        Ok(Bond {
            atom1,
            atom2,
            order,
            stereo,
            topology,
            v3000_id: Some(v3000_id),
            reacting_center,
        })
    }

    /// Parses the SGROUP block.
    fn parse_sgroup_block(&mut self, id_to_index: &HashMap<u32, usize>) -> Result<Vec<SGroup>> {
        let mut sgroups = Vec::new();

        loop {
            let Some(line) = self.read_v30_line()? else {
                break;
            };

            if line.starts_with("END SGROUP") {
                break;
            }

            if let Some(sgroup) = self.parse_sgroup_line(&line, id_to_index)? {
                sgroups.push(sgroup);
            }
        }

        Ok(sgroups)
    }

    /// Parses a single SGROUP line.
    fn parse_sgroup_line(
        &self,
        line: &str,
        id_to_index: &HashMap<u32, usize>,
    ) -> Result<Option<SGroup>> {
        // Format: index type [ATOMS=(n a1 a2...)] [BONDS=(n b1 b2...)] [key=value...]
        let parts: Vec<&str> = line.split_whitespace().collect();

        if parts.len() < 2 {
            return Ok(None);
        }

        let id: u32 = match parts[0].parse() {
            Ok(v) => v,
            Err(_) => return Ok(None),
        };

        let sgroup_type = match SGroupType::parse(parts[1]) {
            Some(t) => t,
            None => return Ok(None),
        };

        let mut sgroup = SGroup::new(id, sgroup_type);

        // Parse the rest of the line for key=value pairs
        let rest = line[parts[0].len()..].trim();
        let rest = &rest[parts[1].len()..].trim();

        // Parse ATOMS=(n a1 a2...)
        if let Some(start) = rest.find("ATOMS=(") {
            if let Some(end) = rest[start..].find(')') {
                let inner = &rest[start + 7..start + end];
                let nums: Vec<&str> = inner.split_whitespace().collect();
                if !nums.is_empty() {
                    for num in nums.iter().skip(1) {
                        if let Ok(atom_id) = num.parse::<u32>() {
                            if let Some(&idx) = id_to_index.get(&atom_id) {
                                sgroup.atoms.push(idx);
                            }
                        }
                    }
                }
            }
        }

        // Parse LABEL="..."
        if let Some(start) = rest.find("LABEL=\"") {
            if let Some(end) = rest[start + 7..].find('"') {
                sgroup.label = Some(rest[start + 7..start + 7 + end].to_string());
            }
        }

        // Parse SUBTYPE for polymer connectivity
        if let Some(start) = rest.find("CONNECT=") {
            let value_start = start + 8;
            let value_end = rest[value_start..]
                .find(|c: char| c.is_whitespace())
                .map(|i| value_start + i)
                .unwrap_or(rest.len());
            sgroup.connectivity = Some(rest[value_start..value_end].to_string());
        }

        Ok(Some(sgroup))
    }

    /// Parses the COLLECTION block.
    fn parse_collection_block(
        &mut self,
        id_to_index: &HashMap<u32, usize>,
    ) -> Result<Vec<Collection>> {
        let mut collections = Vec::new();
        let mut current_id: u32 = 0;

        loop {
            let Some(line) = self.read_v30_line()? else {
                break;
            };

            if line.starts_with("END COLLECTION") {
                break;
            }

            // Parse collection type from line
            if let Some(coll) = self.parse_collection_line(&line, &mut current_id, id_to_index)? {
                collections.push(coll);
            }
        }

        Ok(collections)
    }

    /// Parses a single COLLECTION line.
    fn parse_collection_line(
        &self,
        line: &str,
        current_id: &mut u32,
        id_to_index: &HashMap<u32, usize>,
    ) -> Result<Option<Collection>> {
        // Collections have various formats
        // MDLV30/STEABS ATOMS=(n a1 a2...)
        // MDLV30/STEREL1 ATOMS=(n a1 a2...)

        if !line.contains('/') {
            return Ok(None);
        }

        *current_id += 1;

        let coll_type =
            if line.contains("STEABS") || line.contains("STEREL") || line.contains("STERAC") {
                CollectionType::AtomList // Use for stereo markers
            } else if line.contains("HILITE") {
                CollectionType::Highlight
            } else {
                return Ok(None);
            };

        let mut collection = Collection::new(coll_type, *current_id);

        // Parse ATOMS=(n a1 a2...)
        if let Some(start) = line.find("ATOMS=(") {
            if let Some(end) = line[start..].find(')') {
                let inner = &line[start + 7..start + end];
                let nums: Vec<&str> = inner.split_whitespace().collect();
                if !nums.is_empty() {
                    for num in nums.iter().skip(1) {
                        if let Ok(atom_id) = num.parse::<u32>() {
                            if let Some(&idx) = id_to_index.get(&atom_id) {
                                collection.atoms.push(idx);
                            }
                        }
                    }
                }
            }
        }

        Ok(Some(collection))
    }

    /// Parses the enhanced stereochemistry block.
    fn parse_stereo_block(
        &mut self,
        id_to_index: &HashMap<u32, usize>,
    ) -> Result<Vec<StereoGroup>> {
        let mut stereogroups = Vec::new();

        loop {
            let Some(line) = self.read_v30_line()? else {
                break;
            };

            if line.starts_with("END STEREO") {
                break;
            }

            if let Some(sg) = self.parse_stereo_line(&line, id_to_index)? {
                stereogroups.push(sg);
            }
        }

        Ok(stereogroups)
    }

    /// Parses a single stereo line.
    fn parse_stereo_line(
        &self,
        line: &str,
        id_to_index: &HashMap<u32, usize>,
    ) -> Result<Option<StereoGroup>> {
        // Format: type [n] ATOMS=(n a1 a2...)
        // Examples: ABS ATOMS=(2 1 2), OR1 ATOMS=(1 3), AND1 ATOMS=(1 4)

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            return Ok(None);
        }

        let type_str = parts[0];
        let (group_type, group_number) = if type_str == "ABS" {
            (StereoGroupType::Absolute, 0)
        } else if let Some(suffix) = type_str.strip_prefix("OR") {
            let num: u32 = suffix.parse().unwrap_or(1);
            (StereoGroupType::Or, num)
        } else if let Some(suffix) = type_str.strip_prefix("AND") {
            let num: u32 = suffix.parse().unwrap_or(1);
            (StereoGroupType::And, num)
        } else {
            return Ok(None);
        };

        let mut atoms = Vec::new();

        // Parse ATOMS=(n a1 a2...)
        if let Some(start) = line.find("ATOMS=(") {
            if let Some(end) = line[start..].find(')') {
                let inner = &line[start + 7..start + end];
                let nums: Vec<&str> = inner.split_whitespace().collect();
                if !nums.is_empty() {
                    for num in nums.iter().skip(1) {
                        if let Ok(atom_id) = num.parse::<u32>() {
                            if let Some(&idx) = id_to_index.get(&atom_id) {
                                atoms.push(idx);
                            }
                        }
                    }
                }
            }
        }

        Ok(Some(StereoGroup::new(group_type, group_number, atoms)))
    }

    /// Reads until M  END or $$$$.
    fn read_until_end(&mut self) -> Result<()> {
        loop {
            if !self.read_line()? {
                break;
            }
            if self.current_line.starts_with("M  END") || self.current_line.starts_with("$$$$") {
                break;
            }
        }
        Ok(())
    }

    /// Parses data block properties after M  END until $$$$.
    pub fn parse_properties(&mut self) -> Result<HashMap<String, String>> {
        let mut properties = HashMap::new();
        let mut current_property_name: Option<String> = None;
        let mut current_property_value = String::new();

        loop {
            if !self.read_line()? {
                break;
            }

            if self.current_line.starts_with("$$$$") {
                // Save any pending property
                if let Some(prop_name) = current_property_name.take() {
                    properties.insert(prop_name, current_property_value.trim().to_string());
                }
                break;
            }

            if self.current_line.starts_with("> ") || self.current_line.starts_with(">  ") {
                // Save previous property if exists
                if let Some(prop_name) = current_property_name.take() {
                    properties.insert(prop_name, current_property_value.trim().to_string());
                }
                current_property_value.clear();

                // Extract property name from angle brackets
                if let Some(start) = self.current_line.find('<') {
                    if let Some(end) = self.current_line[start + 1..].find('>') {
                        let prop_name = self.current_line[start + 1..start + 1 + end].to_string();
                        current_property_name = Some(prop_name);
                    }
                }
            } else if current_property_name.is_some() && !self.current_line.is_empty() {
                // Accumulate property value
                if !current_property_value.is_empty() {
                    current_property_value.push('\n');
                }
                current_property_value.push_str(&self.current_line);
            }
        }

        Ok(properties)
    }
}

/// Type alias for molecule header info: (name, program_line, comment)
type MolHeader = (String, Option<String>, Option<String>);

/// Iterator over molecules in a V3000 SDF file.
pub struct SdfV3000Iterator<R> {
    parser: SdfV3000Parser<R>,
    finished: bool,
}

impl<R: BufRead> SdfV3000Iterator<R> {
    /// Creates a new iterator from a buffered reader.
    pub fn new(reader: R) -> Self {
        Self {
            parser: SdfV3000Parser::new(reader),
            finished: false,
        }
    }

    /// Parses the header (first 4 lines) and returns molecule info.
    fn parse_header(&mut self) -> Result<Option<MolHeader>> {
        // Line 1: Molecule name
        if !self.parser.read_line()? {
            return Ok(None);
        }
        let name = self.parser.current_line.trim().to_string();

        // Line 2: Program/timestamp line
        if !self.parser.read_line()? {
            return Err(SdfError::MissingSection("header".to_string()));
        }
        let program_line = if self.parser.current_line.trim().is_empty() {
            None
        } else {
            Some(self.parser.current_line.clone())
        };

        // Line 3: Comment
        if !self.parser.read_line()? {
            return Err(SdfError::MissingSection("header".to_string()));
        }
        let comment = if self.parser.current_line.trim().is_empty() {
            None
        } else {
            Some(self.parser.current_line.clone())
        };

        // Line 4: Counts line (should indicate V3000)
        if !self.parser.read_line()? {
            return Err(SdfError::MissingSection("counts line".to_string()));
        }

        // V3000 counts line: "  0  0  0     0  0            999 V3000"
        if !self.parser.current_line.contains("V3000") {
            return Err(SdfError::InvalidV3000Block(
                "Expected V3000 in counts line".to_string(),
            ));
        }

        Ok(Some((name, program_line, comment)))
    }
}

impl<R: BufRead> Iterator for SdfV3000Iterator<R> {
    type Item = Result<Molecule>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        match self.parse_header() {
            Ok(Some((name, program_line, comment))) => {
                match self.parser.parse_molecule_body(name, program_line, comment) {
                    Ok(mut mol) => {
                        // Parse properties
                        match self.parser.parse_properties() {
                            Ok(props) => {
                                mol.properties = props;
                                Some(Ok(mol))
                            }
                            Err(e) => {
                                self.finished = true;
                                Some(Err(e))
                            }
                        }
                    }
                    Err(e) => {
                        self.finished = true;
                        Some(Err(e))
                    }
                }
            }
            Ok(None) => {
                self.finished = true;
                None
            }
            Err(e) => {
                self.finished = true;
                Some(Err(e))
            }
        }
    }
}

/// Parses a single V3000 molecule from a string.
pub fn parse_sdf_v3000_string(content: &str) -> Result<Molecule> {
    let cursor = std::io::Cursor::new(content);
    let reader = std::io::BufReader::new(cursor);
    let mut iter = SdfV3000Iterator::new(reader);

    iter.next().ok_or(SdfError::EmptyFile)?
}

/// Parses all V3000 molecules from a string.
pub fn parse_sdf_v3000_string_multi(content: &str) -> Result<Vec<Molecule>> {
    let cursor = std::io::Cursor::new(content);
    let reader = std::io::BufReader::new(cursor);
    let iter = SdfV3000Iterator::new(reader);

    iter.collect()
}

/// Parses a single V3000 molecule from a file.
pub fn parse_sdf_v3000_file<P: AsRef<std::path::Path>>(path: P) -> Result<Molecule> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let mut iter = SdfV3000Iterator::new(reader);

    iter.next().ok_or(SdfError::EmptyFile)?
}

/// Parses all V3000 molecules from a file.
pub fn parse_sdf_v3000_file_multi<P: AsRef<std::path::Path>>(path: P) -> Result<Vec<Molecule>> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let iter = SdfV3000Iterator::new(reader);

    iter.collect()
}

/// Returns an iterator over V3000 molecules in a file.
pub fn iter_sdf_v3000_file<P: AsRef<std::path::Path>>(
    path: P,
) -> Result<SdfV3000Iterator<std::io::BufReader<std::fs::File>>> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    Ok(SdfV3000Iterator::new(reader))
}

#[cfg(test)]
mod tests {
    use super::*;

    const V3000_METHANE: &str = r#"methane
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 2 H 0.6289 0.6289 0.6289 0
M  V30 3 H -0.6289 -0.6289 0.6289 0
M  V30 4 H -0.6289 0.6289 -0.6289 0
M  V30 5 H 0.6289 -0.6289 -0.6289 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;

    #[test]
    fn test_parse_v3000_methane() {
        let mol = parse_sdf_v3000_string(V3000_METHANE).unwrap();

        assert_eq!(mol.name, "methane");
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.bond_count(), 4);
        assert_eq!(mol.formula(), "CH4");
        assert_eq!(mol.format_version, SdfFormat::V3000);

        // Check carbon
        let carbon = &mol.atoms[0];
        assert_eq!(carbon.element, "C");
        assert_eq!(carbon.v3000_id, Some(1));

        // Check bonds
        for bond in &mol.bonds {
            assert_eq!(bond.order, BondOrder::Single);
        }
    }

    #[test]
    fn test_parse_v3000_with_charge() {
        let content = r#"charged
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 0.0000 0.0000 0.0000 0 CHG=1
M  V30 2 O 1.5000 0.0000 0.0000 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;
        let mol = parse_sdf_v3000_string(content).unwrap();

        assert_eq!(mol.atoms[0].formal_charge, 1);
        assert_eq!(mol.atoms[1].formal_charge, -1);
        assert_eq!(mol.total_charge(), 0);
    }

    #[test]
    fn test_parse_v3000_aromatic() {
        let content = r#"benzene
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 1.4000 0.0000 0
M  V30 2 C 1.2124 0.7000 0.0000 0
M  V30 3 C 1.2124 -0.7000 0.0000 0
M  V30 4 C 0.0000 -1.4000 0.0000 0
M  V30 5 C -1.2124 -0.7000 0.0000 0
M  V30 6 C -1.2124 0.7000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2
M  V30 2 4 2 3
M  V30 3 4 3 4
M  V30 4 4 4 5
M  V30 5 4 5 6
M  V30 6 4 6 1
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;
        let mol = parse_sdf_v3000_string(content).unwrap();

        assert_eq!(mol.name, "benzene");
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
        assert_eq!(mol.formula(), "C6");

        // All bonds should be aromatic
        for bond in &mol.bonds {
            assert_eq!(bond.order, BondOrder::Aromatic);
        }
        assert!(mol.has_aromatic_bonds());
    }
}
