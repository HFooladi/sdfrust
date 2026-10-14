//! Element data for molecular calculations.
//!
//! Provides atomic weights (IUPAC 2021 standard) and monoisotopic masses
//! for common elements found in chemical structures.

use std::sync::LazyLock;

use super::element_table::{ELEMENT_ROWS, ISOTOPE_MASSES};

/// Per-element data: weights, masses and valence rules.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElementData {
    /// Atomic number.
    pub atomic_number: u8,
    /// Element symbol.
    pub symbol: &'static str,
    /// Standard atomic weight.
    pub atomic_weight: f64,
    /// Monoisotopic mass (mass of the most abundant isotope).
    pub monoisotopic_mass: f64,
    /// Number of outer-shell (valence) electrons.
    pub outer_electrons: u8,
    /// Allowed valences in increasing order; empty means any valence
    /// (no implicit hydrogens are added, e.g. for metals).
    pub valences: &'static [u8],
}

/// All 118 elements plus the hydrogen isotope symbols D and T.
///
/// Data is generated from RDKit's periodic table (IUPAC/NIST values) by
/// `scripts/gen_element_table.py`, so masses agree with RDKit's descriptors.
static ELEMENTS: LazyLock<Vec<ElementData>> = LazyLock::new(|| {
    let mut table: Vec<ElementData> = ELEMENT_ROWS
        .iter()
        .map(
            |&(z, symbol, weight, mono, outer, valences, _)| ElementData {
                atomic_number: z,
                symbol,
                atomic_weight: weight,
                monoisotopic_mass: mono,
                outer_electrons: outer,
                valences,
            },
        )
        .collect();
    let hydrogen = table[0];
    for (symbol, mass_number) in [("D", 2), ("T", 3)] {
        let mass = isotope_mass(1, mass_number).unwrap_or(mass_number as f64);
        table.push(ElementData {
            symbol,
            atomic_weight: mass,
            monoisotopic_mass: mass,
            ..hydrogen
        });
    }
    table
});

/// Get element data by symbol (case-insensitive).
///
/// Returns `None` if the element is not found in the database.
///
/// # Example
///
/// ```rust
/// use sdfrust::descriptors::get_element;
///
/// let carbon = get_element("C").unwrap();
/// assert_eq!(carbon.atomic_number, 6);
/// assert!((carbon.atomic_weight - 12.011).abs() < 0.001);
/// assert_eq!(get_element("U").unwrap().atomic_number, 92);
/// ```
pub fn get_element(symbol: &str) -> Option<&'static ElementData> {
    let normalized = normalize_symbol(symbol);
    ELEMENTS.iter().find(|e| e.symbol == normalized)
}

/// Get element data by atomic number (1-118).
pub fn element_by_atomic_number(atomic_number: u8) -> Option<&'static ElementData> {
    if atomic_number as usize > ELEMENT_ROWS.len() {
        return None; // entries past the periodic table are the D/T aliases
    }
    ELEMENTS.get((atomic_number as usize).checked_sub(1)?)
}

/// Exact mass of a specific isotope, e.g. `isotope_mass(6, 13)` for carbon-13.
///
/// Covers all naturally occurring isotopes plus common radiolabels (³H, ¹¹C,
/// ¹⁴C, ¹⁸F, ¹²⁵I, ...). Returns `None` for isotopes not in the table.
pub fn isotope_mass(atomic_number: u8, mass_number: u16) -> Option<f64> {
    ISOTOPE_MASSES
        .binary_search_by(|&(z, a, _)| (z, a).cmp(&(atomic_number, mass_number)))
        .ok()
        .map(|i| ISOTOPE_MASSES[i].2)
}

/// Get atomic weight by element symbol.
///
/// Returns `None` if the element is not found.
///
/// # Example
///
/// ```rust
/// use sdfrust::descriptors::atomic_weight;
///
/// let weight = atomic_weight("O").unwrap();
/// assert!((weight - 15.999).abs() < 0.001);
/// ```
pub fn atomic_weight(symbol: &str) -> Option<f64> {
    get_element(symbol).map(|e| e.atomic_weight)
}

/// Get monoisotopic mass by element symbol.
///
/// Returns `None` if the element is not found.
///
/// # Example
///
/// ```rust
/// use sdfrust::descriptors::monoisotopic_mass;
///
/// let mass = monoisotopic_mass("C").unwrap();
/// assert!((mass - 12.0).abs() < 0.001);
/// ```
pub fn monoisotopic_mass(symbol: &str) -> Option<f64> {
    get_element(symbol).map(|e| e.monoisotopic_mass)
}

/// Covalent radius in Angstroms.
///
/// Source: Cordero et al. (2008), Dalton Trans., 2832-2838, for common
/// elements; other elements fall back to RDKit's covalent radii.
/// Returns `None` if the element is not found.
///
/// # Example
///
/// ```rust
/// use sdfrust::descriptors::covalent_radius;
///
/// let r = covalent_radius("C").unwrap();
/// assert!((r - 0.77).abs() < 0.01);
/// ```
pub fn covalent_radius(symbol: &str) -> Option<f64> {
    let normalized = normalize_symbol(symbol);
    COVALENT_RADII
        .iter()
        .find(|(s, _)| *s == normalized)
        .map(|(_, r)| *r)
        .or_else(|| {
            ELEMENT_ROWS
                .iter()
                .find(|row| row.1 == normalized)
                .map(|row| row.6)
        })
}

/// Covalent radii in Angstroms (Cordero et al. 2008).
static COVALENT_RADII: &[(&str, f64)] = &[
    ("H", 0.31),
    ("He", 0.28),
    ("Li", 1.28),
    ("Be", 0.96),
    ("B", 0.84),
    ("C", 0.77),
    ("N", 0.71),
    ("O", 0.66),
    ("F", 0.57),
    ("Ne", 0.58),
    ("Na", 1.66),
    ("Mg", 1.41),
    ("Al", 1.21),
    ("Si", 1.11),
    ("P", 1.07),
    ("S", 1.05),
    ("Cl", 1.02),
    ("Ar", 1.06),
    ("K", 2.03),
    ("Ca", 1.76),
    ("Mn", 1.39),
    ("Fe", 1.32),
    ("Co", 1.26),
    ("Ni", 1.24),
    ("Cu", 1.32),
    ("Zn", 1.22),
    ("As", 1.19),
    ("Se", 1.20),
    ("Br", 1.20),
    ("I", 1.39),
    ("Pt", 1.36),
    ("Au", 1.36),
    ("D", 0.31),
    ("T", 0.31),
];

/// Normalize element symbol (capitalize first letter, lowercase rest).
fn normalize_symbol(symbol: &str) -> String {
    let symbol = symbol.trim();
    let mut chars = symbol.chars();
    match chars.next() {
        None => String::new(),
        Some(first) => first.to_uppercase().collect::<String>() + &chars.as_str().to_lowercase(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_element_carbon() {
        let c = get_element("C").unwrap();
        assert_eq!(c.atomic_number, 6);
        assert_eq!(c.symbol, "C");
        assert!((c.atomic_weight - 12.011).abs() < 0.001);
        assert!((c.monoisotopic_mass - 12.0).abs() < 0.001);
    }

    #[test]
    fn test_get_element_hydrogen() {
        let h = get_element("H").unwrap();
        assert_eq!(h.atomic_number, 1);
        assert!((h.atomic_weight - 1.008).abs() < 0.001);
    }

    #[test]
    fn test_get_element_case_insensitive() {
        assert!(get_element("c").is_some());
        assert!(get_element("C").is_some());
        assert!(get_element("cL").is_some());
        assert!(get_element("Cl").is_some());
        assert!(get_element("CL").is_some());
    }

    #[test]
    fn test_get_element_unknown() {
        assert!(get_element("Xx").is_none());
        assert!(get_element("").is_none());
    }

    #[test]
    fn test_atomic_weight() {
        assert!((atomic_weight("O").unwrap() - 15.999).abs() < 0.001);
        assert!((atomic_weight("N").unwrap() - 14.007).abs() < 0.001);
    }

    #[test]
    fn test_monoisotopic_mass() {
        // Carbon-12 is exactly 12.0 by definition
        assert!((monoisotopic_mass("C").unwrap() - 12.0).abs() < 0.0001);
    }

    #[test]
    fn test_deuterium_tritium() {
        let d = get_element("D").unwrap();
        assert_eq!(d.atomic_number, 1);
        assert!(d.atomic_weight > 2.0);

        let t = get_element("T").unwrap();
        assert_eq!(t.atomic_number, 1);
        assert!(t.atomic_weight > 3.0);
    }

    #[test]
    fn test_common_organic_elements() {
        // Verify all common organic elements are present
        for symbol in &["H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"] {
            assert!(
                get_element(symbol).is_some(),
                "Element {} should be in database",
                symbol
            );
        }
    }

    #[test]
    fn test_common_metals() {
        // Verify common metals are present
        for symbol in &["Na", "K", "Ca", "Mg", "Fe", "Cu", "Zn"] {
            assert!(
                get_element(symbol).is_some(),
                "Element {} should be in database",
                symbol
            );
        }
    }

    #[test]
    fn test_covalent_radius_carbon() {
        let r = covalent_radius("C").unwrap();
        assert!((r - 0.77).abs() < 0.01);
    }

    #[test]
    fn test_covalent_radius_hydrogen() {
        let r = covalent_radius("H").unwrap();
        assert!((r - 0.31).abs() < 0.01);
    }

    #[test]
    fn test_covalent_radius_case_insensitive() {
        assert!(covalent_radius("c").is_some());
        assert!(covalent_radius("CL").is_some());
        assert!(covalent_radius("br").is_some());
    }

    #[test]
    fn test_full_periodic_table() {
        for z in 1..=118u8 {
            let e = element_by_atomic_number(z).unwrap();
            assert_eq!(e.atomic_number, z);
            assert_eq!(get_element(e.symbol).unwrap().atomic_number, z);
            assert!(covalent_radius(e.symbol).is_some());
        }
        assert_eq!(get_element("U").unwrap().atomic_number, 92);
        assert!(element_by_atomic_number(0).is_none());
        assert!(element_by_atomic_number(119).is_none());
    }

    #[test]
    fn test_isotope_mass() {
        assert!((isotope_mass(6, 13).unwrap() - 13.003355).abs() < 1e-5);
        assert!((isotope_mass(1, 2).unwrap() - 2.014102).abs() < 1e-5);
        assert!((isotope_mass(9, 18).unwrap() - 18.000938).abs() < 1e-5);
        assert!(isotope_mass(6, 99).is_none());
    }

    #[test]
    fn test_covalent_radius_unknown() {
        assert!(covalent_radius("Xx").is_none());
        assert!(covalent_radius("").is_none());
    }
}
