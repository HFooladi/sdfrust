//! Element data for molecular calculations.
//!
//! Provides atomic weights (IUPAC 2021 standard) and monoisotopic masses
//! for common elements found in chemical structures.

/// Element data with atomic weights and monoisotopic masses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElementData {
    /// Atomic number.
    pub atomic_number: u8,
    /// Element symbol.
    pub symbol: &'static str,
    /// Standard atomic weight (IUPAC 2021).
    pub atomic_weight: f64,
    /// Monoisotopic mass (mass of most abundant isotope).
    pub monoisotopic_mass: f64,
}

/// Element data table for common elements.
/// Data sources:
/// - Atomic weights: IUPAC 2021 (https://www.ciaaw.org/atomic-weights.htm)
/// - Monoisotopic masses: NIST (https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl)
static ELEMENTS: &[ElementData] = &[
    ElementData {
        atomic_number: 1,
        symbol: "H",
        atomic_weight: 1.008,
        monoisotopic_mass: 1.00782503207,
    },
    ElementData {
        atomic_number: 2,
        symbol: "He",
        atomic_weight: 4.002602,
        monoisotopic_mass: 4.00260325413,
    },
    ElementData {
        atomic_number: 3,
        symbol: "Li",
        atomic_weight: 6.94,
        monoisotopic_mass: 7.01600343426,
    },
    ElementData {
        atomic_number: 4,
        symbol: "Be",
        atomic_weight: 9.0121831,
        monoisotopic_mass: 9.01218306500,
    },
    ElementData {
        atomic_number: 5,
        symbol: "B",
        atomic_weight: 10.81,
        monoisotopic_mass: 11.00930536000,
    },
    ElementData {
        atomic_number: 6,
        symbol: "C",
        atomic_weight: 12.011,
        monoisotopic_mass: 12.00000000000,
    },
    ElementData {
        atomic_number: 7,
        symbol: "N",
        atomic_weight: 14.007,
        monoisotopic_mass: 14.00307400443,
    },
    ElementData {
        atomic_number: 8,
        symbol: "O",
        atomic_weight: 15.999,
        monoisotopic_mass: 15.99491461960,
    },
    ElementData {
        atomic_number: 9,
        symbol: "F",
        atomic_weight: 18.998403163,
        monoisotopic_mass: 18.99840316273,
    },
    ElementData {
        atomic_number: 10,
        symbol: "Ne",
        atomic_weight: 20.1797,
        monoisotopic_mass: 19.99244017617,
    },
    ElementData {
        atomic_number: 11,
        symbol: "Na",
        atomic_weight: 22.98976928,
        monoisotopic_mass: 22.98976928200,
    },
    ElementData {
        atomic_number: 12,
        symbol: "Mg",
        atomic_weight: 24.305,
        monoisotopic_mass: 23.98504169700,
    },
    ElementData {
        atomic_number: 13,
        symbol: "Al",
        atomic_weight: 26.9815384,
        monoisotopic_mass: 26.98153853000,
    },
    ElementData {
        atomic_number: 14,
        symbol: "Si",
        atomic_weight: 28.085,
        monoisotopic_mass: 27.97692653465,
    },
    ElementData {
        atomic_number: 15,
        symbol: "P",
        atomic_weight: 30.973761998,
        monoisotopic_mass: 30.97376199800,
    },
    ElementData {
        atomic_number: 16,
        symbol: "S",
        atomic_weight: 32.06,
        monoisotopic_mass: 31.97207117400,
    },
    ElementData {
        atomic_number: 17,
        symbol: "Cl",
        atomic_weight: 35.45,
        monoisotopic_mass: 34.96885268200,
    },
    ElementData {
        atomic_number: 18,
        symbol: "Ar",
        atomic_weight: 39.948,
        monoisotopic_mass: 39.96238312463,
    },
    ElementData {
        atomic_number: 19,
        symbol: "K",
        atomic_weight: 39.0983,
        monoisotopic_mass: 38.96370668200,
    },
    ElementData {
        atomic_number: 20,
        symbol: "Ca",
        atomic_weight: 40.078,
        monoisotopic_mass: 39.96259098000,
    },
    ElementData {
        atomic_number: 25,
        symbol: "Mn",
        atomic_weight: 54.938043,
        monoisotopic_mass: 54.93804391000,
    },
    ElementData {
        atomic_number: 26,
        symbol: "Fe",
        atomic_weight: 55.845,
        monoisotopic_mass: 55.93493633000,
    },
    ElementData {
        atomic_number: 27,
        symbol: "Co",
        atomic_weight: 58.933194,
        monoisotopic_mass: 58.93319429000,
    },
    ElementData {
        atomic_number: 28,
        symbol: "Ni",
        atomic_weight: 58.6934,
        monoisotopic_mass: 57.93534241000,
    },
    ElementData {
        atomic_number: 29,
        symbol: "Cu",
        atomic_weight: 63.546,
        monoisotopic_mass: 62.92959772000,
    },
    ElementData {
        atomic_number: 30,
        symbol: "Zn",
        atomic_weight: 65.38,
        monoisotopic_mass: 63.92914201000,
    },
    ElementData {
        atomic_number: 33,
        symbol: "As",
        atomic_weight: 74.921595,
        monoisotopic_mass: 74.92159457000,
    },
    ElementData {
        atomic_number: 34,
        symbol: "Se",
        atomic_weight: 78.971,
        monoisotopic_mass: 79.91652176000,
    },
    ElementData {
        atomic_number: 35,
        symbol: "Br",
        atomic_weight: 79.904,
        monoisotopic_mass: 78.91833710000,
    },
    ElementData {
        atomic_number: 53,
        symbol: "I",
        atomic_weight: 126.90447,
        monoisotopic_mass: 126.90447190000,
    },
    ElementData {
        atomic_number: 78,
        symbol: "Pt",
        atomic_weight: 195.084,
        monoisotopic_mass: 194.96479120000,
    },
    ElementData {
        atomic_number: 79,
        symbol: "Au",
        atomic_weight: 196.966570,
        monoisotopic_mass: 196.96656879000,
    },
    // Deuterium (special handling as D)
    ElementData {
        atomic_number: 1,
        symbol: "D",
        atomic_weight: 2.014101778,
        monoisotopic_mass: 2.01410177812,
    },
    // Tritium (special handling as T)
    ElementData {
        atomic_number: 1,
        symbol: "T",
        atomic_weight: 3.01604928,
        monoisotopic_mass: 3.01604928199,
    },
];

/// Get element data by symbol.
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
/// ```
pub fn get_element(symbol: &str) -> Option<&'static ElementData> {
    // Normalize symbol: uppercase first letter, lowercase rest
    let normalized = normalize_symbol(symbol);
    ELEMENTS.iter().find(|e| e.symbol == normalized)
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

/// Normalize element symbol (capitalize first letter, lowercase rest).
fn normalize_symbol(symbol: &str) -> String {
    let symbol = symbol.trim();
    let mut chars = symbol.chars();
    match chars.next() {
        None => String::new(),
        Some(first) => {
            first.to_uppercase().collect::<String>() + &chars.as_str().to_lowercase()
        }
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
}
