use sdfrust::{SdfError, parse_sdf_file};
use std::collections::HashMap;
use std::path::PathBuf;
use std::time::Instant;

const YEAR_RANGES: &[&str] = &["1981-2000", "2001-2010", "2011-2020", "2021-2023", "demo"];

fn categorize_error(err: &SdfError) -> &'static str {
    match err {
        SdfError::Io(_) => "Io",
        SdfError::Parse { .. } => "Parse",
        SdfError::AtomCountMismatch { .. } => "AtomCountMismatch",
        SdfError::BondCountMismatch { .. } => "BondCountMismatch",
        SdfError::InvalidAtomIndex { .. } => "InvalidAtomIndex",
        SdfError::InvalidBondOrder(_) => "InvalidBondOrder",
        SdfError::InvalidCountsLine(_) => "InvalidCountsLine",
        SdfError::MissingSection(_) => "MissingSection",
        SdfError::EmptyFile => "EmptyFile",
        SdfError::InvalidCoordinate(_) => "InvalidCoordinate",
        SdfError::InvalidCharge(_) => "InvalidCharge",
        SdfError::InvalidV3000Block(_) => "InvalidV3000Block",
        SdfError::InvalidV3000AtomLine { .. } => "InvalidV3000AtomLine",
        SdfError::InvalidV3000BondLine { .. } => "InvalidV3000BondLine",
        SdfError::AtomIdNotFound { .. } => "AtomIdNotFound",
        SdfError::UnsupportedV3000Feature(_) => "UnsupportedV3000Feature",
        SdfError::GzipNotEnabled => "GzipNotEnabled",
        SdfError::BondInferenceError { .. } => "BondInferenceError",
    }
}

/// Large-scale SDF parsing benchmark against PDBbind 2024 dataset.
///
/// This test walks ~27,670 ligand SDF files from the PDBbind 2024 dataset
/// and measures parsing success/failure rates, error categories, and
/// basic molecule statistics.
///
/// Run with:
/// ```bash
/// cargo test --release pdbbind_benchmark -- --ignored --nocapture
/// ```
#[test]
#[ignore]
fn pdbbind_benchmark() {
    let dataset_dir = match std::env::var("PDBBIND_2024_DIR") {
        Ok(dir) => dir,
        Err(_) => {
            println!("SKIP: Set PDBBIND_2024_DIR to run this benchmark");
            return;
        }
    };
    let dataset_path = PathBuf::from(&dataset_dir);
    if !dataset_path.exists() {
        println!("SKIP: PDBbind 2024 dataset not found at {}", dataset_dir);
        return;
    }

    // Discover all *_ligand.sdf files
    let mut files_by_range: HashMap<String, Vec<PathBuf>> = HashMap::new();
    for &year_range in YEAR_RANGES {
        let range_dir = dataset_path.join(year_range);
        if !range_dir.is_dir() {
            continue;
        }
        let entries = match std::fs::read_dir(&range_dir) {
            Ok(e) => e,
            Err(_) => continue,
        };
        let mut paths = Vec::new();
        for entry in entries.flatten() {
            if !entry.path().is_dir() {
                continue;
            }
            let pdb_id = entry.file_name();
            let sdf_file = entry
                .path()
                .join(format!("{}_ligand.sdf", pdb_id.to_string_lossy()));
            if sdf_file.exists() {
                paths.push(sdf_file);
            }
        }
        paths.sort();
        files_by_range.insert(year_range.to_string(), paths);
    }

    let total_files: usize = files_by_range.values().map(|v| v.len()).sum();
    println!();
    println!("{}", "=".repeat(70));
    println!("PDBbind 2024 SDF Parsing Benchmark");
    println!("{}", "=".repeat(70));
    println!("Dataset: {}", dataset_dir);
    println!("Total ligand SDF files discovered: {}", total_files);
    println!();

    if total_files == 0 {
        println!("SKIP: No SDF files found in the dataset directory.");
        return;
    }

    // Parse all files and collect statistics
    let mut total_success = 0usize;
    let mut total_failure = 0usize;
    let mut error_counts: HashMap<&str, usize> = HashMap::new();
    let mut failed_files: Vec<(String, String)> = Vec::new(); // (path, error)
    let mut range_stats: HashMap<String, (usize, usize)> = HashMap::new(); // (success, failure)

    // Molecule stats
    let mut min_atoms = usize::MAX;
    let mut max_atoms = 0usize;
    let mut total_atoms = 0u64;
    let mut min_bonds = usize::MAX;
    let mut max_bonds = 0usize;
    let mut total_bonds = 0u64;
    let mut element_freq: HashMap<String, u64> = HashMap::new();

    let start = Instant::now();

    for (year_range, files) in &files_by_range {
        let mut range_success = 0usize;
        let mut range_failure = 0usize;

        for file_path in files {
            match parse_sdf_file(file_path) {
                Ok(mol) => {
                    total_success += 1;
                    range_success += 1;

                    let ac = mol.atom_count();
                    let bc = mol.bond_count();
                    min_atoms = min_atoms.min(ac);
                    max_atoms = max_atoms.max(ac);
                    total_atoms += ac as u64;
                    min_bonds = min_bonds.min(bc);
                    max_bonds = max_bonds.max(bc);
                    total_bonds += bc as u64;

                    for atom in mol.atoms.iter() {
                        *element_freq.entry(atom.element.clone()).or_insert(0) += 1;
                    }
                }
                Err(e) => {
                    total_failure += 1;
                    range_failure += 1;
                    let category = categorize_error(&e);
                    *error_counts.entry(category).or_insert(0) += 1;
                    if failed_files.len() < 20 {
                        failed_files
                            .push((file_path.to_string_lossy().to_string(), format!("{}", e)));
                    }
                }
            }
        }

        range_stats.insert(year_range.clone(), (range_success, range_failure));
    }

    let elapsed = start.elapsed();
    let total = total_success + total_failure;

    // Print results
    println!("{}", "-".repeat(70));
    println!("OVERALL RESULTS");
    println!("{}", "-".repeat(70));
    println!(
        "Total parsed:  {} / {} ({:.2}% success)",
        total_success,
        total,
        if total > 0 {
            total_success as f64 / total as f64 * 100.0
        } else {
            0.0
        }
    );
    println!(
        "Total failed:  {} / {} ({:.2}% failure)",
        total_failure,
        total,
        if total > 0 {
            total_failure as f64 / total as f64 * 100.0
        } else {
            0.0
        }
    );
    println!(
        "Time:          {:.2}s ({:.0} files/sec)",
        elapsed.as_secs_f64(),
        if elapsed.as_secs_f64() > 0.0 {
            total as f64 / elapsed.as_secs_f64()
        } else {
            0.0
        }
    );

    // Per year-range breakdown
    println!();
    println!("{}", "-".repeat(70));
    println!("PER YEAR-RANGE BREAKDOWN");
    println!("{}", "-".repeat(70));
    println!(
        "{:<15} {:>8} {:>8} {:>8} {:>10}",
        "Range", "Total", "Success", "Failed", "Rate"
    );
    let mut sorted_ranges: Vec<_> = range_stats.iter().collect();
    sorted_ranges.sort_by_key(|(k, _)| (*k).clone());
    for (range, (success, failure)) in &sorted_ranges {
        let range_total = success + failure;
        println!(
            "{:<15} {:>8} {:>8} {:>8} {:>9.2}%",
            range,
            range_total,
            success,
            failure,
            if range_total > 0 {
                *success as f64 / range_total as f64 * 100.0
            } else {
                0.0
            }
        );
    }

    // Molecule statistics
    if total_success > 0 {
        println!();
        println!("{}", "-".repeat(70));
        println!("MOLECULE STATISTICS (successfully parsed)");
        println!("{}", "-".repeat(70));
        let avg_atoms = total_atoms as f64 / total_success as f64;
        let avg_bonds = total_bonds as f64 / total_success as f64;
        println!(
            "Atoms:  min={}, max={}, avg={:.1}",
            min_atoms, max_atoms, avg_atoms
        );
        println!(
            "Bonds:  min={}, max={}, avg={:.1}",
            min_bonds, max_bonds, avg_bonds
        );

        // Top 20 elements
        let mut elements: Vec<_> = element_freq.iter().collect();
        elements.sort_by(|a, b| b.1.cmp(a.1));
        println!();
        println!("Top elements:");
        for (elem, count) in elements.iter().take(20) {
            println!("  {:<4} {:>10}", elem, count);
        }
    }

    // Error breakdown
    if !error_counts.is_empty() {
        println!();
        println!("{}", "-".repeat(70));
        println!("ERROR BREAKDOWN");
        println!("{}", "-".repeat(70));
        let mut errors: Vec<_> = error_counts.iter().collect();
        errors.sort_by(|a, b| b.1.cmp(a.1));
        for (category, count) in &errors {
            println!("  {:<30} {:>6}", category, count);
        }
    }

    // Sample failed files
    if !failed_files.is_empty() {
        println!();
        println!("{}", "-".repeat(70));
        println!("SAMPLE FAILED FILES (first {})", failed_files.len());
        println!("{}", "-".repeat(70));
        for (path, err) in &failed_files {
            println!("  {} -> {}", path, err);
        }
    }

    println!();
    println!("{}", "=".repeat(70));
    println!("Benchmark complete.");
}
