pub mod mol2;
pub mod sdf;
pub mod sdf_v3000;
pub mod xyz;

#[cfg(feature = "gzip")]
pub mod compression;
#[cfg(feature = "gzip")]
pub use compression::{MaybeGzReader, is_gzip_path, open_maybe_gz, read_maybe_gz_to_string};

pub use sdf::{
    AutoIterator, FileFormat, SdfIterator, SdfParser, detect_format, detect_sdf_format,
    iter_auto_file, iter_sdf_file, parse_auto_file, parse_auto_file_multi, parse_auto_string,
    parse_auto_string_multi, parse_sdf_auto_file, parse_sdf_auto_file_multi, parse_sdf_auto_string,
    parse_sdf_auto_string_multi, parse_sdf_file, parse_sdf_file_multi, parse_sdf_string,
    parse_sdf_string_multi,
};

pub use mol2::{
    Mol2Iterator, Mol2Parser, iter_mol2_file, parse_mol2_file, parse_mol2_file_multi,
    parse_mol2_string, parse_mol2_string_multi,
};

pub use sdf_v3000::{
    SdfV3000Iterator, SdfV3000Parser, iter_sdf_v3000_file, parse_sdf_v3000_file,
    parse_sdf_v3000_file_multi, parse_sdf_v3000_string, parse_sdf_v3000_string_multi,
};

pub use xyz::{
    XyzIterator, XyzParser, iter_xyz_file, parse_xyz_file, parse_xyz_file_multi, parse_xyz_string,
    parse_xyz_string_multi,
};
