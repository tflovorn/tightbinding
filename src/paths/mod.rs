use std::env;
use std::path::PathBuf;

mod errors {
    error_chain!{}
}

use self::errors::*;

/// Return the directory under which all DFT/Wannier90 data to be considered
/// is expected to be found. This directory is given by the environment
/// variable `TB_W90_WORK`.
pub fn get_work_base() -> Result<String> {
    let work_err_msg = "Could not find environment variable TB_W90_WORK,\
                        which should be the data root directory.";
    let work_base = env::var("TB_W90_WORK").chain_err(|| work_err_msg);

    work_base
}

/// Construct the path to the work directory for an individual calculation.
/// This path has the form `{work_base}/{subdir}/{prefix}`.
/// Subdir may be omitted.
pub fn build_work(work_base: &str, subdir: Option<&str>, prefix: &str) -> PathBuf {
    let mut work = PathBuf::new();
    work.push(work_base);

    if let Some(ref subdir) = subdir {
        work.push(subdir);
    }

    work.push(prefix);

    work
}

/// Construct the path to the SCF calculation output `data-file.xml`.
/// This is assumed to be found in `{work}/scf/{prefix}.save/data-file.xml`.
pub fn get_scf_path(work: &PathBuf, prefix: &str) -> PathBuf {
    let mut scf_path = work.clone();
    scf_path.push("scf");
    scf_path.push(format!("{}.save", prefix));
    scf_path.push("data-file.xml");

    scf_path
}

/// Construct the path to the Wannier90 tight-binding model
/// file `{prefix}_hr.dat`.
/// This is assumed to be found in `{work}/wannier/{prefix}_hr.dat`.
pub fn get_hr_path(work: &PathBuf, prefix: &str) -> PathBuf {
    let mut hr_path = work.clone();
    hr_path.push("wannier");
    hr_path.push(format!("{}_hr.dat", prefix));

    hr_path
}
