use std::path::Path;
use std::fs::File;
use std::io;
use std::io::Read;
use rulinalg::matrix::Matrix;
use sxd_document::parser::parse;
use sxd_document::dom::Document;
use sxd_xpath::evaluate_xpath;

pub struct Scf {
    /// A matrix with columns giving the lattice vectors in Cartesian
    /// coordinates in units of Angstroms.
    pub d: Matrix<f64>,
    /// The Fermi energy in units of eV.
    pub fermi: f64,
}

/// Extract relevant data describing a DFT calculation and its results from
/// the file data-file.xml produced by Quantum Espresso.
///
/// # Arguments
/// * `scf_path` - path to the SCF data-file.xml file, produced by Quantum
/// Espresso as (prefix)/data-file.xml during the SCF calculation.
impl Scf {
    pub fn new<P: AsRef<Path>>(scf_path: P) -> Result<Scf, io::Error> {
        let mut contents = String::new();

        {
            let mut f = File::open(scf_path)?;
            f.read_to_string(&mut contents).expect(
                "error reading SCF data-file.xml file",
            );
        }

        let package = parse(&contents).unwrap();
        let doc = package.as_document();

        let angstrom_per_bohr = 0.52917721067;
        let d = extract_d_bohr(&doc) * angstrom_per_bohr;

        let ev_per_hartree = 27.21138602;
        let fermi = extract_fermi_hartree(&doc) * ev_per_hartree;

        Ok(Scf { d, fermi })
    }
}

fn extract_d_bohr(doc: &Document) -> Matrix<f64> {
    let base = "/Root/CELL/DIRECT_LATTICE_VECTORS";
    let units_path = format!("{}/UNITS_FOR_DIRECT_LATTICE_VECTORS/@UNITS", base);
    let units = evaluate_xpath(doc, &units_path).unwrap().into_string();
    assert_eq!(units, "Bohr");

    let mut d = Matrix::zeros(3, 3);

    for i in 0..3 {
        let a_path = format!("{}/a{}", base, i + 1);
        let a_text = evaluate_xpath(doc, &a_path).unwrap().into_string();
        let a_parts = a_text
            .trim()
            .split(" ")
            .filter(|x| x.len() > 0)
            .collect::<Vec<&str>>();

        for c in 0..3 {
            d[[c, i]] = a_parts[c].parse().unwrap();
        }
    }

    d
}

fn extract_fermi_hartree(doc: &Document) -> f64 {
    0.0
}
