use std::path::Path;
use std::fs::File;
use std::io;
use std::io::Read;
use ndarray::Array2;
use sxd_document::parser::parse;
use sxd_document::dom::Document;
use sxd_xpath::evaluate_xpath;

use units::{ANGSTROM_PER_BOHR, EV_PER_HARTREE};

pub struct Scf {
    /// A matrix with columns giving the lattice vectors in Cartesian
    /// coordinates in units of Angstroms.
    pub d: Array2<f64>,
    /// The Fermi energy in units of eV.
    pub fermi: f64,
    /// The lattice parameter (in Bohr) used by Quantum Espresso for units of some
    /// quantities.
    pub alat: f64,
}

impl Scf {
    /// Extract relevant data describing a DFT calculation and its results from
    /// the file `data-file.xml` produced by Quantum Espresso.
    ///
    /// # Arguments
    /// * `scf_path` - path to the SCF data-file.xml file, produced by Quantum
    /// Espresso as (prefix)/data-file.xml during the SCF calculation.
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

        let d = extract_d_bohr(&doc) * ANGSTROM_PER_BOHR;

        let fermi = extract_fermi_hartree(&doc) * EV_PER_HARTREE;

        let alat = extract_alat_bohr(&doc);

        Ok(Scf { d, fermi, alat })
    }
}

fn extract_d_bohr(doc: &Document) -> Array2<f64> {
    let base = "/Root/CELL/DIRECT_LATTICE_VECTORS";
    let units_path = format!("{}/UNITS_FOR_DIRECT_LATTICE_VECTORS/@UNITS", base);
    let units = evaluate_xpath(doc, &units_path).unwrap().into_string();
    assert_eq!(units, "Bohr");

    let mut d = Array2::zeros((3, 3));

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
    let base = "/Root/BAND_STRUCTURE_INFO";
    let units_path = format!("{}/UNITS_FOR_ENERGIES/@UNITS", base);
    let units = evaluate_xpath(doc, &units_path).unwrap().into_string();
    assert_eq!(units, "Hartree");

    let fermi_path = format!("{}/FERMI_ENERGY", base);
    let fermi_text = evaluate_xpath(doc, &fermi_path).unwrap().into_string();
    let fermi = fermi_text.trim().parse().unwrap();

    fermi
}

fn extract_alat_bohr(doc: &Document) -> f64 {
    let alat_path = "/Root/CELL/LATTICE_PARAMETER";
    let units_path = format!("{}/@UNITS", alat_path);
    let units = evaluate_xpath(doc, &units_path).unwrap().into_string();
    assert_eq!(units, "Bohr");

    let alat_text = evaluate_xpath(doc, alat_path).unwrap().into_string();
    let alat = alat_text.trim().parse().unwrap();

    alat
}
