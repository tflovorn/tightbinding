extern crate serde_json;
extern crate tightbinding;

use std::path::PathBuf;
use std::fs::File;
use tightbinding::float::is_near_float;
use tightbinding::dos::DosValues;

pub fn write_dos_out(dos: &DosValues, test_prefix: &str) {
    let out_path = PathBuf::from(format!("plot/test_{}_dos.json", test_prefix));
    let out_file = File::create(out_path).expect(&format!(
        "Error creating output file: {}_dos expects to be run in crate root.",
        test_prefix
    ));
    serde_json::to_writer(&out_file, &dos).expect("Error writing output file");
}

pub fn read_dos(test_prefix: &str) -> DosValues {
    let expected_path = PathBuf::from(format!(
        "test_data/{}/test_{}_dos.json",
        test_prefix,
        test_prefix
    ));
    let expected_file = File::open(expected_path).expect(&format!(
        "Error opening expected value file: {}_dos expects to be run in crate root.",
        test_prefix
    ));
    let expected_dos: DosValues =
        serde_json::from_reader(&expected_file).expect("Error reading expected data file");

    expected_dos
}

pub fn check_dos(
    dos: &DosValues,
    expected_dos: &DosValues,
    eps_abs_e: f64,
    eps_rel_e: f64,
    eps_abs_dos: f64,
    eps_rel_dos: f64,
) {
    for (e_out, e_expected) in dos.es.iter().zip(expected_dos.es.iter()) {
        assert!(is_near_float(*e_out, *e_expected, eps_abs_e, eps_rel_e));
    }

    for (band_dos_out, band_dos_expected) in
        dos.orbital_dos.iter().zip(expected_dos.orbital_dos.iter())
    {
        for (dos_out, dos_expected) in band_dos_out.iter().zip(band_dos_expected.iter()) {
            assert!(is_near_float(
                *dos_out,
                *dos_expected,
                eps_abs_dos,
                eps_rel_dos,
            ));
        }
    }

    for (dos_out, dos_expected) in dos.total_dos.iter().zip(expected_dos.total_dos.iter()) {
        assert!(is_near_float(
            *dos_out,
            *dos_expected,
            eps_abs_dos,
            eps_rel_dos,
        ));
    }
}
