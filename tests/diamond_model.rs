extern crate num_complex;
extern crate rulinalg;
extern crate tightbinding;

use num_complex::Complex64;
use rulinalg::matrix::Matrix;
use tightbinding::Model;
use tightbinding::w90::W90Model;
use tightbinding::qe::Scf;

mod common;
use common::is_near;

#[test]
fn diamond_model() {
    let scf_path = "test_data/diamond/scf.data-file.xml";
    let hr_path = "test_data/diamond/diamond_hr.dat";

    let scf_data = Scf::new(scf_path).unwrap();
    let d = scf_data.d;

    let expected_d_data = vec![
        -1.6139904925435,
        0.0,
        -1.6139904925435,
        0.0,
        1.6139904925435,
        1.6139904925435,
        1.6139904925435,
        1.6139904925435,
        0.0,
    ];
    check_d(&d, expected_d_data);

    let m = W90Model::new(hr_path, d).unwrap();

    assert_eq!(m.bands(), 4);
    assert_eq!(m.hrs().len(), 93);

    let expected_rs = vec![
        [-3, 1, 1],
        [-3, 1, 1],
        [3, -1, -1],
        [3, -1, -1],
        [-2, -2, 2],
        [-2, -2, 2],
    ];
    let expected_indexes = vec![[0, 0], [2, 1], [0, 0], [2, 1], [0, 0], [2, 3]];
    let expected_vals = vec![
        Complex64::new(0.007378 / 4.0, 0.0),
        Complex64::new(-0.008540 / 4.0, 0.0),
        Complex64::new(0.007378 / 4.0, 0.0),
        Complex64::new(-0.008540 / 4.0, 0.0),
        Complex64::new(0.011647 / 6.0, 0.0),
        Complex64::new(-0.003502 / 6.0, 0.0),
    ];

    check_hrs(&m, expected_rs, expected_indexes, expected_vals);
}

fn check_d(d: &Matrix<f64>, expected_d_data: Vec<f64>) {
    let eps_abs = 1e-12; // Angstrom
    let eps_rel = 1e-12;

    for (x, y) in expected_d_data.iter().zip(d.data()) {
        let diff = (x - y).abs();

        assert!(diff < eps_abs || diff < eps_rel * x.abs().max(y.abs()));
    }
}

fn check_hrs(
    m: &W90Model,
    expected_rs: Vec<[i32; 3]>,
    expected_indexes: Vec<[usize; 2]>,
    expected_vals: Vec<Complex64>,
) {
    let eps_abs = 1e-12; // eV
    let eps_rel = 1e-12;

    for ((r, indexes), val) in expected_rs.iter().zip(expected_indexes).zip(expected_vals) {
        assert!(is_near(m.hrs()[r][indexes], val, eps_abs, eps_rel));
    }
}