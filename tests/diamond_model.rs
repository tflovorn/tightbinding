extern crate num_complex;
extern crate rulinalg;
extern crate tightbinding;

use num_complex::Complex64;
use rulinalg::matrix::Matrix;
use tightbinding::Model;
use tightbinding::w90::W90Model;

mod common;
use common::is_near;

#[test]
fn diamond_model() {
    let hr_path = "test_data/diamond/diamond_hr.dat";

    // TODO extract d from scf.out
    let d = Matrix::identity(3);

    let m = W90Model::new(hr_path, d).unwrap();

    assert_eq!(m.bands(), 4);
    assert_eq!(m.hrs().len(), 93);

    let eps_abs = 1e-12; // eV
    let eps_rel = 1e-12;

    assert!(is_near(
        m.hrs()[&[-3, 1, 1]][[0, 0]],
        Complex64::new(0.007378 / 4.0, 0.0),
        eps_abs,
        eps_rel,
    ));
    assert!(is_near(
        m.hrs()[&[-3, 1, 1]][[2, 1]],
        Complex64::new(-0.008540 / 4.0, 0.0),
        eps_abs,
        eps_rel,
    ));

    assert!(is_near(
        m.hrs()[&[3, -1, -1]][[0, 0]],
        Complex64::new(0.007378 / 4.0, 0.0),
        eps_abs,
        eps_rel,
    ));
    assert!(is_near(
        m.hrs()[&[3, -1, -1]][[2, 1]],
        Complex64::new(-0.008540 / 4.0, 0.0),
        eps_abs,
        eps_rel,
    ));

    assert!(is_near(
        m.hrs()[&[-2, -2, 2]][[0, 0]],
        Complex64::new(0.011647 / 6.0, 0.0),
        eps_abs,
        eps_rel,
    ));
    assert!(is_near(
        m.hrs()[&[-2, -2, 2]][[2, 3]],
        Complex64::new(-0.003502 / 6.0, 0.0),
        eps_abs,
        eps_rel,
    ));
}
