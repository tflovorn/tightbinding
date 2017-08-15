extern crate num_complex;
extern crate tightbinding;

use num_complex::Complex64;
use tightbinding::float::{is_near_float, is_near_complex};

#[test]
fn is_near() {
    let eps_abs = 1e-12;
    let eps_rel = 1e-12;

    assert!(is_near_float(1.0, 1.0, eps_abs, eps_rel));
    assert!(is_near_float(1.0, 1.0 + 1e-13, eps_abs, eps_rel));
    assert!(!is_near_float(1.0, 1.0 + 1e-11, eps_abs, eps_rel));

    let x = Complex64::new(1.0, 0.0);
    let d1 = Complex64::new(1e-13, 0.0);
    let d2 = Complex64::new(1e-11, 0.0);

    assert!(is_near_complex(x, x, eps_abs, eps_rel));
    assert!(is_near_complex(x, x + d1, eps_abs, eps_rel));
    assert!(!is_near_complex(x, x + d2, eps_abs, eps_rel));

    assert!(is_near_complex(
        x,
        x + Complex64::i() * d1,
        eps_abs,
        eps_rel,
    ));
    assert!(!is_near_complex(
        x,
        x + Complex64::i() * d2,
        eps_abs,
        eps_rel,
    ));
}
