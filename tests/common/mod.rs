extern crate num_complex;

use num_complex::Complex64;

pub fn is_near(x: Complex64, y: Complex64, eps_abs: f64, eps_rel: f64) -> bool {
    let diff = (x - y).norm();

    if diff < eps_abs {
        return true;
    }

    diff < eps_rel * x.norm().max(y.norm())
}
