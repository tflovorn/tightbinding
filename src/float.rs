use num_complex::Complex64;
use std::cmp::Ordering;

pub fn is_near_float(x: f64, y: f64, eps_abs: f64, eps_rel: f64) -> bool {
    let diff = (x - y).abs();

    if diff < eps_abs {
        return true;
    }

    diff < eps_rel * x.abs().max(y.abs())
}

pub fn is_near_complex(x: Complex64, y: Complex64, eps_abs: f64, eps_rel: f64) -> bool {
    let diff = (x - y).norm();

    if diff < eps_abs {
        return true;
    }

    diff < eps_rel * x.norm().max(y.norm())
}

/// Floating point number which panics if compared to a NaN.
/// Implementation taken from:
/// https://stackoverflow.com/questions/28247990/
/// how-to-do-a-binary-search-on-a-vec-of-floats/28248065#28248065
#[derive(PartialEq, PartialOrd, Copy, Clone)]
pub struct NonNan(f64);

impl NonNan {
    pub fn new(val: f64) -> Option<NonNan> {
        if val.is_nan() {
            None
        } else {
            Some(NonNan(val))
        }
    }

    pub fn val(&self) -> f64 {
        self.0
    }
}

impl Eq for NonNan {}

impl Ord for NonNan {
    fn cmp(&self, other: &NonNan) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}
