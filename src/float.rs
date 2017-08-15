use num_complex::Complex64;

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
