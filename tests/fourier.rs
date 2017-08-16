extern crate num_complex;
extern crate rulinalg;
extern crate tightbinding;

use std::f64::consts::PI;
use std::collections::HashMap;
use num_complex::Complex64;
use rulinalg::matrix::Matrix;
use tightbinding::float::is_near_complex;
use tightbinding::Model;
use tightbinding::fourier::hk_cart;

/// One-band tight-binding model on the cubic lattice with uniform
/// nearest-neighbor hopping.
/// Has the spectrum
///     \epsilon(k) = -2 * t * (cos(k_x a) + cos(k_y a) + cos(k_z a))
struct CubicNNModel {
    hrs: HashMap<[i32; 3], Matrix<Complex64>>,
    d: Matrix<f64>,
}

impl CubicNNModel {
    pub fn new(t: f64, a: f64) -> CubicNNModel {
        let mut hrs = HashMap::new();

        let mt = Complex64::new(-t, 0.0);

        hrs.insert([1, 0, 0], Matrix::identity(1) * mt);
        hrs.insert([-1, 0, 0], Matrix::identity(1) * mt);
        hrs.insert([0, 1, 0], Matrix::identity(1) * mt);
        hrs.insert([0, -1, 0], Matrix::identity(1) * mt);
        hrs.insert([0, 0, 1], Matrix::identity(1) * mt);
        hrs.insert([0, 0, -1], Matrix::identity(1) * mt);

        let d = Matrix::identity(1) * a;

        CubicNNModel { hrs, d }
    }
}

impl Model for CubicNNModel {
    fn hrs(&self) -> &HashMap<[i32; 3], Matrix<Complex64>> {
        &self.hrs
    }

    fn bands(&self) -> usize {
        1
    }

    fn d(&self) -> &Matrix<f64> {
        &self.d
    }
}

#[test]
fn cubic_nn() {
    let t = 1.0;
    let a = 1.0;

    let m = CubicNNModel::new(t, a);

    let k_cart = [PI / a, PI / a, PI / a];

    let hk = hk_cart(&m, &k_cart);
    let expected = Complex64::new(6.0 * t, 0.0);

    let eps_abs = 1e-12 * t;
    let eps_rel = 1e-12;

    assert!(is_near_complex(hk[[0, 0]], expected, eps_abs, eps_rel));
}
