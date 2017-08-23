use std::f64::consts::PI;
use std::collections::HashMap;
use num_complex::Complex64;
use ndarray::Array2;
use tightbinding::Model;

/// One-band tight-binding model on the simple cubic lattice with uniform
/// nearest-neighbor hopping.
/// Has the spectrum
///     \epsilon(k) = -2 * t * (cos(k_x a) + cos(k_y a) + cos(k_z a))
#[derive(Clone)]
pub struct SimpleCubicNNModel {
    hrs: HashMap<[i32; 3], Array2<Complex64>>,
    d: Array2<f64>,
}

impl SimpleCubicNNModel {
    pub fn new(t: f64, a: f64) -> SimpleCubicNNModel {
        let mut hrs = HashMap::new();

        let mt = Complex64::new(-t, 0.0);

        hrs.insert([1, 0, 0], Array2::eye(1) * mt);
        hrs.insert([-1, 0, 0], Array2::eye(1) * mt);
        hrs.insert([0, 1, 0], Array2::eye(1) * mt);
        hrs.insert([0, -1, 0], Array2::eye(1) * mt);
        hrs.insert([0, 0, 1], Array2::eye(1) * mt);
        hrs.insert([0, 0, -1], Array2::eye(1) * mt);

        let d = Array2::eye(3) * a;

        SimpleCubicNNModel { hrs, d }
    }

    pub fn epsilon(t: f64, k_lat: &[f64; 3]) -> f64 {
        -2.0 * t *
            ((2.0 * PI * k_lat[0]).cos() + (2.0 * PI * k_lat[1]).cos() +
                (2.0 * PI * k_lat[2]).cos())
    }
}

impl Model for SimpleCubicNNModel {
    fn hrs(&self) -> &HashMap<[i32; 3], Array2<Complex64>> {
        &self.hrs
    }

    fn bands(&self) -> usize {
        1
    }

    fn d(&self) -> &Array2<f64> {
        &self.d
    }
}
