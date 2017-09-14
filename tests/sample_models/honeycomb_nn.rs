extern crate num_complex;
extern crate ndarray;
extern crate tightbinding;

use std::f64::consts::PI;
use std::collections::BTreeMap;
use num_complex::Complex64;
use ndarray::{Array2, arr2};
use tightbinding::Model;

/// One-band tight-binding model on the honeycomb lattice with uniform
/// nearest-neighbor hopping.
#[derive(Clone)]
pub struct HoneycombNNModel {
    hrs: BTreeMap<[i32; 3], Array2<Complex64>>,
    d: Array2<f64>,
}

impl HoneycombNNModel {
    #[allow(dead_code)]
    pub fn new(t: f64, a: f64) -> HoneycombNNModel {
        let mut hrs = BTreeMap::new();

        let zero = Complex64::new(0.0, 0.0);
        let mt = Complex64::new(-t, 0.0);

        hrs.insert([0, 0, 0], arr2(&[[zero, mt], [mt, zero]]));
        hrs.insert([0, 1, 0], arr2(&[[zero, zero], [mt, zero]]));
        hrs.insert([1, 1, 0], arr2(&[[zero, zero], [mt, zero]]));
        hrs.insert([0, -1, 0], arr2(&[[zero, mt], [zero, zero]]));
        hrs.insert([-1, -1, 0], arr2(&[[zero, mt], [zero, zero]]));

        let d = arr2(
            &[
                [(1.0 / 2.0) * a, (1.0 / 2.0) * a, 0.0],
                [(-3.0_f64.sqrt() / 2.0) * a, (3.0_f64.sqrt() / 2.0) * a, 0.0],
                [0.0, 0.0, 10.0 * a],
            ],
        );

        HoneycombNNModel { hrs, d }
    }

    #[allow(dead_code)]
    pub fn hk_lat(t: f64, k_lat: &[f64; 3]) -> Array2<Complex64> {
        let mut hk = Array2::zeros((2, 2));

        hk[[0, 1]] = t *
            (1.0 + (Complex64::i() * 2.0 * PI * k_lat[1]).exp() +
                 (Complex64::i() * 2.0 * PI * (k_lat[0] + k_lat[1])).exp());
        hk[[1, 0]] = hk[[0, 1]].conj();

        hk
    }
}

impl Model for HoneycombNNModel {
    fn hrs(&self) -> &BTreeMap<[i32; 3], Array2<Complex64>> {
        &self.hrs
    }

    fn bands(&self) -> usize {
        2
    }

    fn d(&self) -> &Array2<f64> {
        &self.d
    }
}
