extern crate num_complex;
extern crate ndarray;
extern crate tightbinding;

use std::f64::consts::PI;
use std::collections::HashMap;
use num_complex::Complex64;
use ndarray::Array2;
use tightbinding::float::{is_near_complex, is_near_float};
use tightbinding::Model;
use tightbinding::fourier::hk_cart;
use tightbinding::tetra::{EnergyGrid, EvecCache, grid_index, grid_k};
use tightbinding::tetra::find_fermi;

/// One-band tight-binding model on the cubic lattice with uniform
/// nearest-neighbor hopping.
/// Has the spectrum
///     \epsilon(k) = -2 * t * (cos(k_x a) + cos(k_y a) + cos(k_z a))
#[derive(Clone)]
struct CubicNNModel {
    hrs: HashMap<[i32; 3], Array2<Complex64>>,
    d: Array2<f64>,
}

impl CubicNNModel {
    pub fn new(t: f64, a: f64) -> CubicNNModel {
        let mut hrs = HashMap::new();

        let mt = Complex64::new(-t, 0.0);

        hrs.insert([1, 0, 0], Array2::eye(1) * mt);
        hrs.insert([-1, 0, 0], Array2::eye(1) * mt);
        hrs.insert([0, 1, 0], Array2::eye(1) * mt);
        hrs.insert([0, -1, 0], Array2::eye(1) * mt);
        hrs.insert([0, 0, 1], Array2::eye(1) * mt);
        hrs.insert([0, 0, -1], Array2::eye(1) * mt);

        let d = Array2::eye(1) * a;

        CubicNNModel { hrs, d }
    }

    pub fn epsilon(t: f64, a: f64, k_cart: &[f64]) -> f64 {
        assert!(k_cart.len() == 3);
        -2.0 * t * ((k_cart[0] * a).cos() + (k_cart[1] * a).cos() + (k_cart[2] * a).cos())
    }
}

impl Model for CubicNNModel {
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

#[test]
fn cubic_nn() {
    let t = 1.0;
    let a = 1.0;

    let m = CubicNNModel::new(t, a);

    let k_cart = [PI / a, PI / a, PI / a];

    let hk = hk_cart(&m, &k_cart);
    let expected = Complex64::new(CubicNNModel::epsilon(t, a, &k_cart), 0.0);

    let eps_abs = 1e-12 * t;
    let eps_rel = 1e-12;

    assert!(is_near_complex(hk[[0, 0]], expected, eps_abs, eps_rel));

    let dims = [4, 4, 4];
    let k_start = [0.0, 0.0, 0.0];
    let k_stop = [1.0, 1.0, 1.0];

    let cache = EvecCache::new(m.clone(), dims, k_start, k_stop);

    for i0 in 0..dims[0] + 1 {
        for i1 in 0..dims[1] + 1 {
            for i2 in 0..dims[2] + 1 {
                let point = [i0, i1, i2];
                let index = grid_index(&point, &dims);
                let k_cart = grid_k(&point, &dims, &k_start, &k_stop)
                    .iter()
                    .map(|ki| (2.0 * PI / a) * ki)
                    .collect::<Vec<f64>>();

                assert!(is_near_float(
                    cache.energy(index)[0],
                    CubicNNModel::epsilon(t, a, &k_cart),
                    eps_abs,
                    eps_rel,
                ))
            }
        }
    }

    let mid_energy = 0.0;
    let mid_occupation = 0.5;
    let fermi_mid = find_fermi(&cache, mid_occupation);
    println!("fermi {:?}", fermi_mid);
    assert!(is_near_float(fermi_mid, mid_energy, eps_abs, eps_rel));
}
