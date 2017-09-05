extern crate num_complex;
extern crate ndarray;
extern crate serde_json;
extern crate tightbinding;

use std::f64::consts::PI;
use num_complex::Complex64;
use ndarray::Array2;
use tightbinding::float::{is_near_complex, is_near_float};
use tightbinding::Model;
use tightbinding::fourier::{hk_lat, hk_cart};
use tightbinding::tetra::{KGrid, EnergyGrid, EvecCache};
use tightbinding::tetra::find_fermi;
use tightbinding::dos::dos_from_num;

mod sample_models;
use sample_models::SimpleCubicNNModel;

mod dos_util;
use dos_util::{write_dos_out, read_dos, check_dos};

#[test]
fn cubic_nn() {
    let t = 1.0;
    let a = 1.0;

    let m = SimpleCubicNNModel::new(t, a);

    let k_lat = [1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0];
    let k_cart = [PI / a, PI / a, PI / a];

    let hk = hk_cart(&m, &k_cart);
    let expected = Complex64::new(SimpleCubicNNModel::epsilon(t, &k_lat), 0.0);

    let eps_abs = 1e-12 * t;
    let eps_rel = 1e-12;

    assert!(is_near_complex(hk[[0, 0]], expected, eps_abs, eps_rel));

    let dims = [4, 4, 4];
    let k_start = [0.0, 0.0, 0.0];
    let k_stop = [1.0, 1.0, 1.0];

    let hk_fn = |k| hk_lat(&m, &k);

    let cache = EvecCache::new(hk_fn, m.bands(), dims, k_start, k_stop);

    for (grid_index, k) in cache.ks().iter().enumerate() {
        assert!(is_near_float(
            cache.energy()[[grid_index, 0]],
            SimpleCubicNNModel::epsilon(t, k),
            eps_abs,
            eps_rel,
        ))
    }

    let mid_energy = 0.0;
    let mid_occupation = 0.5;
    let fermi_mid = find_fermi(&cache, mid_occupation);
    assert!(is_near_float(fermi_mid, mid_energy, eps_abs, eps_rel));

    let num_energies = 10;
    let dos = dos_from_num(&cache, num_energies);

    let expected_dos = vec![
        0.018518518518518514,
        0.060185185185185154,
        0.12268518518518509,
        0.17129629629629625,
        0.17129629629629647,
        0.12268518518518522,
        0.060185185185185175,
        0.01851851851851851,
    ];

    let eps_abs_dos = 1e-12 / t;
    let eps_rel_dos = 1e-12;

    for (x, y) in dos.total_dos.iter().zip(expected_dos) {
        assert!(is_near_float(*x, y, eps_abs_dos, eps_rel_dos));
    }
}

#[test]
#[ignore]
fn sc_nn_dos() {
    let t = 1.0;
    let bands = 1;

    let eps_abs_e = 1e-12 * t;
    let eps_abs_dos = 1e-12 / t;
    let eps_rel = 1e-12;

    let dims = [32, 32, 32];
    let k_start = [0.0, 0.0, 0.0];
    let k_stop = [1.0, 1.0, 1.0];

    let hk_fn = |k| Array2::eye(1) * SimpleCubicNNModel::epsilon(t, &k);

    let cache = EvecCache::new(hk_fn, bands, dims, k_start, k_stop);

    let num_energies = 1001;
    let dos = dos_from_num(&cache, num_energies);

    let prefix = "sc_nn";
    write_dos_out(&dos, &prefix);
    let expected_dos = read_dos(&prefix);

    check_dos(
        &dos,
        &expected_dos,
        eps_abs_e,
        eps_rel,
        eps_abs_dos,
        eps_rel,
    );
}
