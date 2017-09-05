extern crate num_complex;
extern crate ndarray;
extern crate linxal;
extern crate tightbinding;

use ndarray::Axis;
use linxal::eigenvalues::SymEigen;
use linxal::types::Symmetric;
use tightbinding::float::is_near_float;
use tightbinding::Model;
use tightbinding::fourier::hk_lat;
use tightbinding::tetra::{KGrid, EnergyGrid, EvecCache};
use tightbinding::dos::dos_from_num;

mod sample_models;
use sample_models::HoneycombNNModel;

mod dos_util;
use dos_util::{write_dos_out, read_dos, check_dos};

#[test]
fn honeycomb_nn() {
    let t = 1.0;
    let a = 1.0;

    let m = HoneycombNNModel::new(t, a);

    let eps_abs = 1e-12 * t;
    let eps_rel = 1e-12;

    let dims = [4, 4, 1];
    let k_start = [0.0, 0.0, 0.0];
    let k_stop = [1.0, 1.0, 1.0];

    let hk_fn = |k| hk_lat(&m, &k);

    let cache = EvecCache::new(hk_fn, m.bands(), dims, k_start, k_stop);

    for (grid_index, k_lat) in cache.ks().iter().enumerate() {
        let hk = HoneycombNNModel::hk_lat(t, k_lat);

        let solution = SymEigen::compute(&hk, Symmetric::Upper, true).unwrap();
        let es = solution.values.to_vec();

        for (e_grid, e_expected) in
            cache.energy().subview(Axis(0), grid_index).iter().zip(
                es.iter(),
            )
        {
            assert!(is_near_float(*e_grid, *e_expected, eps_abs, eps_rel));
        }
    }
}

#[test]
#[ignore]
fn honeycomb_nn_dos() {
    let t = 1.0;
    let bands = 2;

    let eps_abs_e = 1e-12 * t;
    let eps_abs_dos = 1e-12 / t;
    let eps_rel = 1e-12;

    let dims = [32, 32, 1];
    let k_start = [0.0, 0.0, 0.0];
    let k_stop = [1.0, 1.0, 1.0];

    let hk_fn = |k| HoneycombNNModel::hk_lat(t, &k);

    let cache = EvecCache::new(hk_fn, bands, dims, k_start, k_stop);

    let num_energies = 1001;
    let dos = dos_from_num(&cache, num_energies);

    let prefix = "honeycomb_nn";
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
