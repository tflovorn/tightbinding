extern crate num_complex;
extern crate ndarray;
extern crate serde_json;
extern crate tightbinding;

use std::path::PathBuf;
use std::fs::File;
use std::f64::consts::PI;
use num_complex::Complex64;
use ndarray::Array2;
use tightbinding::float::{is_near_complex, is_near_float};
use tightbinding::Model;
use tightbinding::sample_models::SimpleCubicNNModel;
use tightbinding::fourier::{hk_lat, hk_cart};
use tightbinding::tetra::{EnergyGrid, EvecCache, grid_index, grid_k};
use tightbinding::tetra::find_fermi;
use tightbinding::dos::{DosValues, dos_from_num};

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

    for i2 in 0..dims[2] + 1 {
        for i1 in 0..dims[1] + 1 {
            for i0 in 0..dims[0] + 1 {
                let point = [i0, i1, i2];
                let index = grid_index(&point, &dims);
                let k_lat = grid_k(&point, &dims, &k_start, &k_stop);

                assert!(is_near_float(
                    cache.energy(index)[0],
                    SimpleCubicNNModel::epsilon(t, &k_lat),
                    eps_abs,
                    eps_rel,
                ))
            }
        }
    }

    let mid_energy = 0.0;
    let mid_occupation = 0.5;
    let fermi_mid = find_fermi(&cache, mid_occupation);
    assert!(is_near_float(fermi_mid, mid_energy, eps_abs, eps_rel));

    let use_curvature_correction = true;
    let num_energies = 10;
    let dos = dos_from_num(&cache, num_energies, use_curvature_correction);

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
fn sc_nn_full_dos() {
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

    let use_curvature_correction = false;
    let num_energies = 1001;
    let dos = dos_from_num(&cache, num_energies, use_curvature_correction);

    let out_path = PathBuf::from("plot/test_sc_nn_dos.json");
    let out_file = File::create(out_path).expect(
        "Error creating output file: sc_nn_full_dos expects to be run in crate root.",
    );
    serde_json::to_writer(&out_file, &dos).expect("Error writing output file");

    let expected_path = PathBuf::from("test_data/sc_nn/test_sc_nn_dos.json");
    let expected_file = File::open(expected_path).expect(
        "Error creating output file: sc_nn_full_dos expects to be run in crate root.",
    );
    let expected_dos: DosValues =
        serde_json::from_reader(&expected_file).expect("Error reading expected data file");

    for (e_out, e_expected) in dos.es.iter().zip(expected_dos.es.iter()) {
        assert!(is_near_float(*e_out, *e_expected, eps_abs_e, eps_rel));
    }

    for (band_dos_out, band_dos_expected) in
        dos.orbital_dos.iter().zip(expected_dos.orbital_dos.iter())
    {
        for (dos_out, dos_expected) in band_dos_out.iter().zip(band_dos_expected.iter()) {
            assert!(is_near_float(*dos_out, *dos_expected, eps_abs_dos, eps_rel));
        }
    }

    for (dos_out, dos_expected) in dos.total_dos.iter().zip(expected_dos.total_dos.iter()) {
        assert!(is_near_float(*dos_out, *dos_expected, eps_abs_dos, eps_rel));
    }
}
