extern crate tightbinding;

use std::f64::consts::PI;
use tightbinding::float::is_near_float;
use tightbinding::qe::DftBands;

#[test]
fn heusler_bands() {
    let bands_path = "test_data/NiHfSn/NiHfSn_bulk_soc_bands.dat";

    let alat = 11.4746062404244;

    let bands = DftBands::new(bands_path, alat).unwrap();

    assert_eq!(bands.ks().len(), 241);
    assert_eq!(bands.emks().len(), 110);

    let eps_abs_k = 1e-12; // Angstrom^{-1}
    let eps_rel_k = 1e-12;

    let expected_k_indices = [1, 239];
    let expected_ks = [[0.0, 0.05, 0.0], [0.0125, 1.0, 0.0125]];

    for (k_index, expected_k) in expected_k_indices.iter().zip(expected_ks.iter()) {
        for (c, kc) in bands.ks()[*k_index].iter().enumerate() {
            assert!(is_near_float(
                *kc,
                (2.0 * PI / alat) * expected_k[c],
                eps_abs_k,
                eps_rel_k,
            ));
        }
    }

    let eps_abs_e = 1e-12; // eV
    let eps_rel_e = 1e-12;

    let m_index = 31;
    let expected_e_k_indices = [0, 1, 239, 240];
    let expected_es = [13.038, 13.016, 12.873, 12.872];

    for (k_index, expected_e) in expected_e_k_indices.iter().zip(expected_es.iter()) {
        assert!(is_near_float(
            bands.emks()[m_index][*k_index],
            *expected_e,
            eps_abs_e,
            eps_rel_e,
        ));
    }
}
