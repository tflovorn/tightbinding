use ndarray::Array2;
use num_complex::Complex64;
use std::f64::consts::PI;

use model::Model;

/// Calculates the Fourier transform H(k) of the tight-binding model H(R),
/// for a particular k given in lattice coordinates.
///
/// # Arguments
///
/// * `k_lat` - k in lattice coordinates
pub fn hk_lat<M: Model>(m: &M, k_lat: &[f64; 3]) -> Array2<Complex64> {
    let mut hk = Array2::<Complex64>::zeros((m.bands(), m.bands()));

    for (r_lat, hr) in m.hrs() {
        let k_dot_r = 2.0
            * PI
            * k_lat
                .iter()
                .zip(r_lat)
                .map(|(&ki, &ri)| ki * (ri as f64))
                .sum::<f64>();

        let coeff = Complex64::new(0.0, k_dot_r).exp();
        hk = hk + hr * coeff;
    }

    hk
}

/// Calculates the Fourier transform H(k) of the tight-binding model H(R),
/// for a particular k given in Cartesian coordinates.
///
/// # Arguments
///
/// * `k_cart` - k in Cartesian coordinates. The units of k_cart entries
/// are the inverse of the units of m.D entries.
pub fn hk_cart<M: Model>(m: &M, k_cart: &[f64; 3]) -> Array2<Complex64> {
    let k_cart_mat = Array2::<f64>::from_shape_vec((1, 3), k_cart.to_vec()).unwrap();
    let k_lat_mat = (1.0 / (2.0 * PI)) * k_cart_mat.dot(m.d());
    let k_lat = k_lat_mat.as_slice().unwrap();

    hk_lat(m, &[k_lat[0], k_lat[1], k_lat[2]])
}
