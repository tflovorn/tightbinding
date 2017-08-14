use std::f64::consts::PI;
use num_complex::Complex64;
use rulinalg::matrix::Matrix;

use model::Model;

/// Calculates the Fourier transform H(k) of the tight-binding model H(R),
/// for a particular k given in lattice coordinates.
///
/// # Arguments
///
/// * `k_lat` - k in lattice coordinates
pub fn hk_lat<M>(m: M, k_lat: &[f64]) -> Matrix<Complex64>
    where M: Model
{
    let mut hk = Matrix::<Complex64>::zeros(m.bands(), m.bands());

    for (r_lat, hr) in m.hrs() {
        let k_dot_r = 2.0 * PI *
            k_lat
                .iter()
                .zip(r_lat)
                .map(|(&ki, &ri)| ki * (ri as f64))
                .sum::<f64>();

        let coeff = Complex64::new(0.0, k_dot_r).exp();
        hk += hr * coeff;
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
pub fn hk_cart<M>(m: M, k_cart: &[f64]) -> Matrix<Complex64>
    where M: Model
{
    let k_cart_mat = Matrix::<f64>::new(1, 3, k_cart);
    let k_lat_mat = m.d() * k_cart_mat * (1.0 / (2.0 * PI));
    let k_lat = k_lat_mat.data();

    hk_lat(m, &k_lat[..])
}
