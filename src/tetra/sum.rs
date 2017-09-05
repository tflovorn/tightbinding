use ndarray::{Array1, Array2};
use tetra::grid::EvecGrid;

/// Calculate the contribution of one band to the expectation value of x.
///
/// <x_n> = \sum_k x_{kn} w_{kn}
///
/// TODO - may want to use Kahan sum here. For very large number of k-points,
/// expect loss of precision.
pub fn band_expectation_value(x_k: &Array1<f64>, w_k: &Array1<f64>) -> f64 {
    x_k.iter().zip(w_k.iter()).map(|(x, w)| x * w).sum()
}

/// Calculate the expectation value of x.
///
/// <x> = \sum_{kn} x_{kn} w_{kn}
pub fn expectation_value(x_kn: &Array2<f64>, w_kn: &Array2<f64>) -> f64 {
    x_kn.iter().zip(w_kn.iter()).map(|(x, w)| x * w).sum()
}

/// Calculate the orbital-resolved contributions to the expectation value of x.
///
/// <x_i> = \sum_{kn} |U_{in}^k|^2 x_{kn} w_{kn}
///
/// TODO - may want to use Kahan sum here. For very large number of k-points,
/// expect loss of precision.
pub fn orbital_expectation_values<G: EvecGrid>(
    grid: &G,
    x_kn: &Array2<f64>,
    w_kn: &Array2<f64>,
) -> Array1<f64> {
    let mut evs = Array1::zeros([grid.bands()]);

    for grid_index in 0..grid.points().len() {
        for orbital_index in 0..grid.bands() {
            for band_index in 0..grid.bands() {
                let u = grid.evec()[[grid_index, orbital_index, band_index]];

                evs[orbital_index] += u.norm_sqr() * x_kn[[grid_index, band_index]] *
                    w_kn[[grid_index, band_index]];
            }
        }
    }

    evs
}

/// Operator elements n_{kn} = 1 for the number operator.
fn num_x(num_ks: usize, num_bands: usize) -> Array2<f64> {
    Array2::from_elem([num_ks, num_bands], 1.0)
}

/// Calculate the contributions of each orbital to the total electron density <n>.
///
/// <n_i> = \sum_{kn} |U_{in}^k|^2 w_{kn}
pub fn orbital_number<G: EvecGrid>(grid: &G, weights: &Array2<f64>) -> Array1<f64> {
    let (num_ks, num_bands) = (weights.shape()[0], weights.shape()[1]);
    orbital_expectation_values(grid, &num_x(num_ks, num_bands), &weights)
}

/// Calculate the total electron density <n>.
///
/// <n> = \sum_{kn} w_{kn}
pub fn total_number(weights: &Array2<f64>) -> f64 {
    let (num_ks, num_bands) = (weights.shape()[0], weights.shape()[1]);
    expectation_value(&num_x(num_ks, num_bands), &weights)
}
