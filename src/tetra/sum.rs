use tetra::grid::EvecGrid;

/// Calculate the contribution of one band to the expectation value of x.
///
/// <x_n> = \sum_k x_{nk} w_{nk}
pub fn band_expectation_value(x_k: &Vec<f64>, band_weights: &Vec<f64>) -> f64 {
    x_k.iter()
        .zip(band_weights)
        .map(|(x, weight)| x * weight)
        .sum()
}

/// Calculate the expectation value of x.
///
/// <x> = \sum_{nk} x_{nk} w_{nk}
pub fn expectation_value(x_nk: &Vec<Vec<f64>>, weights: &Vec<Vec<f64>>) -> f64 {
    x_nk.iter()
        .zip(weights)
        .map(|(band_x, band_weights)| {
            band_expectation_value(band_x, band_weights)
        })
        .sum()
}

/// Calculate the orbital-resolved contributions to the expectation value of x.
///
/// <x_i> = \sum_{nk} |U_{in}^k|^2 x_{nk} w_{nk}
pub fn orbital_expectation_values<G: EvecGrid>(
    grid: &G,
    x_nk: &Vec<Vec<f64>>,
    weights: &Vec<Vec<f64>>,
) -> Vec<f64> {
    let mut evs = vec![0.0; grid.bands()];

    for (band_index, (x_ks, w_ks)) in x_nk.iter().zip(weights).enumerate() {
        for (grid_index, (x, w)) in x_ks.iter().zip(w_ks).enumerate() {
            for orbital_index in 0..grid.bands() {
                let u = grid.evec(grid_index)[[orbital_index, band_index]];

                evs[orbital_index] += u.norm_sqr() * x * w;
            }
        }
    }

    evs
}

/// Operator elements n_{nk} = 1 for the number operator.
fn num_x(weights: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let mut x_nk = Vec::new();

    for band_weights in weights {
        x_nk.push(vec![1.0; band_weights.len()]);
    }

    x_nk
}

/// Calculate the contributions of each orbital to the total electron density <n>.
///
/// <n_i> = \sum_{nk} |U_{in}^k|^2 w_{nk}
pub fn orbital_number<G: EvecGrid>(grid: &G, weights: &Vec<Vec<f64>>) -> Vec<f64> {
    orbital_expectation_values(grid, &num_x(&weights), &weights)
}

/// Calculate the total electron density <n>.
///
/// <n> = \sum_{nk} w_{nk}
pub fn total_number(weights: &Vec<Vec<f64>>) -> f64 {
    expectation_value(&num_x(&weights), &weights)
}
