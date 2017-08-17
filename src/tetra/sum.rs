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

/// Calculate the contribution of one band to the total electron density <n>.
///
/// <n_n> = \sum_{k} w_{nk}
pub fn band_number(band_weights: &Vec<f64>) -> f64 {
    let x_k = vec![1.0; band_weights.len()];

    band_expectation_value(&x_k, &band_weights)
}

/// Calculate the total electron density <n>.
///
/// <n> = \sum_{nk} w_{nk}
pub fn total_number(weights: &Vec<Vec<f64>>) -> f64 {
    let mut x_nk = Vec::new();

    for band_weights in weights {
        x_nk.push(vec![1.0; band_weights.len()]);
    }

    expectation_value(&x_nk, &weights)
}
