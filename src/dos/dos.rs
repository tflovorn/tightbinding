use ndarray::Array1;
use rayon::prelude::*;

use vec_util::transpose_vecs;
use tetra::{EvecGrid, all_weights, orbital_number};

pub fn dos_from_num<G: Sync + EvecGrid>(
    grid: &G,
    num_energies: usize,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let (min_e, max_e) = grid.energy_bounds();

    let es = Array1::linspace(min_e, max_e, num_energies)
        .as_slice()
        .unwrap()
        .to_vec();
    let use_curvature_correction = true;

    let weights = es.par_iter()
        .map(|e| all_weights(grid, *e, use_curvature_correction))
        .collect::<Vec<Vec<Vec<f64>>>>();

    let orbital_nums = transpose_vecs(
        &(weights
              .par_iter()
              .map(|w| orbital_number(grid, w))
              .collect()),
    );

    // DOS(E_i) \approx (n(E_i) - n(E_{i - 1})) / (E_i - E_{i - 1})
    let mut orbital_dos = Vec::with_capacity(grid.bands());
    for num in orbital_nums.iter() {
        let mut dos = Vec::with_capacity(num_energies - 2);

        for (i, e_i) in es.iter().enumerate() {
            if i == 0 || i == num_energies - 1 {
                continue;
            }
            let delta: f64 = e_i - es[i - 1];

            dos.push(0.5 * (num[i + 1] - num[i - 1]) / delta);
        }

        orbital_dos.push(dos);
    }

    (es[1..num_energies - 1].to_vec(), orbital_dos)
}
