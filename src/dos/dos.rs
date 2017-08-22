use ndarray::Array1;
use rayon::prelude::*;

use vec_util::transpose_vecs;
use tetra::{EvecGrid, all_weights, orbital_number};

pub fn dos_from_num<G: Sync + EvecGrid>(
    grid: &G,
    num_energies: usize,
    use_curvature_correction: bool,
) -> (Vec<f64>, Vec<Vec<f64>>, Vec<f64>) {
    let (min_e, max_e) = grid.energy_bounds();

    let es = Array1::linspace(min_e, max_e, num_energies)
        .as_slice()
        .unwrap()
        .to_vec();

    let orbital_nums = transpose_vecs(
        &(es.par_iter()
              .map(|e| {
            let w = all_weights(grid, *e, use_curvature_correction);
            orbital_number(grid, &w)
        })
              .collect()),
    );

    // Use lowest-order central difference for DOS:
    // DOS(E_i) \approx 0.5 * (n(E_{i + 1}) - n(E_{i - 1})) / (E_i - E_{i - 1}).
    // Since we use this, we have two fewer DOS points than energies.
    let num_dos = num_energies - 2;

    let mut orbital_dos = Vec::with_capacity(grid.bands());
    for num in orbital_nums.iter() {
        let mut dos = Vec::with_capacity(num_dos);

        for (i, e_i) in es.iter().enumerate() {
            if i == 0 || i == num_energies - 1 {
                continue;
            }
            let delta: f64 = e_i - es[i - 1];

            dos.push(0.5 * (num[i + 1] - num[i - 1]) / delta);
        }

        orbital_dos.push(dos);
    }

    let mut total_dos = vec![0.0; num_dos];
    for band_dos in orbital_dos.iter() {
        for (e_index, e_dos) in band_dos.iter().enumerate() {
            total_dos[e_index] += *e_dos;
        }
    }

    (es[1..num_energies - 1].to_vec(), orbital_dos, total_dos)
}
