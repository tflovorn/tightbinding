use ndarray::Array1;
use rayon::prelude::*;

use vec_util::transpose_vecs;
use model::Model;
use tetra::{EvecCache, EnergyGrid, all_weights, orbital_number};

pub fn dos_from_num<M: Model + Sync>(
    m: &M,
    num_energies: usize,
    dims: [usize; 3],
    k_start: [f64; 3],
    k_stop: [f64; 3],
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let cache = EvecCache::new(m.clone(), dims, k_start, k_stop);
    let (min_e, max_e) = cache.energy_bounds();

    let es = Array1::linspace(min_e, max_e, num_energies)
        .as_slice()
        .unwrap()
        .to_vec();
    let use_curvature_correction = true;

    let weights = es.par_iter()
        .map(|e| all_weights(&cache, *e, use_curvature_correction))
        .collect::<Vec<Vec<Vec<f64>>>>();

    let orbital_nums = transpose_vecs(
        &(weights
              .par_iter()
              .map(|w| orbital_number(&cache, w))
              .collect()),
    );

    // DOS(E_i) \approx (n(E_i) - n(E_{i - 1})) / (E_i - E_{i - 1})
    let mut orbital_dos = Vec::with_capacity(m.bands());
    for num in orbital_nums.iter() {
        let mut dos = Vec::with_capacity(num_energies - 1);

        for (i, e_i) in es.iter().enumerate() {
            if i == 0 {
                continue;
            }
            let delta: f64 = e_i - es[i - 1];

            dos.push((num[i] - num[i - 1]) / delta);
        }

        orbital_dos.push(dos);
    }

    (es[1..num_energies].to_vec(), orbital_dos)
}
