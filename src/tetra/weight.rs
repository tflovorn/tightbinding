use ndarray::Array2;
use std::f64::MAX;
use tetra::dos::dos_contrib;
use tetra::grid::EnergyGrid;

/// Return an array with elements `w[[grid_index, band_index]]` which gives the
/// weight w_{kn} of each (k-point, band) pair to the expectation value
///
/// <X_n> = \sum_k X_n(k) w_{kn}
pub fn all_weights<G: EnergyGrid>(
    grid: &G,
    fermi: f64,
    use_curvature_correction: bool,
) -> Array2<f64> {
    let mut weights_kn = Array2::zeros([grid.points().len(), grid.bands()]);
    let energy = grid.energy();
    let (tetra_indices, vertex_diffs) = tetra_indices();

    for subcell in get_subcells(&grid.dims()) {
        // Iterate over all tetrahedra in this subcell.
        for tetra in tetra_indices.iter() {
            let vertices = get_tetra_vertices(&tetra, &vertex_diffs, &subcell);

            // Collect band energies at the vertices of this tetrahedron.
            let vertex_indices = [
                grid.grid_index(&vertices[0]),
                grid.grid_index(&vertices[1]),
                grid.grid_index(&vertices[2]),
                grid.grid_index(&vertices[3]),
            ];

            for band_index in 0..grid.bands() {
                let band_energies = [
                    energy[[vertex_indices[0], band_index]],
                    energy[[vertex_indices[1], band_index]],
                    energy[[vertex_indices[2], band_index]],
                    energy[[vertex_indices[3], band_index]],
                ];

                // For each band, sort vertices according to band energy.
                let (sorted_es, sorted_indices) = sort_vertices(band_energies, &vertex_indices);

                // Get the weight contribution for each vertex.
                let tetra_weights =
                    weight_contrib(grid, fermi, use_curvature_correction, &sorted_es);

                // Add the weight contributions to the total.
                for (weight, grid_index) in tetra_weights.iter().zip(sorted_indices.iter()) {
                    weights_kn[[*grid_index, band_index]] += *weight;
                }
            }
        }
    }

    weights_kn
}

fn get_subcells(dims: &[usize; 3]) -> Vec<[usize; 3]> {
    let mut subcells = Vec::with_capacity(dims[0] * dims[1] * dims[2]);

    for s_i2 in 0..dims[2] {
        for s_i1 in 0..dims[1] {
            for s_i0 in 0..dims[0] {
                subcells.push([s_i0, s_i1, s_i2]);
            }
        }
    }

    subcells
}

/// Return the vertices of each tetrahedron in a subcell.
///
/// Each cubic subcell is divided into six tetrahedra. Here we construct a list of these
/// tetrahedra, identified by their vertices (four of the eight corners of the cubic subcell).
/// The vertices are specified by their relation to one corner of the subcell.
fn tetra_indices() -> ([[usize; 4]; 6], [[usize; 3]; 8]) {
    let tetras = [
        [0, 1, 2, 5],
        [0, 2, 4, 5],
        [2, 4, 5, 6],
        [2, 5, 6, 7],
        [2, 3, 5, 7],
        [1, 2, 3, 5],
    ];
    let vertices = [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [1, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 1],
    ];

    (tetras, vertices)
}

/// Give the point coordinates `[i0, i1, i2]` for each of the vertices of the
/// given tetrahedron located inside the given subcell.
fn get_tetra_vertices(
    tetra: &[usize; 4],
    vertex_diffs: &[[usize; 3]; 8],
    subcell: &[usize; 3],
) -> [[usize; 3]; 4] {
    let mut vertices = [[0; 3]; 4];

    for (vertex_index, diff_index) in tetra.iter().enumerate() {
        let vertex_diff = vertex_diffs[*diff_index];

        vertices[vertex_index] = [
            subcell[0] + vertex_diff[0],
            subcell[1] + vertex_diff[1],
            subcell[2] + vertex_diff[2],
        ];
    }

    vertices
}

/// Given a set of vertices, return their energies in sorted order and
/// those vertices sorted in the same order.
fn sort_vertices(band_energies: [f64; 4], vertex_indices: &[usize; 4]) -> ([f64; 4], [usize; 4]) {
    let order = get_order(band_energies);

    let sorted_es = [
        band_energies[order[0]],
        band_energies[order[1]],
        band_energies[order[2]],
        band_energies[order[3]],
    ];

    let sorted_indices = [
        vertex_indices[order[0]],
        vertex_indices[order[1]],
        vertex_indices[order[2]],
        vertex_indices[order[3]],
    ];

    (sorted_es, sorted_indices)
}

/// Return the order in which xs is sorted. See tests::check_get_order() for
/// an example.
fn get_order(mut xs: [f64; 4]) -> [usize; 4] {
    let mut order = [0; 4];

    for i in 0..4 {
        let mut min = MAX;
        let mut min_index = 0;
        for j in 0..4 {
            if xs[j] < min {
                min_index = j;
                min = xs[j];
            }
        }
        order[i] = min_index;
        xs[min_index] = MAX;
    }

    order
}

/// Return the weight contributions of the tetrahedron vertices with energies given
/// in ascending order by `sorted_es`.
fn weight_contrib<G: EnergyGrid>(
    grid: &G,
    fermi: f64,
    use_curvature_correction: bool,
    sorted_es: &[f64; 4],
) -> [f64; 4] {
    let fac = grid.tetra_volume() / 4.0;
    let (e1, e2, e3, e4) = (sorted_es[0], sorted_es[1], sorted_es[2], sorted_es[3]);
    let mut ws = [0.0; 4];

    if fermi <= e1 {
        ws = [0.0; 4];
    } else if e1 <= fermi && fermi <= e2 {
        let c = fac * (fermi - e1).powi(3) / ((e2 - e1) * (e3 - e1) * (e4 - e1));

        ws[0] = c * (4.0 - (fermi - e1) * (1.0 / (e2 - e1) + 1.0 / (e3 - e1) + 1.0 / (e4 - e1)));
        ws[1] = c * (fermi - e1) / (e2 - e1);
        ws[2] = c * (fermi - e1) / (e3 - e1);
        ws[3] = c * (fermi - e1) / (e4 - e1);
    } else if e2 <= fermi && fermi <= e3 {
        let c_1 = fac * (fermi - e1).powi(2) / ((e4 - e1) * (e3 - e1));
        let c_2 =
            fac * (fermi - e1) * (fermi - e2) * (e3 - fermi) / ((e4 - e1) * (e3 - e2) * (e3 - e1));
        let c_3 = fac * (fermi - e2).powi(2) * (e4 - fermi) / ((e4 - e2) * (e3 - e2) * (e4 - e1));

        ws[0] = c_1
            + (c_1 + c_2) * (e3 - fermi) / (e3 - e1)
            + (c_1 + c_2 + c_3) * (e4 - fermi) / (e4 - e1);
        ws[1] = c_1
            + c_2
            + c_3
            + (c_2 + c_3) * (e3 - fermi) / (e3 - e2)
            + c_3 * (e4 - fermi) / (e4 - e2);
        ws[2] = (c_1 + c_2) * (fermi - e1) / (e3 - e1) + (c_2 + c_3) * (fermi - e2) / (e3 - e2);
        ws[3] = (c_1 + c_2 + c_3) * (fermi - e1) / (e4 - e1) + c_3 * (fermi - e2) / (e4 - e2);
    } else if e3 <= fermi && fermi <= e4 {
        let c = fac * (e4 - fermi).powi(3) / ((e4 - e1) * (e4 - e2) * (e4 - e3));

        ws[0] = fac - c * (e4 - fermi) / (e4 - e1);
        ws[1] = fac - c * (e4 - fermi) / (e4 - e2);
        ws[2] = fac - c * (e4 - fermi) / (e4 - e3);
        ws[3] =
            fac - c * (4.0 - (1.0 / (e4 - e1) + 1.0 / (e4 - e2) + 1.0 / (e4 - e3)) * (e4 - fermi));
    } else {
        // fermi > e4
        ws[0] = fac;
        ws[1] = fac;
        ws[2] = fac;
        ws[3] = fac;
    }

    if use_curvature_correction {
        let ccs = curvature_correction(grid, fermi, sorted_es);

        [
            ws[0] + ccs[0],
            ws[1] + ccs[1],
            ws[2] + ccs[2],
            ws[3] + ccs[3],
        ]
    } else {
        ws
    }
}

/// Return the curvature corrections for the tetrahedron vertices with
/// energies given in ascending order by `sorted_es`.
fn curvature_correction<G: EnergyGrid>(grid: &G, fermi: f64, sorted_es: &[f64; 4]) -> [f64; 4] {
    let dos_tetra = dos_contrib(grid, fermi, sorted_es);
    let sum_es: f64 = sorted_es.iter().sum();

    [
        (dos_tetra / 40.0) * (sum_es - 4.0 * sorted_es[0]),
        (dos_tetra / 40.0) * (sum_es - 4.0 * sorted_es[1]),
        (dos_tetra / 40.0) * (sum_es - 4.0 * sorted_es[2]),
        (dos_tetra / 40.0) * (sum_es - 4.0 * sorted_es[3]),
    ]
}

#[cfg(test)]
mod tests {
    use super::get_order;

    #[test]
    fn check_get_order() {
        let xs = [0.0, -5.0, 3.0, 2.0];
        assert_eq!(get_order(xs), [1, 0, 3, 2]);
    }
}
