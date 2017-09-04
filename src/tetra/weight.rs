use std::f64::MAX;
use vec_util::transpose_vecs;
use tetra::grid::EnergyGrid;
use tetra::dos::dos_contrib;

/// Return a list with elements `w[band_index][grid_index]` which gives the
/// weight w_{nk} of each (band, k-point) pair to the expectation value
///
/// <X_n> = \sum_k X_n(k) w_{nk}
pub fn all_weights<G: EnergyGrid>(
    grid: &G,
    fermi: f64,
    use_curvature_correction: bool,
) -> Vec<Vec<f64>> {
    let mut weights_kn = Vec::new();

    for point in grid.points() {
        weights_kn.push(band_weights(grid, fermi, use_curvature_correction, &point));
    }

    transpose_vecs(&weights_kn)
}

/// Return a list of the weights w_{nk}, which give the contribution of each band
/// at the given k-point to the expectation value
///
/// <X_n> = \sum_k X_n(k) w_{nk}
pub fn band_weights<G: EnergyGrid>(
    grid: &G,
    fermi: f64,
    use_curvature_correction: bool,
    point: &[usize; 3],
) -> Vec<f64> {
    let (tetra_indices, vertex_diffs) = tetra_indices();
    let neighbors = subcell_neighbors(point, &grid.dims());

    let mut weights = vec![0.0; grid.bands()];

    // Iterate over all subcells neighboring the point specified by `point`.
    for neighbor in neighbors {
        // Iterate over all tetrahedra in this subcell.
        for tetra in tetra_indices.iter() {
            let vertices = get_tetra_vertices(&tetra, &vertex_diffs, &neighbor);

            // Skip tetrahedra that do not include `point` as a vertex.
            if !contains_point(&vertices, &point) {
                continue;
            }

            // Collect band energies at the vertices of this tetrahedron.
            let vertex_energies = [
                grid.energy(grid.grid_index(&vertices[0])),
                grid.energy(grid.grid_index(&vertices[1])),
                grid.energy(grid.grid_index(&vertices[2])),
                grid.energy(grid.grid_index(&vertices[3])),
            ];

            for band_index in 0..grid.bands() {
                // For each band, sort vertices according to band energy.
                let (sorted_es, sorted_vs) = sort_vertices(&vertex_energies, &vertices, band_index);

                // Get the weight contribution for each vertex.
                let tetra_weights =
                    weight_contrib(grid, fermi, use_curvature_correction, &sorted_es);

                // Add the weight contribution from the vertex `point` to the total.
                let point_index = find_point_index(&sorted_vs, &point).unwrap();
                weights[band_index] += tetra_weights[point_index];
            }
        }
    }

    weights
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

/// Return a list of points identifying each cubic subcell which neighbors the subcell
/// at the given point.
fn subcell_neighbors(point: &[usize; 3], dims: &[usize; 3]) -> Vec<[usize; 3]> {
    let mut neighbors = Vec::new();

    let bounds = vec![
        [dims[0], dims[1], dims[2]],
        [0, dims[1], dims[2]],
        [dims[0], 0, dims[2]],
        [0, 0, dims[2]],
        [dims[0], dims[1], 0],
        [0, dims[1], 0],
        [dims[0], 0, 0],
        [0, 0, 0],
    ];
    // usize can't be negative, but we can subtract one usize from another.
    // Give list of (-difference) values here, which are all non-negative.
    let mdiffs = vec![
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [1, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 1],
    ];

    for (bound, mdiff) in bounds.iter().zip(mdiffs) {
        if (0..3).all(|i| point[i] != bound[i]) {
            let neighbor = (0..3).map(|i| point[i] - mdiff[i]).collect::<Vec<usize>>();

            neighbors.push([neighbor[0], neighbor[1], neighbor[2]]);
        }
    }

    neighbors
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

/// Return true iff the tetrahedron identified by tetra, inside the given subcell,
/// contains the given point as a vertex.
fn contains_point(vertices: &[[usize; 3]; 4], point: &[usize; 3]) -> bool {
    for vertex in vertices {
        if point.iter().zip(vertex).all(|(pi, vi)| pi == vi) {
            return true;
        }
    }

    false
}

/// Given a set of vertices, return their energies in sorted order and
/// those vertices sorted in the same order.
fn sort_vertices(
    vertex_energies: &[&Vec<f64>; 4],
    vertices: &[[usize; 3]; 4],
    band_index: usize,
) -> ([f64; 4], [[usize; 3]; 4]) {
    let band_energies = [
        vertex_energies[0][band_index],
        vertex_energies[1][band_index],
        vertex_energies[2][band_index],
        vertex_energies[3][band_index],
    ];

    let order = get_order(band_energies);

    let sorted_es = [
        band_energies[order[0]],
        band_energies[order[1]],
        band_energies[order[2]],
        band_energies[order[3]],
    ];

    let sorted_vs = [
        vertices[order[0]],
        vertices[order[1]],
        vertices[order[2]],
        vertices[order[3]],
    ];

    (sorted_es, sorted_vs)
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

/// Return Some(v_index), where sorted_vs[v_index] == point, if such a v_index exists;
/// otherwise return None.
fn find_point_index(sorted_vs: &[[usize; 3]; 4], point: &[usize; 3]) -> Option<usize> {
    for (v_index, vertex) in sorted_vs.iter().enumerate() {
        if point.iter().zip(vertex.iter()).all(|(pi, vi)| pi == vi) {
            return Some(v_index);
        }
    }

    None
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
        let c_2 = fac * (fermi - e1) * (fermi - e2) * (e3 - fermi) /
            ((e4 - e1) * (e3 - e2) * (e3 - e1));
        let c_3 = fac * (fermi - e2).powi(2) * (e4 - fermi) / ((e4 - e2) * (e3 - e2) * (e4 - e1));

        ws[0] = c_1 + (c_1 + c_2) * (e3 - fermi) / (e3 - e1) +
            (c_1 + c_2 + c_3) * (e4 - fermi) / (e4 - e1);
        ws[1] = c_1 + c_2 + c_3 + (c_2 + c_3) * (e3 - fermi) / (e3 - e2) +
            c_3 * (e4 - fermi) / (e4 - e2);
        ws[2] = (c_1 + c_2) * (fermi - e1) / (e3 - e1) + (c_2 + c_3) * (fermi - e2) / (e3 - e2);
        ws[3] = (c_1 + c_2 + c_3) * (fermi - e1) / (e4 - e1) + c_3 * (fermi - e2) / (e4 - e2);
    } else if e3 <= fermi && fermi <= e4 {
        let c = fac * (e4 - fermi).powi(3) / ((e4 - e1) * (e4 - e2) * (e4 - e3));

        ws[0] = fac - c * (e4 - fermi) / (e4 - e1);
        ws[1] = fac - c * (e4 - fermi) / (e4 - e2);
        ws[2] = fac - c * (e4 - fermi) / (e4 - e3);
        ws[3] = fac -
            c * (4.0 - (1.0 / (e4 - e1) + 1.0 / (e4 - e2) + 1.0 / (e4 - e3)) * (e4 - fermi));
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
