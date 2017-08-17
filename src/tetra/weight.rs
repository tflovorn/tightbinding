use float::NonNan;
use vec_util::transpose_vecs;
use tetra::grid::{EnergyGrid, grid_index};
use tetra::dos::dos_contrib;

/// Return a list with elements `w[band_index][grid_index]` which gives the
/// weight w_{nk} of each (band, k-point) pair to the expectation value
///
/// <X_n> = \sum_k X_n(k) w_{nk}
pub fn all_weights<G: EnergyGrid>(grid: &G, fermi: f64) -> Vec<Vec<f64>> {
    let mut weights_kn = Vec::new();
    let dims = grid.dims();

    for i0 in 0..dims[0] {
        for i1 in 0..dims[1] {
            for i2 in 0..dims[2] {
                let point = [i0, i1, i2];

                weights_kn.push(band_weights(grid, fermi, &point));
            }
        }
    }

    transpose_vecs(&weights_kn)
}

/// Return a list of the weights w_{nk}, which give the contribution of each band
/// at the given k-point to the expectation value
///
/// <X_n> = \sum_k X_n(k) w_{nk}
pub fn band_weights<G: EnergyGrid>(grid: &G, fermi: f64, point: &[usize; 3]) -> Vec<f64> {
    let (tetra_indices, vertex_diffs) = tetra_indices();
    let neighbors = subcell_neighbors(point, &grid.dims());

    let mut weights = vec![0.0; grid.bands()];

    // Iterate over all subcells neighboring the subcell specified by `point`.
    for neighbor in neighbors {
        // Iterate over all tetrahedra in this subcell.
        for tetra in tetra_indices.iter() {
            let vertices = get_tetra_vertices(&tetra, &vertex_diffs, &neighbor);

            // Skip tetrahedra that do not include `point` as a vertex.
            if !contains_point(&vertices, &point) {
                continue;
            }

            // Collect band energies at the vertices of this tetrahedron.
            let vertex_energies: Vec<&Vec<f64>> = vertices
                .iter()
                .map(|v| grid_index(v, &grid.dims()))
                .map(|i| grid.energy(i))
                .collect();

            for band_index in 0..grid.bands() {
                // For each band, sort vertices according to band energy.
                let mut sorted_es_vs: Vec<(f64, &[usize; 3])> = vertex_energies
                    .iter()
                    .map(|es| es[band_index])
                    .zip(&vertices)
                    .collect();
                sorted_es_vs.sort_by_key(|&(e, _)| NonNan::new(e).unwrap());

                let sorted_es = sorted_es_vs.iter().map(|&(e, _)| e).collect();
                let sorted_vs = sorted_es_vs.iter().map(|&(_, v)| v).collect();

                // Get the weight contribution for each vertex.
                let tetra_weights = weight_contrib(grid, fermi, &sorted_es);

                // Add the weight contribution from the vertex `point` to the total.
                // TODO - for very dense grid, Kahan summation may be useful.
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
fn tetra_indices() -> (Vec<[usize; 4]>, Vec<[usize; 3]>) {
    let tetras = vec![
        [0, 1, 2, 5],
        [0, 2, 4, 5],
        [2, 4, 5, 6],
        [2, 5, 6, 7],
        [2, 3, 5, 7],
        [1, 2, 3, 5],
    ];
    let vertices = vec![
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
    vertex_diffs: &Vec<[usize; 3]>,
    subcell: &[usize; 3],
) -> Vec<[usize; 3]> {
    let mut vertices = Vec::new();

    for diff_index in tetra {
        let vertex_diff = vertex_diffs[*diff_index];

        let vertex = (0..3)
            .map(|i| subcell[i] + vertex_diff[i])
            .collect::<Vec<usize>>();

        vertices.push([vertex[0], vertex[1], vertex[2]]);
    }

    vertices
}

/// Return true iff the tetrahedron identified by tetra, inside the given subcell,
/// contains the given point as a vertex.
fn contains_point(vertices: &Vec<[usize; 3]>, point: &[usize; 3]) -> bool {
    for vertex in vertices {
        if point.iter().zip(vertex).all(|(pi, vi)| pi == vi) {
            return true;
        }
    }

    false
}

/// Return Some(v_index), where sorted_vs[v_index] == point, if such a v_index exists;
/// otherwise return None.
fn find_point_index(sorted_vs: &Vec<&[usize; 3]>, point: &[usize; 3]) -> Option<usize> {
    for (v_index, vertex) in sorted_vs.iter().enumerate() {
        if point.iter().zip(*vertex).all(|(pi, vi)| pi == vi) {
            return Some(v_index);
        }
    }

    None
}

/// Return the weight contributions of the tetrahedron vertices with energies given
/// by `sorted_es`.
fn weight_contrib<G: EnergyGrid>(grid: &G, fermi: f64, sorted_es: &Vec<f64>) -> Vec<f64> {
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

    let ccs = curvature_correction(grid, fermi, sorted_es);

    ws.iter().zip(ccs).map(|(w, cc)| w + cc).collect()
}

fn curvature_correction<G: EnergyGrid>(grid: &G, fermi: f64, sorted_es: &Vec<f64>) -> Vec<f64> {
    let dos_tetra = dos_contrib(grid, fermi, sorted_es);
    let sum_es: f64 = sorted_es.iter().sum();

    sorted_es
        .iter()
        .map(|e| (dos_tetra / 40.0) * (sum_es - e))
        .collect()
}
