use num_complex::Complex64;
use ndarray::Array2;
use linxal::eigenvalues::SymEigen;
use linxal::types::Symmetric;
use itertools::multizip;

use model::Model;
use fourier::hk_lat;

/// Converts discrete k-point grid coordinates to associated energy values.
pub trait EnergyGrid {
    /// Band energies associated with the k-point at grid_index.
    fn energy(&self, grid_index: usize) -> &Vec<f64>;

    /// Dimensions of the grid along [k1, k2, k3] directions.
    ///
    /// The full grid ranges in 0..dims[i]+1 in each direction;
    /// the final point, at dims[i], is equivalent to the first point,
    /// at 0, when the grid covers the full Brillouin zone.
    fn dims(&self) -> [usize; 3];

    /// Volume of a tetrahedron as a fraction of the Brillouin zone volume.
    ///
    /// When the grid covers the full Brillouin zone, this is
    /// (1 / number of tetrahedra) = (1 / (dims[0] * dims[1] * dims[2])).
    fn tetra_volume(&self) -> f64;
}

/// Converts discrete k-point grid coordinates to associated eigenvectors.
/// The eigenvectors are sorted in order of the corresponding eigenenergies.
pub trait EvecGrid: EnergyGrid {
    /// Eigenvectors associated with the k-point at grid_index.
    fn evec(&self, grid_index: usize) -> &Array2<Complex64>;
}

/// Compute the linearized grid index associated with the given k-point index.
///
/// TODO - would passing copy of point, dims be more efficient?
pub fn grid_index(point: &[usize; 3], dims: &[usize; 3]) -> usize {
    point[0] + (dims[0] + 1) * (point[1] + point[2] * (dims[1] + 1))
}

pub fn grid_k(
    point: &[usize; 3],
    dims: &[usize; 3],
    k_start: &[f64; 3],
    k_stop: &[f64; 3],
) -> [f64; 3] {
    let mut k = [0.0, 0.0, 0.0];

    for (dim_index, (i, dim, start, stop)) in multizip((point, dims, k_start, k_stop)).enumerate() {
        let step = (stop - start) / (*dim as f64);
        k[dim_index] = start + (*i as f64) * step;
    }

    k
}

pub struct EvecCache<M: Model> {
    m: M,
    dims: [usize; 3],
    k_start: [f64; 3],
    k_stop: [f64; 3],
    energy: Vec<Vec<f64>>,
    evec: Vec<Array2<Complex64>>,
}

impl<M: Model> EvecCache<M> {
    /// Construct a new EvecCache from the given model, with the given number of k-points in
    /// each reciprocal lattice direction and occupying the given region of the Brillouin zone
    /// (to sample the full Brillouin zone, set k_start = [0.0, 0.0, 0.0] and
    /// k_stop = [1.0, 1.0, 1.0]).
    ///
    /// TODO convert k from original coordinates to coordinates with sign / order of G permuted to
    /// minimize tetrahedron length: "In order to minimize interpolation distances the shortest
    /// main diagonal is chosen." Ensure that k_start, k_stop are used correctly when this is done:
    /// permute elements / change sign of k_start, k_stop in the same way as reciprocal lattice
    /// vectors are permuted / changed sign.
    ///
    /// TODO consider symmetry operations which make some k-points redundant.
    /// Since Wannier90 makes no guarantees about preserving symmetry
    /// (unless symmetry-preservation is enabled, which forbids use of inner window),
    /// would we be able to take advantage of symmetry operations if they were implemented?
    ///
    /// TODO implement FFT calculation to compute all H(k) values at once.
    pub fn new(m: M, dims: [usize; 3], k_start: [f64; 3], k_stop: [f64; 3]) -> EvecCache<M> {
        let mut energy = Vec::new();
        let mut evec = Vec::new();

        for i0 in 0..dims[0] + 1 {
            for i1 in 0..dims[1] + 1 {
                for i2 in 0..dims[2] + 1 {
                    let k = grid_k(&[i0, i1, i2], &dims, &k_start, &k_stop);

                    let hk = hk_lat(&m, &k);

                    let solution = SymEigen::compute(&hk, Symmetric::Upper, true).unwrap();

                    energy.push(solution.values.to_vec());
                    evec.push(solution.right_vectors.unwrap());
                }
            }
        }

        EvecCache {
            m,
            dims,
            k_start,
            k_stop,
            energy,
            evec,
        }
    }
}

impl<M: Model> EnergyGrid for EvecCache<M> {
    fn energy(&self, grid_index: usize) -> &Vec<f64> {
        &self.energy[grid_index]
    }

    fn dims(&self) -> [usize; 3] {
        self.dims
    }

    fn tetra_volume(&self) -> f64 {
        let k_volume: f64 = self.k_start.iter().zip(&self.k_stop).map(|(start, stop)| stop - start).product();
        let num_tetra: f64 = self.dims.iter().map(|x| *x as f64).sum();

        1.0 / (k_volume * num_tetra)
    }
}

impl<M: Model> EvecGrid for EvecCache<M> {
    fn evec(&self, grid_index: usize) -> &Array2<Complex64> {
        &self.evec[grid_index]
    }
}
