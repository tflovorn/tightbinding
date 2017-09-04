use std::f64;
use num_complex::Complex64;
use ndarray::{Array2, Array3, Axis};
use linxal::eigenvalues::SymEigen;
use linxal::eigenvalues::types::Solution;
use linxal::types::Symmetric;
use itertools::multizip;
use rayon::prelude::*;

use float::NonNan;

pub trait KGrid {
    /// Dimensions of the grid along [k1, k2, k3] directions.
    ///
    /// The full grid ranges in 0..dims[i]+1 in each direction;
    /// the final point, at dims[i], is equivalent to the first point,
    /// at 0, when the grid covers the full Brillouin zone.
    fn dims(&self) -> [usize; 3];

    /// Return a reference to the list of grid points.
    fn points(&self) -> &Vec<[usize; 3]>;

    /// Return the linearized grid index associated with the given k-point index.
    ///
    /// TODO - would passing copy of point, dims be more efficient?
    ///
    /// # Panics
    ///
    /// Panics if point is not a valid point on the grid, i.e. if point is not a value
    /// stored in points().
    fn grid_index(&self, point: &[usize; 3]) -> usize;

    /// Return a reference to the list of k-points, which are given in reciprocal
    /// lattice coordinates.
    fn ks(&self) -> &Vec<[f64; 3]>;

    /// Smallest and largest value of k included in each reciprocal lattice direction.
    ///
    /// When k_range = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), the full Brillouin zone is included.
    fn k_range(&self) -> ([f64; 3], [f64; 3]);

    /// Volume of a tetrahedron as a fraction of the sampled Brillouin zone volume.
    ///
    /// Equal to:
    /// (1 / number of tetrahedra) = 1 / (6 * (dims[0] * dims[1] * dims[2])).
    fn tetra_volume(&self) -> f64;
}

/// Converts discrete k-point grid coordinates to associated energy values.
pub trait EnergyGrid: KGrid {
    /// Band energies e[[grid_index, band_index]].
    fn energy(&self) -> &Array2<f64>;

    /// Largest and smallest energy values contained in the grid.
    fn energy_bounds(&self) -> (f64, f64) {
        let mut min = NonNan::new(f64::MAX).unwrap();
        let mut max = NonNan::new(f64::MIN).unwrap();

        for e_f64 in self.energy().iter() {
            let e = NonNan::new(*e_f64).unwrap();
            if e < min {
                min = e;
            }
            if e > max {
                max = e;
            }
        }

        (min.val(), max.val())
    }

    /// Number of bands at each k-point.
    fn bands(&self) -> usize;
}

/// Converts discrete k-point grid coordinates to associated eigenvectors.
/// The eigenvectors are sorted in order of the corresponding eigenenergies.
pub trait EvecGrid: EnergyGrid {
    /// Eigenvectors evec[[grid_index, orbital_index, band_index]].
    fn evec(&self) -> &Array3<Complex64>;
}

/// Compute the linearized grid index associated with the given k-point index.
///
/// TODO - would passing copy of point, dims be more efficient?
fn get_grid_index(point: &[usize; 3], dims: &[usize; 3]) -> usize {
    point[0] + (dims[0] + 1) * (point[1] + point[2] * (dims[1] + 1))
}

/// Convert the grid point to the corresponding k-point.
fn get_grid_k(
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

/// Cache storing the eigenvalues and eigenvectors of a model computed
/// on the k-point grid bounded by `k_start` and `k_stop` and with number of
/// k-points in each direction given by `dims`.
struct KCache {
    dims: [usize; 3],
    k_start: [f64; 3],
    k_stop: [f64; 3],
    tetra_volume: f64,
    points: Vec<[usize; 3]>,
    ks: Vec<[f64; 3]>,
}

impl KCache {
    /// Construct a new KCache with the given number of k-points in each reciprocal lattice
    /// direction and occupying the given region of the Brillouin zone.
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
    /// # Arguments
    ///
    /// * `dims` - the number of points to sample inside the Brillouin zone
    /// along each (lattice coordinate) direction.
    ///
    /// * `k_start` - smallest value of k (in lattice coordinates) to sample
    /// in each direction.
    ///
    /// * `k_start` - largest value of k (in lattice coordinates) to sample in
    /// each direction. To sample the whole Brillouin zone, give
    /// k_start = [0.0, 0.0, 0.0] and k_stop = [1.0, 1.0, 1.0].
    pub fn new(dims: [usize; 3], k_start: [f64; 3], k_stop: [f64; 3]) -> KCache {
        let end_grid_index = (dims[0] + 1) * (dims[1] + 1) * (dims[2] + 1);

        let mut points = Vec::with_capacity(end_grid_index);
        let mut ks = Vec::with_capacity(end_grid_index);

        let mut grid_index = 0;
        for i2 in 0..dims[2] + 1 {
            for i1 in 0..dims[1] + 1 {
                for i0 in 0..dims[0] + 1 {
                    let point = [i0, i1, i2];

                    points.push(point);
                    ks.push(get_grid_k(&point, &dims, &k_start, &k_stop));

                    grid_index += 1;
                }
            }
        }

        assert_eq!(grid_index, end_grid_index);

        let num_tetra = 6.0 * dims.iter().map(|x| *x as f64).product::<f64>();
        let tetra_volume = 1.0 / num_tetra;

        KCache {
            dims,
            k_start,
            k_stop,
            tetra_volume,
            points,
            ks,
        }
    }
}

/// Cache storing the eigenvalues and eigenvectors of a model computed
/// on the k-point grid defined by kcache.
pub struct EvecCache {
    kcache: KCache,
    bands: usize,
    energy: Array2<f64>,
    evec: Array3<Complex64>,
}

impl EvecCache {
    /// Construct a new EvecCache from the given model, with the given number of k-points in
    /// each reciprocal lattice direction and occupying the given region of the Brillouin zone.
    ///
    /// TODO prefer to pass hk_fn as trait giving hk_lat()?
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
    /// Can implement alternate new() function to construct EvecCache this way.
    ///
    /// # Arguments
    ///
    /// * `hk_fn` - a function H(k_lat) giving the Hamiltonian at a particular
    /// k which is given in lattice coordinates.
    ///
    /// * `bands` - the number of bands in H(k).
    ///
    /// * `dims` - the number of points to sample inside the Brillouin zone
    /// along each (lattice coordinate) direction.
    ///
    /// * `k_start` - smallest value of k (in lattice coordinates) to sample
    /// in each direction.
    ///
    /// * `k_start` - largest value of k (in lattice coordinates) to sample in
    /// each direction. To sample the whole Brillouin zone, give
    /// k_start = [0.0, 0.0, 0.0] and k_stop = [1.0, 1.0, 1.0].
    pub fn new<F: Sync + Fn([f64; 3]) -> Array2<Complex64>>(
        hk_fn: F,
        bands: usize,
        dims: [usize; 3],
        k_start: [f64; 3],
        k_stop: [f64; 3],
    ) -> EvecCache {
        let kcache = KCache::new(dims, k_start, k_stop);

        let solutions: Vec<Solution<Complex64, f64>> = kcache
            .ks
            .par_iter()
            .map(|k| {
                let hk = hk_fn(*k);
                // TODO could return Err if this assertion fails instead of panicking.
                assert!(hk.dim() == (bands, bands));

                SymEigen::compute(&hk, Symmetric::Upper, true).unwrap()
            })
            .collect();

        // TODO - could create these as uninitialized: will be visiting all elements.
        // For simplicity, avoid unsafe behavior here.
        let mut energy = Array2::zeros([kcache.points.len(), bands]);
        let mut evec = Array3::zeros([kcache.points.len(), bands, bands]);

        for (k_index, sol) in solutions.into_iter().enumerate() {
            energy.subview_mut(Axis(0), k_index).assign(&sol.values);
            evec.subview_mut(Axis(0), k_index).assign(
                &sol.right_vectors
                    .unwrap(),
            );
        }

        EvecCache {
            kcache,
            bands,
            energy,
            evec,
        }
    }
}

impl KGrid for EvecCache {
    fn dims(&self) -> [usize; 3] {
        self.kcache.dims
    }

    fn points(&self) -> &Vec<[usize; 3]> {
        &self.kcache.points
    }

    fn grid_index(&self, point: &[usize; 3]) -> usize {
        get_grid_index(point, &self.kcache.dims)
    }

    fn ks(&self) -> &Vec<[f64; 3]> {
        &self.kcache.ks
    }

    fn k_range(&self) -> ([f64; 3], [f64; 3]) {
        (self.kcache.k_start, self.kcache.k_stop)
    }

    fn tetra_volume(&self) -> f64 {
        self.kcache.tetra_volume
    }
}

impl EnergyGrid for EvecCache {
    fn energy(&self) -> &Array2<f64> {
        &self.energy
    }

    fn bands(&self) -> usize {
        self.bands
    }
}

impl EvecGrid for EvecCache {
    fn evec(&self) -> &Array3<Complex64> {
        &self.evec
    }
}
