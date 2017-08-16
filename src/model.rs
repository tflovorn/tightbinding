use std::collections::HashMap;
use num_complex::Complex64;
use ndarray::Array2;

/// Tight-binding model specified by (R, H(R)) pairs.
pub trait Model {
    /// A collection of (displacement vector, hopping matrix) pairs, (R, H(R)).
    fn hrs(&self) -> &HashMap<[i32; 3], Array2<Complex64>>;

    /// The number of bands in the model. Each matrix value in hrs is
    /// nbands x nbands.
    fn bands(&self) -> usize;

    /// A matrix with columns giving the lattice vectors in Cartesian
    /// coordinates.
    fn d(&self) -> &Array2<f64>;
}
