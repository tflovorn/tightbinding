use rulinalg::matrix::Matrix;

struct Scf {
    d: Matrix<f64>,
}

impl Scf {
    pub fn new<P: AsRef<Path>>(scf_path: P) -> Scf {
        let mut d = Matrix::<f64>::zeros(3, 3);

        Scf { d }
    }
}
