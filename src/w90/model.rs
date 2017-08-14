use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io::Read;
use num_complex::Complex64;
use rulinalg::matrix::Matrix;

use model::Model;

pub struct W90Model {
    hr: HashMap<[i32; 3], Matrix<Complex64>>,
    bands: usize,
    d: Matrix<f64>,
}

impl W90Model {
    pub fn new<P: AsRef<Path>>(hr_path: P, d: Matrix<f64>) -> W90Model {
        // TODO should handle this with Result: can happen due to user error.
        let mut f = File::open(hr_path).expect("hr file not found");

        let mut contents = String::new();
        f.read_to_string(&mut contents).expect("error reading hr file");

        let header = extract_hr_header(&contents);

        let hr = extract_hr_model(&contents, &header);

        W90Model { hr, bands: header.bands, d }
    }
}

impl Model for W90Model {
    fn hrs(&self) -> &HashMap<[i32; 3], Matrix<Complex64>> {
        &self.hr
    }

    fn bands(&self) -> usize {
        self.bands
    }

    fn d(&self) -> &Matrix<f64> {
        &self.d
    }
}

struct HrHeader {
    comment_line: String,
    bands: usize,
    rs: usize,
    degen: Vec<u32>,
    start_hr: usize,
}

fn extract_hr_header(contents: &str) -> HrHeader {
    let comment_line = String::new();
    let bands = 0;
    let rs = 0;
    let degen = vec![];
    let start_hr = 0;

    HrHeader { comment_line, bands, rs, degen, start_hr }
}

fn extract_hr_model(contents: &str, header: &HrHeader) -> HashMap<[i32; 3], Matrix<Complex64>> {
    let mut hr = HashMap::new();

    hr
}
