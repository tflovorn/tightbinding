use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io;
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
    pub fn new<P: AsRef<Path>>(hr_path: P, d: Matrix<f64>) -> Result<W90Model, io::Error> {
        let mut contents = String::new();

        {
            let mut f = File::open(hr_path)?;
            f.read_to_string(&mut contents).expect("error reading hr file");
        }

        let header = extract_hr_header(&contents);

        let hr = extract_hr_model(&contents, &header);

        Ok(W90Model { hr, bands: header.bands, d })
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

/// Wannier90 hr.dat file header has the format:
/// comment line
/// number of bands
/// number of displacement vectors (rs)
/// list of degen values, 15 per line, total number equal to rs
fn extract_hr_header(contents: &str) -> HrHeader {
    let mut lines = contents.lines();

    let comment_line = lines.next().unwrap().to_string();
    let bands = lines.next().unwrap().trim().parse::<usize>().unwrap();
    let rs = lines.next().unwrap().trim().parse::<usize>().unwrap();

    let mut start_hr = 3;
    let mut degen = vec![];
    while degen.len() < rs {
        let mut degen_line = lines.next().unwrap().trim().split(" ").map(|d| d.parse::<u32>().unwrap()).collect();
        degen.append(&mut degen_line);

        start_hr += 1;
    }

    HrHeader { comment_line, bands, rs, degen, start_hr }
}

fn extract_hr_model(contents: &str, header: &HrHeader) -> HashMap<[i32; 3], Matrix<Complex64>> {
    let mut hr = HashMap::new();

    hr
}
