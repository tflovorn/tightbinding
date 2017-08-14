use std::str;
use std::fmt;
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

/// Parse the header lines from the hr.dat file given by contents.
///
/// The Wannier90 hr.dat file header has the format:
///
/// comment line
/// number of bands
/// number of displacement vectors (rs)
/// list of degen values, 15 per line, total number equal to rs
///
/// # Arguments
///
/// * `contents` - the contents of the hr.dat file
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

/// Get the next value from the Split xs, parsing it as type F.
fn next_value<F>(xs: &mut str::Split<&str>) -> F
where
    F: str::FromStr,
    F::Err : fmt::Debug
{
    xs.next().unwrap().parse().unwrap()
}

/// Parse the tight-binding Hamiltonian from the hr.dat file given by contents.
///
/// The tight-binding Hamiltonian is specified in the Wannier90 hr.dat file
/// one matrix element per line, with each line having the format:
///
/// Ra  Rb  Rc  ip  i  Re{[H(R)]_{ip, i}}*degen[R]  Im{[H(R)]_{ip, i}}*degen[R]
///
/// ip and i are tight-binding basis indices, and R is the displacement vector.
/// ip varies the fastest, then i, then R.
///
/// # Arguments
/// * `contents` - the contents of the hr.dat file
/// * `header` - information in the header of the hr.dat file, obtained from
/// extract_hr_header().
fn extract_hr_model(contents: &str, header: &HrHeader) -> HashMap<[i32; 3], Matrix<Complex64>> {
    let mut lines = contents.lines().skip(header.start_hr);
    let mut hr = HashMap::new();

    let mut r_index = 0;
    let mut i = 0;
    let mut ip = 0;

    while r_index < header.rs {
        let mut hr_entry = Matrix::<Complex64>::zeros(header.bands, header.bands);
        let mut r = [0, 0, 0];

        while i < header.bands {
            while ip < header.bands {
                let mut line_contents = lines.next().unwrap().trim().split(" ");
                let ra: i32 = next_value(&mut line_contents);
                let rb: i32 = next_value(&mut line_contents);
                let rc: i32 = next_value(&mut line_contents);

                let ip_from_line = next_value::<usize>(&mut line_contents) - 1;
                let i_from_line = next_value::<usize>(&mut line_contents) - 1;

                assert_eq!(ip, ip_from_line);
                assert_eq!(i, i_from_line);

                if ip == 0 && i == 0 {
                    r = [ra, rb, rc];
                } else {
                    assert_eq!(r, [ra, rb, rc]);
                }

                let val_re: f64 = next_value(&mut line_contents);
                let val_im: f64 = next_value(&mut line_contents);

                let degen = header.degen[r_index] as f64;

                hr_entry[[ip, i]] = Complex64::new(val_re / degen, val_im / degen);
                ip += 1;
            }

            i += 1;
        }

        hr.insert(r, hr_entry);
        r_index += 1;
    }

    hr
}
