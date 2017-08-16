use std::str;
use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io;
use std::io::Read;
use num_complex::Complex64;
use ndarray::Array2;

use model::Model;

#[derive(Clone)]
pub struct W90Model {
    hrs: HashMap<[i32; 3], Array2<Complex64>>,
    bands: usize,
    d: Array2<f64>,
}

impl W90Model {
    pub fn new<P: AsRef<Path>>(hr_path: P, d: Array2<f64>) -> Result<W90Model, io::Error> {
        let mut contents = String::new();

        {
            let mut f = File::open(hr_path)?;
            f.read_to_string(&mut contents).expect(
                "error reading hr.dat file",
            );
        }

        let header = extract_hr_header(&contents);

        let hrs = extract_hr_model(&contents, &header);

        Ok(W90Model {
            hrs,
            bands: header.bands,
            d,
        })
    }
}

/// Tight-binding model specified by (R, H(R)) pairs, as extracted from
/// Wannier90 calculation.
impl Model for W90Model {
    /// A collection of (displacement vector, hopping matrix) pairs, (R, H(R)).
    /// H(R) matrix elements are given in units of eV.
    fn hrs(&self) -> &HashMap<[i32; 3], Array2<Complex64>> {
        &self.hrs
    }

    /// The number of bands in the model. Each matrix value in hrs is
    /// nbands x nbands.
    fn bands(&self) -> usize {
        self.bands
    }

    /// A matrix with columns giving the lattice vectors in Cartesian
    /// coordinates in units of Angstroms.
    fn d(&self) -> &Array2<f64> {
        &self.d
    }
}

struct HrHeader {
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
    // Skip comment line.
    let mut lines = contents.lines().skip(1);

    let bands = lines.next().unwrap().trim().parse::<usize>().unwrap();
    let rs = lines.next().unwrap().trim().parse::<usize>().unwrap();

    let mut start_hr = 3;
    let mut degen = vec![];
    while degen.len() < rs {
        let mut degen_line = lines
            .next()
            .unwrap()
            .trim()
            .split(" ")
            .filter(|d| d.len() > 0)
            .map(|d| d.parse::<u32>().unwrap())
            .collect();
        degen.append(&mut degen_line);

        start_hr += 1;
    }

    HrHeader {
        bands,
        rs,
        degen,
        start_hr,
    }
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
/// ip and i are given in the file starting at 1, not 0. We store them here
/// starting at 0.
///
/// # Arguments
/// * `contents` - the contents of the hr.dat file
/// * `header` - information in the header of the hr.dat file, obtained from
/// extract_hr_header().
fn extract_hr_model(contents: &str, header: &HrHeader) -> HashMap<[i32; 3], Array2<Complex64>> {
    let mut lines = contents.lines().skip(header.start_hr);
    let mut hrs = HashMap::new();

    for r_index in 0..header.rs {
        let mut hrs_entry = Array2::<Complex64>::zeros((header.bands, header.bands));
        let mut r = [0, 0, 0];

        for i in 0..header.bands {
            for ip in 0..header.bands {
                let line_contents = lines
                    .next()
                    .unwrap()
                    .trim()
                    .split(" ")
                    .filter(|d| d.len() > 0)
                    .collect::<Vec<&str>>();
                let ra: i32 = line_contents[0].parse().unwrap();
                let rb: i32 = line_contents[1].parse().unwrap();
                let rc: i32 = line_contents[2].parse().unwrap();

                let ip_from_line = line_contents[3].parse::<usize>().unwrap() - 1;
                let i_from_line = line_contents[4].parse::<usize>().unwrap() - 1;

                assert_eq!(ip, ip_from_line);
                assert_eq!(i, i_from_line);

                if ip == 0 && i == 0 {
                    r = [ra, rb, rc];
                } else {
                    assert_eq!(r, [ra, rb, rc]);
                }

                let val_re: f64 = line_contents[5].parse().unwrap();
                let val_im: f64 = line_contents[6].parse().unwrap();

                let degen = header.degen[r_index] as f64;

                hrs_entry[[ip, i]] = Complex64::new(val_re / degen, val_im / degen);
            }

        }

        hrs.insert(r, hrs_entry);
    }

    hrs
}
