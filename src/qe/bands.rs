use std::str;
use std::path::Path;
use std::fs::File;
use std::io;
use std::io::Read;

use vec_util::transpose_vecs;

pub struct DftBands {
    ks: Vec<[f64; 3]>,
    Emks: Vec<Vec<f64>>,
}

impl DftBands {
    pub fn new<P: AsRef<Path>>(bands_path: P) -> Result<DftBands, io::Error> {
        let mut contents = String::new();

        {
            let mut f = File::open(bands_path)?;
            f.read_to_string(&mut contents).expect(
                "error reading bands.dat file",
            );
        }

        let (bands, num_ks) = extract_bands_header(&contents);

        println!("bands: {}", bands);
        println!("num_ks: {}", num_ks);

        let (ks, Ekms) = extract_bands(&contents, bands, num_ks);

        let Emks = transpose_vecs(&Ekms);

        Ok(DftBands { ks, Emks })
    }
    /// A list of the k-points in the calculation, given in Cartesian
    /// coordinates in units of Angstrom^{-1}.
    pub fn ks(&self) -> &Vec<[f64; 3]> {
        &self.ks
    }

    /// A list of band energies for each band: Emks()[m][k] gives the
    /// energy of the m'th band at the k-point corresponding to the
    /// k'th entry of ks.
    ///
    /// The number of bands in the calculation is given by the length of Emks().
    /// Each Emks()[m] has the same length, equal to the length of ks().
    pub fn Emks(&self) -> &Vec<Vec<f64>> {
        &self.Emks
    }
}

fn extract_bands_header(contents: &str) -> (usize, usize) {
    let header_line = contents.lines().next().unwrap();
    
    // Header line has the format:
    // 
    // &plot nbnd= 110, nks=   241 /
    //
    // To extract nbnd and nks values, we take the value between the first "=" and ","
    // and between the second "=" and "/".
    
    let bands = header_line.split("=").skip(1).next().unwrap().split(",").next().unwrap().trim().parse().unwrap();

    let num_ks = header_line.split("=").skip(2).next().unwrap().split("/").next().unwrap().trim().parse().unwrap();

    (bands, num_ks)
}

fn extract_bands(contents: &str, bands: usize, num_ks: usize) -> (Vec<[f64; 3]>, Vec<Vec<f64>>) {
    let mut ks = Vec::new();
    let mut Ekms = Vec::new();

    (ks, Ekms)
}
