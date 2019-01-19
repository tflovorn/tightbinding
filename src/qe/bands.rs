use std::cmp::min;
use std::f64::consts::PI;
use std::fs::File;
use std::io;
use std::io::Read;
use std::path::Path;
use std::str;

use vec_util::transpose_vecs;

pub struct DftBands {
    ks: Vec<[f64; 3]>,
    emks: Vec<Vec<f64>>,
}

impl DftBands {
    pub fn new<P: AsRef<Path>>(bands_path: P, alat: f64) -> Result<DftBands, io::Error> {
        let mut contents = String::new();

        {
            let mut f = File::open(bands_path)?;
            f.read_to_string(&mut contents)
                .expect("error reading bands.dat file");
        }

        let (num_bands, num_ks) = extract_bands_header(&contents);

        let (ks, ekms) = extract_bands(&contents, num_bands, num_ks, alat);

        let emks = transpose_vecs(&ekms);

        Ok(DftBands { ks, emks })
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
    pub fn emks(&self) -> &Vec<Vec<f64>> {
        &self.emks
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

    let num_bands = header_line
        .split("=")
        .skip(1)
        .next()
        .unwrap()
        .split(",")
        .next()
        .unwrap()
        .trim()
        .parse()
        .unwrap();

    let num_ks = header_line
        .split("=")
        .skip(2)
        .next()
        .unwrap()
        .split("/")
        .next()
        .unwrap()
        .trim()
        .parse()
        .unwrap();

    (num_bands, num_ks)
}

fn extract_bands(
    contents: &str,
    num_bands: usize,
    num_ks: usize,
    alat: f64,
) -> (Vec<[f64; 3]>, Vec<Vec<f64>>) {
    // Skip header line.
    let mut lines = contents.lines().skip(1);

    let mut ks = Vec::new();
    let mut ekms: Vec<Vec<f64>> = Vec::new();

    // Band entries have the format:
    //   kx  ky  kz
    // E0   E1   E2  ... E8   E9
    // E10  E11  E12 ... E18  E19
    // ...
    //
    // kx, ky, kz are in units of 2pi/alat.
    // Energies are in units of eV.
    //
    // kx, ky, kz can be safely split on " ".
    // Energies may take the full space and may not be split on " " - the exact
    // column values must be used.

    let bytes_per_e = 9;
    let es_per_line = 10;

    for _ in 0..num_ks {
        let k_line = lines.next().unwrap().trim();
        let this_ks = k_line
            .split(" ")
            .filter(|d| d.len() > 0)
            .map(|d| (2.0 * PI / alat) * d.parse::<f64>().unwrap())
            .collect::<Vec<f64>>();
        assert_eq!(this_ks.len(), 3);
        ks.push([this_ks[0], this_ks[1], this_ks[2]]);

        let mut ek = Vec::new();

        while ek.len() < num_bands {
            let e_line = lines.next().unwrap();
            for n in 0..min(es_per_line, num_bands - ek.len()) {
                let (e_start, e_stop) = (bytes_per_e * n, bytes_per_e * (n + 1));

                let e = e_line[e_start..e_stop].trim().parse().unwrap();

                ek.push(e);
            }
        }

        ekms.push(ek);
    }

    (ks, ekms)
}
