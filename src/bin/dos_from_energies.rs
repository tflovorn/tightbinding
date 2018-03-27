extern crate clap;
extern crate serde_json;
extern crate tightbinding;

use std::io::Read;
use std::path::PathBuf;
use std::fs::File;
use clap::{App, Arg};
use tightbinding::tetra::EnergyCache;
use tightbinding::dos::total_dos_from_num;

fn main() {
    let args = App::new("Plot DOS from sampled energies on Brillouin zone grid")
        .version("0.1.0")
        .arg(
            Arg::with_name("energy_file_path")
                .long("energy_file_path")
                .takes_value(true)
                .help(
                    "Path to a JSON file containing a list of energies Es[k0][k1][k2][m], \
                     where (k0, k1, k2) are indices of k-points that tile the Brillouin zone \
                     in reciprocal lattice coordinates (in the range [0, 1) along each \
                     direction) and m is the band index.",
                ),
        )
        .arg(
            Arg::with_name("num_dos_energies")
                .long("num_dos_energies")
                .takes_value(true)
                .default_value("1000")
                .help("Number of energy values to include in the DOS calculation."),
        )
        .arg(
            Arg::with_name("min_energy")
                .long("min_energy")
                .takes_value(true)
                .help(
                    "Minimum energy value to include in the calculation. \
                     If either min_energy or max_energy is specified, must specify the other.",
                ),
        )
        .arg(
            Arg::with_name("max_energy")
                .long("max_energy")
                .takes_value(true)
                .help(
                    "Maximum energy value to include in the calculation. \
                     If either min_energy or max_energy is specified, must specify the other.",
                ),
        )
        .get_matches();

    let energy_file_path = PathBuf::from(args.value_of("energy_file_path").unwrap());
    let num_dos_energies = args.value_of("num_dos_energies").unwrap().parse().unwrap();

    let (min_e, max_e) = (args.value_of("min_energy"), args.value_of("max_energy"));

    // If min_energy and max_energy are both specified, energy_bounds = (min_e, max_e).
    // Otherwise, energy_bounds = None.
    let energy_bounds = match min_e {
        Some(min_e) => match max_e {
            Some(max_e) => Some((min_e.parse().unwrap(), max_e.parse().unwrap())),
            None => None,
        },
        None => None,
    };

    let mut in_file = File::open(&energy_file_path).expect("Error opening input file");
    let mut contents = String::new();
    in_file.read_to_string(&mut contents).unwrap();
    let energies: Vec<Vec<Vec<Vec<f64>>>> =
        serde_json::from_str(&contents).expect("Error reading input file");

    let cache = EnergyCache::new(&energies);

    let total_dos = total_dos_from_num(&cache, num_dos_energies, energy_bounds);

    let energy_file_stem = energy_file_path.file_stem().unwrap().to_os_string();
    let mut output_file_stem = energy_file_stem.clone();
    output_file_stem.push("_total_dos.json");
    let output_file_path = energy_file_path.with_file_name(output_file_stem);

    let out_file = File::create(output_file_path).unwrap();
    serde_json::to_writer(&out_file, &total_dos).expect("Error writing output file");
}
