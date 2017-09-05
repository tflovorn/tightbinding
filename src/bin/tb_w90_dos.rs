extern crate clap;
extern crate serde_json;
extern crate tightbinding;

use std::fs::File;
use clap::{Arg, App};
use tightbinding::paths::{get_work_base, build_work, get_scf_path, get_hr_path};
use tightbinding::Model;
use tightbinding::w90::W90Model;
use tightbinding::qe::Scf;
use tightbinding::fourier::hk_lat;
use tightbinding::tetra::EvecCache;
use tightbinding::dos::dos_from_num;

fn main() {
    let args = App::new("Plot DOS from Wannier90 tight-binding model")
        .version("0.1.0")
        .arg(Arg::with_name("subdir").long("subdir").takes_value(true))
        .arg(
            Arg::with_name("num_energies")
                .long("num_energies")
                .takes_value(true)
                .default_value("1000"),
        )
        .arg(Arg::with_name("prefix").index(1).required(true))
        .get_matches();

    let prefix = args.value_of("prefix").unwrap();
    let num_energies = args.value_of("num_energies").unwrap().parse().unwrap();

    let work_base = get_work_base().unwrap();
    let work = build_work(&work_base, args.value_of("subdir"), prefix);
    let scf_path = get_scf_path(&work, prefix);
    let hr_path = get_hr_path(&work, prefix);

    let scf_data = Scf::new(scf_path).expect("could not find scf output file");
    let model = W90Model::new(hr_path, scf_data.d).expect("error constructing model");

    // TODO make dims a parameter.
    // Want to allow specifying k_start, k_stop?
    let dims = [8, 8, 8];
    let k_start = [0.0, 0.0, 0.0];
    let k_stop = [1.0, 1.0, 1.0];

    let hk_fn = |k| hk_lat(&model, &k);

    let cache = EvecCache::new(hk_fn, model.bands(), dims, k_start, k_stop);
    let dos = dos_from_num(&cache, num_energies);

    let out_path = format!("{}_dos.json", prefix);

    let out_file = File::create(out_path).expect("Eror creating output file");
    serde_json::to_writer(&out_file, &dos).expect("Error writing output file");
}
