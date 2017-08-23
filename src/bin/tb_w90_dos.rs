extern crate clap;
extern crate serde_json;
extern crate tightbinding;

use std::env;
use std::path::PathBuf;
use std::fs::File;
use clap::{Arg, App};
use tightbinding::Model;
use tightbinding::w90::W90Model;
use tightbinding::qe::Scf;
use tightbinding::fourier::hk_lat;
use tightbinding::tetra::EvecCache;
use tightbinding::dos::dos_from_num;

fn build_work(work_base: &str, subdir: Option<&str>, prefix: &str) -> PathBuf {
    let mut work = PathBuf::new();
    work.push(work_base);

    if let Some(ref subdir) = subdir {
        work.push(subdir);
    }

    work.push(prefix);

    work
}

fn get_scf_path(work: &PathBuf, prefix: &str) -> PathBuf {
    let mut scf_path = work.clone();
    scf_path.push("scf");
    scf_path.push(format!("{}.save", prefix));
    scf_path.push("data-file.xml");

    scf_path
}

fn get_hr_path(work: &PathBuf, prefix: &str) -> PathBuf {
    let mut hr_path = work.clone();
    hr_path.push("wannier");
    hr_path.push(format!("{}_hr.dat", prefix));

    hr_path
}

fn main() {
    let work_err_msg = "Could not find environment variable TB_W90_WORK,\
                        which should be the data root directory.";
    let work_base = env::var("TB_W90_WORK").expect(work_err_msg);

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

    let work = build_work(&work_base, args.value_of("subdir"), prefix);
    let scf_path = get_scf_path(&work, prefix);
    let hr_path = get_hr_path(&work, prefix);

    let scf_data = Scf::new(scf_path).expect("could not find scf output file");
    let model = W90Model::new(hr_path, scf_data.d).expect("error constructing model");

    // TODO make dims a parameter.
    // Want to allow specifying k_start, k_stop?
    let use_curvature_correction = true;
    let dims = [8, 8, 8];
    let k_start = [0.0, 0.0, 0.0];
    let k_stop = [1.0, 1.0, 1.0];

    let hk_fn = |k| hk_lat(&model, &k);

    let cache = EvecCache::new(hk_fn, model.bands(), dims, k_start, k_stop);
    let dos = dos_from_num(&cache, num_energies, use_curvature_correction);

    let out_path = format!("{}_dos.json", prefix);

    let out_file = File::create(out_path).expect("Eror creating output file");
    serde_json::to_writer(&out_file, &dos).expect("Error writing output file");
}
