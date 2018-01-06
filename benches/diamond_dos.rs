#[macro_use]
extern crate bencher;
extern crate tightbinding;

use bencher::Bencher;
use tightbinding::Model;
use tightbinding::w90::W90Model;
use tightbinding::qe::Scf;
use tightbinding::fourier::hk_lat;
use tightbinding::tetra::EvecCache;
use tightbinding::dos::dos_from_num;

fn diamond_dos(bench: &mut Bencher) {
    let cache = diamond_cache();

    let num_energies = 100;
    let energy_bounds = None;

    bench.iter(|| {
        dos_from_num(&cache, num_energies, energy_bounds);
    });
}

fn diamond_cache() -> EvecCache {
    let scf_path = "test_data/diamond/scf.data-file.xml";
    let hr_path = "test_data/diamond/diamond_hr.dat";

    let scf_data = Scf::new(scf_path).unwrap();
    let d = scf_data.d;

    let m = W90Model::new(hr_path, d).unwrap();

    let dims = [8, 8, 8];
    let k_start = [0.0, 0.0, 0.0];
    let k_stop = [1.0, 1.0, 1.0];
    let hk_fn = |k| hk_lat(&m, &k);

    EvecCache::new(hk_fn, m.bands(), dims, k_start, k_stop)
}

benchmark_group!(benches, diamond_dos);
benchmark_main!(benches);
