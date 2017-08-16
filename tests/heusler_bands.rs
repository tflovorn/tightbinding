extern crate tightbinding;

use tightbinding::qe::DftBands;

#[test]
fn heusler_bands() {
    let bands_path = "test_data/NiHfSn/NiHfSn_bulk_soc_bands.dat";

    let bands = DftBands::new(bands_path);
}
