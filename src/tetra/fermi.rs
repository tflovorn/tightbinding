use tetra::grid::EnergyGrid;
use tetra::weight::all_weights;
use tetra::sum::total_number;

pub fn find_fermi<G: EnergyGrid>(grid: &G, occupation: f64) -> f64 {
    let (e_min, e_max) = grid.energy_bounds();

    let occupation_error = |fermi| {
        let weights = all_weights(grid, fermi);

        total_number(&weights) - occupation
    };

    let eps_abs = 1e-12; // occupation has scale ~ 1.
    bisect(occupation_error, e_min, e_max, eps_abs).expect(
        "min and max energy do not bracket Fermi energy",
    )
}

fn bisect<F>(f: F, a: f64, b: f64, eps_abs: f64) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    let f_a = f(a);
    let f_b = f(b);

    // Is root on one of the endpoints?
    if f_a.abs() < eps_abs {
        return Some(a);
    } else if f_b.abs() < eps_abs {
        return Some(b);
    }

    // Do the endpoints bracket a root?
    // If not, no root is present.
    if (f_a < 0.0 && f_b < 0.0) || (f_a > 0.0 && f_b > 0.0) {
        return None;
    }

    // Root is somewhere between the endpoints. Find it by bisection.
    let mut low = a;
    let mut high = b;
    if f_b < 0.0 {
        low = b;
        high = a;
    }

    // TODO - want a maximum iteration number?
    // If f is not continuous, may have no root in bracket. In this case,
    // we will loop forever here.
    loop {
        let mid = (low + high) / 2.0;
        let f_mid = f(mid);

        if f_mid.abs() < eps_abs {
            return Some(mid);
        } else if f_mid > 0.0 {
            // mid is above the root. The root is between low and mid.
            high = mid;
        } else {
            // mid is below the root. The root is between mid and high.
            low = mid;
        }
    }
}
