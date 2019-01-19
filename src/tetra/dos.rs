use tetra::grid::EnergyGrid;

pub fn dos_contrib<G: EnergyGrid>(grid: &G, e: f64, sorted_es: &[f64; 4]) -> f64 {
    let fac = grid.tetra_volume();
    let (e1, e2, e3, e4) = (sorted_es[0], sorted_es[1], sorted_es[2], sorted_es[3]);

    if e <= e1 {
        return 0.0;
    } else if e1 <= e && e <= e2 {
        return fac * 3.0 * (e - e1).powi(2) / ((e2 - e1) * (e3 - e1) * (e4 - e1));
    } else if e2 <= e && e <= e3 {
        return fac
            * (1.0 / ((e3 - e1) * (e4 - e1)))
            * (3.0 * (e2 - e1) + 6.0 * (e - e2)
                - 3.0 * ((e3 - e1) + (e4 - e2)) * (e - e2).powi(2) / ((e3 - e2) * (e4 - e2)));
    } else if e3 <= e && e <= e4 {
        return fac * 3.0 * (e4 - e).powi(2) / ((e4 - e1) * (e4 - e2) * (e4 - e3));
    }
    // e > e4
    0.0
}
