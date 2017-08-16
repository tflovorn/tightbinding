/// Given a vector of vectors with elements `x[i][j]`, return the corresponding
/// vector of vectors `y` where `y[j][i] == x[i][j]`. 
///
/// Each element of x has clone() called on it to construct y.
///
/// Panics if the lengths of each `x[i]` are not all equal.
pub fn transpose_vecs<T: Clone>(x: &Vec<Vec<T>>) -> Vec<Vec<T>> {
    let mut y = Vec::new();

    let ni = x.len();

    if ni == 0 {
        return y;
    }

    for i in 1..ni {
        assert!(x[i].len() == x[0].len());
    }

    let nj = x[0].len();

    for j in 0..nj {
        y.push(Vec::new());
        for i in 0..ni {
            y[j].push(x[i][j].clone());
        }
    }

    y
}
