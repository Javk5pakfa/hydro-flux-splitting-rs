// Burgers' equations solver.
// Author: Dr. Jonathan Zrake 2021.

pub(crate) fn next(u: Vec<f64>, dx: f64, dt: f64) -> Vec<f64> {
    let n = u.len();
    let mut u1 = vec![0.0; n];

    for i in 0..n {

        // neighboring solution values for periodic BC
        let ul = u[(i + n - 1) % n];
        let ur = u[i];

        // upwinded interface fluxes (assume right-going characteristics)
        let fimh = 0.5 * ul * ul;
        let fiph = 0.5 * ur * ur;

        u1[i] = u[i] - (fiph - fimh) * dt / dx
    }
    u1
}

