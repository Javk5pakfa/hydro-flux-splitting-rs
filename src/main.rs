pub mod flux_split;

use std::io::Write;
use std::fs;
pub use crate::flux_split::{Conserved, Primitive, flux_split};


fn main() {
    let num_cells = 1000;
    let dx = 1.0 / num_cells as f64;
    let dt = dx * 0.01;
    let x = (0..num_cells).map(|i| -0.5 + (i as f64 + 0.5) * dx);
    let mut t = 0.0;
    let mut u: Vec<_> = x.clone().map(|x| Conserved(f64::exp((-x * x) / 0.01) + 1.0, 0.0)).collect();

    let tfinal = 0.1;

    while t < tfinal {
        u = flux_split(u, dx, dt);
        t += dt;
    }

    let file = fs::File::create("solution.dat").unwrap();

    for (x, u) in x.zip(u) {
        writeln!(&file, "{:.6} {:.6} {:.6}", x, u.get_density(), u.get_momentum()).unwrap();
    }
}
