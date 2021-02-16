pub mod flux_split;

use std::{io::Write, fs};
use crate::flux_split::*;


fn main() {
    let num_cells = 1000;
    let dx = 1.0 / num_cells as f64;
    let dt = dx * 0.1;
    let x = (0..num_cells).map(|i| -0.5 + (i as f64 + 0.5) * dx);
    let mut t = 0.0;
    //let mut u: Vec<_> = x.clone().map(|x| Conserved(f64::exp((-x * x) / 0.01) + 1.0, 0.0)).collect();
    let mut u: Vec<_> = x.clone().map(|x| Conserved(if x<0.0 {1.0} else {0.5}, 0.0)).collect();

    let tfinal = 0.2;

    while t < tfinal {
        u = flux_split(u, dx, dt);
        t += dt;
    }

    let file = fs::File::create("solution.dat").unwrap();

    for (x, u) in x.zip(u) {
        writeln!(&file, "{:.6} {:.6} {:.6}", x, u.density(), u.momentum()).unwrap();
    }
}
