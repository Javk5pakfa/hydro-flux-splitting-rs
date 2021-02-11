mod burgers;

use std::io::Write;
use std::fs;

fn main() {
    let num_cells = 1000;
    let dx = 1.0 / num_cells as f64;
    let dt = dx * 0.5;
    let x = (0..num_cells).map(|i| -0.5 + (i as f64 + 0.5) * dx);
    let mut t = 0.0;
    let mut u: Vec<_> = x.clone().map(|x| f64::exp(-(x * x) / 0.01)).collect();

    while t < 0.1 {
        u = burgers::next(u, dx, dt);
        t += dt;
    }

    let file = fs::File::create("solution.dat").unwrap();

    for (x, u) in x.zip(u) {
        writeln!(&file, "{:.6} {:.6}", x, u).unwrap();
    }
}
