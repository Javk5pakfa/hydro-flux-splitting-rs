pub mod scheme;

use std::{io::Write, fs};
use crate::scheme::*;

#[allow(dead_code)]
enum Setup {
    SquareWave,
    Gaussian,
    Cavitation
}


fn main() {

    let setup = Setup::SquareWave;
    let num_cells = 10000;

    let dx = 1.0 / num_cells as f64;
    let dt = dx * 0.1;

    let tfinal = 0.2;
    
    let x = (0..num_cells).map(|i| -0.5 + (i as f64 + 0.5) * dx);
    let mut t = 0.0;

    let mut u = match setup {
         Setup::SquareWave => {
             let density_left   = 10.0;
             let density_right  = 1.0;
             let pressure_left  = 13.3;
             let pressure_right = 0.1;
             let epsilon_left   = pressure_left  / density_left  / (GAMMA_LAW_INDEX - 1.0);
             let epsilon_right  = pressure_right / density_right / (GAMMA_LAW_INDEX - 1.0);
             let entropy_left   = Primitive(density_left , 0.0, epsilon_left,  0.0).to_conserved().entropy_density();
             let entropy_right  = Primitive(density_right, 0.0, epsilon_right, 0.0).to_conserved().entropy_density();
             x.clone().map(|x| Conserved(if x < 0.0 {density_left} else {density_right}, 0.0, if x < 0.0 {entropy_left} else {entropy_right}, if x < 0.0 {2.0} else {1.0})).collect()
         },
         Setup::Gaussian   => x.clone().map(|x| Conserved(f64::exp((-x * x) / 0.01) + 1.0, 0.0, f64::exp((-x * x) / 0.01) + 1.0, if x < 0.0 {2.0} else {1.0})).collect(),
         Setup::Cavitation => {
             let density   = 1.0;
             let pressure  = 1.0;
             let epsilon   = pressure  / density  / (GAMMA_LAW_INDEX - 1.0);
             let entropy   = Primitive(density, 0.0, epsilon, 0.0).to_conserved().entropy_density();
             x.clone().map(|x| Conserved(density, if x < 0.0 {-1.0} else {0.0}, entropy, if x < 0.0 {2.0} else {1.0})).collect()
         }
    };

    while t < tfinal {
        u = hllc_solver(u, dx, dt);
        t += dt;
        println!("t = {}", t);
    }

    let file = fs::File::create("solution.dat").unwrap();

    for (x, u) in x.zip(u) {
        writeln!(&file, "{:.6} {:.6} {:.6} {:.6} {:.6}", x, u.density(), u.to_primitive().velocity(), u.to_primitive().specific_internal_energy(), u.to_primitive().passive_scalar()).unwrap();
    }
}
