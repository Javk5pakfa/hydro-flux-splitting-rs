pub static GAMMA_LAW_INDEX: f64 = 5.0 / 3.0;
use std::{io::Write, fs};


/**
 * Conserved hydrodynamic quantities: density, momentum, entropy density
 */
#[derive(Copy, Clone)]
pub struct Conserved(pub f64, pub f64, pub f64);




/**
 * Primitive hydrodynamic quantities: density, velocity, specific internal energy
 */
#[derive(Copy, Clone)]
pub struct Primitive(pub f64, pub f64, pub f64);



// ============================================================================
impl Conserved {

    pub fn density(self) -> f64 {
        self.0
    }

    pub fn momentum(self) -> f64 {
        self.1
    }

    pub fn energy_density(self) -> f64 {
        self.2
    }

    pub fn velocity(self) -> f64 {
        self.momentum() / self.density()
    }

    pub fn specific_internal_energy(self) -> f64 {
        (self.energy_density() - (0.5) * self.density() * self.velocity().powf(2.0)) / self.density()
    }

    pub fn pressure(self) -> f64 {
        self.density() * self.specific_internal_energy() * ( GAMMA_LAW_INDEX - 1.0 )
    }

    pub fn to_primitive(self) -> Primitive {
        Primitive(self.density(), self.velocity(), self.pressure())
    }
}




// ============================================================================
impl Primitive {

    pub fn density(self) -> f64 {
        self.0
    }

    pub fn velocity(self) -> f64 {
        self.1
    }

    pub fn pressure(self) -> f64 {
        self.2
    }

    pub fn specific_internal_energy(self) -> f64 {
        self.pressure() / self.density() / ( GAMMA_LAW_INDEX - 1.0 )
    }

    pub fn to_conserved(self) -> Conserved {
        Conserved(self.density(), self.density() * self.velocity(), self.density() * self.specific_internal_energy() + 0.5 * self.density() * self.velocity().powf(2.0))
    }

    pub fn get_fluxes(self) -> Conserved {
        let mass_flux     = self.density() * self.velocity();
        let momentum_flux = self.density() * self.velocity().powi(2) + self.pressure();
        let energy_flux  = (self.density() * self.specific_internal_energy() + 0.5 * self.density() * self.velocity().powf(2.0) + self.pressure()) * self.velocity();
        Conserved(mass_flux, momentum_flux, energy_flux)
    }

    pub fn sound_speed(self) -> f64 {
        ( GAMMA_LAW_INDEX * self.pressure() / self.density() ).sqrt()
    }

    pub fn eigenval_plus(self) -> f64 {
        self.velocity() + self.sound_speed()
    }

    pub fn eigenval_0(self) -> f64 {
        self.velocity()
    }

    pub fn eigenval_minus(self) -> f64 {
        self.velocity() - self.sound_speed()
    }

}




// ============================================================================
// Eqn 4.5 of Toro (2009).
pub fn f_pstar(pl: Primitive, pr: Primitive, pres: f64) -> f64 {

    faux(pl,pres) + faux(pr, pres) + pr.velocity() - pl.velocity()

}



// Eqn 4.6/4.7 of Toro (2009).
fn faux(prim: Primitive, pres_guess: f64) -> f64 {
    // Eqn 4.8 of Toro (2009).
    let ak = 2.0 / ( ( GAMMA_LAW_INDEX + 1.0 ) * prim.density() );
    let bk = ( prim.pressure() ) * ( GAMMA_LAW_INDEX - 1.0 ) / ( GAMMA_LAW_INDEX + 1.0 );

    if pres_guess > prim.pressure() {
        // Shock.
        ( pres_guess - prim.pressure() ) * ( ak / ( pres_guess + bk ) ).powf(0.5)
    } else {
        // Rarefaction.
        ( 2.0 * ak ) / ( GAMMA_LAW_INDEX - 1.0 ) * ( ( pres_guess / prim.pressure() ).powf((GAMMA_LAW_INDEX - 1.0) / (2.0 * GAMMA_LAW_INDEX)) - 1.0 )
    }
}



// Eqn 4.9 of Toro (2009).
pub fn u_star(pl: Primitive, pr: Primitive, pstar: f64) -> f64 {
    0.5 * ( pl.velocity() + pr.velocity() ) + 0.5 * ( faux(pr, pstar) - faux(pl, pstar) )
}



// Eqn 4.53 of Toro (2009).
pub fn rho_star_left(pl: Primitive, pstar: f64) -> f64 {
    pl.density() * ( pstar / pl.pressure() ).powf(1.0 / GAMMA_LAW_INDEX)
}



// Eqn 4.57 of Toro (2009).
pub fn rho_star_right(pr: Primitive, pstar: f64) -> f64 {
    pr.density() * ( ( pstar / pr.pressure() + ( ( GAMMA_LAW_INDEX - 1.0 ) / ( GAMMA_LAW_INDEX + 1.0 ) ) ) / ( ( ( GAMMA_LAW_INDEX - 1.0 ) / ( GAMMA_LAW_INDEX + 1.0 ) ) * pstar / pr.pressure() + 1.0 ) )
}



// ============================================================================


pub fn shockwave_speed(pr: Primitive, pstar_r: Primitive) -> f64 {
    ( pr.get_fluxes().0 - pstar_r.get_fluxes().0 ) / ( pr.density() - pstar_r.density() )
}


// ============================================================================



fn rarefaction_state(pl: Primitive, x: f64, t: f64) -> Primitive {
    let r_density = pl.density() * ( 2.0 / (GAMMA_LAW_INDEX - 1.0) + (GAMMA_LAW_INDEX - 1.0) / ((GAMMA_LAW_INDEX + 1.0) * pl.sound_speed()) * (pl.velocity() - x/t) ).powf(2.0 / (GAMMA_LAW_INDEX - 1.0));
    let r_velocity = 2.0 / (GAMMA_LAW_INDEX + 1.0) * ( pl.sound_speed() + (GAMMA_LAW_INDEX - 1.0) / 2.0 * pl.velocity() + x/t );
    let r_pressure = pl.pressure() * ( 2.0 / (GAMMA_LAW_INDEX + 1.0) + (GAMMA_LAW_INDEX - 1.0) / ((GAMMA_LAW_INDEX + 1.0) * pl.sound_speed()) * ( pl.velocity() - x/t ) ).powf(2.0 * GAMMA_LAW_INDEX / ( GAMMA_LAW_INDEX - 1.0 ));

    Primitive(r_density, r_velocity, r_pressure)
}



// ============================================================================



fn find_bracket(pl: Primitive, pr: Primitive) -> (f64, f64) {
    let mut plow = 1e-16;

    if f_pstar(pl, pr, plow) > 0.0 {
        panic!("Whoops! plow is not a lower bracket!")
    }

    let mut phigh = 10.0 * plow;

    while f_pstar(pl, pr, phigh) < 0.0 {
        plow = phigh;
        phigh = 10.0 * plow;
    }

    (plow, phigh)
}



// Root-finder algorithm.
fn find_root(pl: Primitive, pr: Primitive) -> f64 {
    let (mut plow, mut phigh) = find_bracket(pl, pr);

    // Percentage.
    let tol = 1e-6;

    let mut rel_err = 1.0;
    while rel_err > tol {

        let pnew = 0.5 * ( plow + phigh );

        if f_pstar(pl, pr, pnew) > 0.0 {
            phigh = pnew;
        } else if f_pstar(pl, pr, pnew) < 0.0 {
            plow = pnew;
        } else {
            phigh = pnew;
            plow = pnew;
        }

        rel_err = (phigh - plow) / (0.5 * ( phigh + plow ));
    }

    0.5 * ( phigh + plow )
}



fn main() {
    let num_cells = 10000;

    let tfinal = 0.2;
    let dx = 1.0 / num_cells as f64;

    let x = (0..num_cells).map(|i| -0.5 + (i as f64 + 0.5) * dx);

    let density_left   = 10.0;
    let density_right  = 1.0;
    let pressure_left  = 13.3;
    let pressure_right = 0.1;

    let pl = Primitive(density_left, 0.0, pressure_left);
    let pr = Primitive(density_right, 0.0, pressure_right);

    let pstar = find_root(pl, pr);
    let ustar = u_star(pl, pr, pstar);
    let rho_str_lft = rho_star_left(pl, pstar);
    let rho_str_rght = rho_star_right(pr, pstar);

    let us = shockwave_speed(pr, Primitive(rho_str_rght, ustar, pstar));
    let vrp = pl.eigenval_minus();
    let vrb = Primitive(rho_str_lft, ustar, pstar).eigenval_minus();

    let mut rho_sol_vec: Vec<f64> = Vec::new();
    let mut vel_sol_vec: Vec<f64> = Vec::new();
    let mut pres_sol_vec: Vec<f64> = Vec::new();

    for i in x {
        // Find region of i.

        // Region 1.
        if i < vrp * tfinal {
            rho_sol_vec.push(pl.density());
            vel_sol_vec.push(pl.velocity());
            pres_sol_vec.push(pl.pressure());
        } else if i < vrb * tfinal { // Region 2.
            let ptemp = rarefaction_state(pl, i, tfinal);
            rho_sol_vec.push(ptemp.density());
            vel_sol_vec.push(ptemp.velocity());
            pres_sol_vec.push(ptemp.pressure());
        } else if i < ustar * tfinal { // Region 3.
            rho_sol_vec.push(rho_str_lft);
            vel_sol_vec.push(ustar);
            pres_sol_vec.push(pstar);
        } else if i < us * tfinal { // Region 4.
            rho_sol_vec.push(rho_str_rght);
            vel_sol_vec.push(ustar);
            pres_sol_vec.push(pstar);
        } else { // Region 5.
            rho_sol_vec.push(pr.density());
            vel_sol_vec.push(pr.velocity());
            pres_sol_vec.push(pr.pressure());
        }
    }

    let file = fs::File::create("solution.dat").unwrap();

    for i in 0..num_cells {
        writeln!(&file, "{:.6} {:.6} {:.6} {:.6}", -0.5 + (i as f64 + 0.5) * dx, rho_sol_vec[i], vel_sol_vec[i], pres_sol_vec[i]).unwrap();
        // println!("{}", i);
    }
}

