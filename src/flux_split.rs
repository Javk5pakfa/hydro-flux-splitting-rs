



/**
 * Conserved hydrodynamic quantities: density, momentum
 */
#[derive(Copy, Clone)]
pub struct Conserved(pub f64, pub f64);




/**
 * Primitive hydrodynamic quantities: density, velocity
 */
#[derive(Copy, Clone)]
pub struct Primitive(pub f64, pub f64);




// ============================================================================
pub static KAPPA: f64 = 1.0;
pub static GAMMA_LAW_INDEX: f64 = 5.0 / 3.0;




// ============================================================================
impl Conserved {

    pub fn density(self) -> f64 {
        self.0
    }

    pub fn momentum(self) -> f64 {
        self.1
    }

    pub fn velocity(self) -> f64 {
        self.momentum() / self.density()
    }

    pub fn to_primitive(self) -> Primitive {
        Primitive(self.density(), self.velocity())
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

    pub fn to_conserved(self) -> Conserved {
        Conserved(self.density(), self.density() * self.velocity())
    }

    pub fn get_fluxes(self) -> Conserved {
        let momentum = self.density() * self.velocity();
        let energy = momentum.powi(2) + KAPPA * self.density().powf(GAMMA_LAW_INDEX);
        Conserved(momentum, energy)
    }

    pub fn sound_speed(self) -> f64 {
        (KAPPA * GAMMA_LAW_INDEX * self.density().powf(GAMMA_LAW_INDEX - 1.0)).sqrt()
    }

    pub fn eigenval_plus(self) -> f64 {
        self.velocity() + self.sound_speed()
    }

    pub fn eigenval_minus(self) -> f64 {
        self.velocity() - self.sound_speed()
    }

    pub fn max_eigenval(self) -> f64 {
        self.eigenval_plus().abs().max(self.eigenval_minus().abs())
    }
}




// NOTE: the implementations below can be auto-generated with the derive_more
// crate.

// ============================================================================
impl std::ops::Add<Conserved> for Conserved { type Output = Self; fn add(self, u: Conserved) -> Conserved { Conserved(self.0 + u.0, self.1 + u.1) } }
impl std::ops::Sub<Conserved> for Conserved { type Output = Self; fn sub(self, u: Conserved) -> Conserved { Conserved(self.0 - u.0, self.1 - u.1) } }
impl std::ops::Mul<Conserved> for f64 { type Output = Conserved; fn mul(self, u: Conserved) -> Conserved { Conserved(self * u.0, self * u.1) } }
impl std::ops::Mul<f64> for Conserved { type Output = Conserved; fn mul(self, a: f64) -> Conserved { Conserved(self.0 * a, self.1 * a) } }
impl std::ops::Div<f64> for Conserved { type Output = Conserved; fn div(self, a: f64) -> Conserved { Conserved(self.0 / a, self.1 / a) } }




// ============================================================================
pub fn flux_split(u: Vec<Conserved>, dx: f64, dt: f64) -> Vec<Conserved> {
    let n = u.len();
    let max_lambda: f64 = u.iter().map(|u| u.to_primitive().max_eigenval()).fold(0.0, |a, b| a.max(b));

    if dx / dt < max_lambda {
        panic!("dx/dt too small!")
    }

    let mut u1 = vec![Conserved(0.0, 0.0); n];

    // Update interior cells.
    for i in 1..n-1 {
        let prim_right = u[i+1].to_primitive();
        let prim_left = u[i-1].to_primitive();
        let fr = prim_right.get_fluxes();
        let fl = prim_left.get_fluxes();
        let flux_cor_term = max_lambda * (u[i-1] - 2.0 * u[i] + u[i+1]);

        u1[i] = u[i] - 0.5 * (dt / dx) * (fr - fl - flux_cor_term);
    }

    // Update boundaries.
    u1[0] = u[0];
    u1[n-1] = u[n-1];
    u1
}
