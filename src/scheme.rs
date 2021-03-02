



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
pub static GAMMA_LAW_INDEX: f64 = 5.0 / 3.0;




// ============================================================================
impl Conserved {

    pub fn density(self) -> f64 {
        self.0
    }

    pub fn momentum(self) -> f64 {
        self.1
    }

    pub fn entropy_density(self) -> f64 {
        self.2
    }

    pub fn velocity(self) -> f64 {
        self.momentum() / self.density()
    }

    pub fn specific_internal_energy(self) -> f64 {
        self.2 / self.0.powf(2.0 - GAMMA_LAW_INDEX) 
    }

    pub fn to_primitive(self) -> Primitive {
        Primitive(self.density(), self.velocity(), self.specific_internal_energy())
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

    pub fn specific_internal_energy(self) -> f64 {
        self.2
    }

    pub fn to_conserved(self) -> Conserved {
        Conserved(self.density(), self.density() * self.velocity(), self.density() * self.specific_internal_energy() * self.density().powf(1.0 - GAMMA_LAW_INDEX) )
    }

    pub fn get_fluxes(self) -> Conserved {
        let mass_flux     = self.density() * self.velocity();
        let momentum_flux = self.density() * self.velocity().powi(2) + self.density() * self.specific_internal_energy() * (GAMMA_LAW_INDEX - 1.0);
        let entropy_flux  = self.to_conserved().entropy_density() * self.velocity();
        Conserved(mass_flux, momentum_flux, entropy_flux)
    }

    pub fn sound_speed(self) -> f64 {
        let pressure = self.density() * self.specific_internal_energy() * (GAMMA_LAW_INDEX - 1.0);
        ( GAMMA_LAW_INDEX * pressure / self.density()).sqrt()
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

    pub fn max_eigenval(self) -> f64 {
        self.eigenval_plus().abs().max(self.eigenval_minus().abs()).max(self.eigenval_0().abs())
    }

    pub fn max_eigenval_signed(self) -> f64 {
        self.eigenval_plus().max(self.eigenval_minus()).max(self.eigenval_0())

    }

    pub fn min_eigenval_signed(self) -> f64 {
        self.eigenval_plus().min(self.eigenval_minus()).min(self.eigenval_0())

    }
}




// NOTE: the implementations below can be auto-generated with the derive_more
// crate.

// ============================================================================
impl std::ops::Add<Conserved> for Conserved { type Output = Self; fn add(self, u: Conserved) -> Conserved { Conserved(self.0 + u.0, self.1 + u.1, self.2 + u.2) } }
impl std::ops::Sub<Conserved> for Conserved { type Output = Self; fn sub(self, u: Conserved) -> Conserved { Conserved(self.0 - u.0, self.1 - u.1, self.2 - u.2) } }
impl std::ops::Mul<Conserved> for f64 { type Output = Conserved; fn mul(self, u: Conserved) -> Conserved { Conserved(self * u.0, self * u.1, self * u.2) } }
impl std::ops::Mul<f64> for Conserved { type Output = Conserved; fn mul(self, a: f64) -> Conserved { Conserved(self.0 * a, self.1 * a, self.2 * a) } }
impl std::ops::Div<f64> for Conserved { type Output = Conserved; fn div(self, a: f64) -> Conserved { Conserved(self.0 / a, self.1 / a, self.2 / a) } }




// ============================================================================
pub fn flux_split(u: Vec<Conserved>, dx: f64, dt: f64) -> Vec<Conserved> {
    let n = u.len();
    let max_lambda: f64 = u.iter().map(|u| u.to_primitive().max_eigenval()).fold(0.0, |a, b| a.max(b));

    if dx / dt < max_lambda {
        panic!("dx/dt too small!")
    }

    let mut u1 = vec![Conserved(0.0, 0.0, 0.0); n];

    // Update interior cells.
    for i in 1..n-1 {
        let prim_right = u[i+1].to_primitive();
        let prim_left  = u[i-1].to_primitive();
        let fr = prim_right.get_fluxes();
        let fl = prim_left.get_fluxes();
        let flux_cor_term = max_lambda * (u[i-1] - 2.0 * u[i] + u[i+1]);

        u1[i] = u[i] - 0.5 * (dt / dx) * (fr - fl - flux_cor_term);
    }

    // Update boundaries.
    u1[0]   = u[0];
    u1[n-1] = u[n-1];
    u1
}

// ============================================================================
pub fn hll_solver(u: Vec<Conserved>, dx: f64, dt: f64) -> Vec<Conserved> {
    let n = u.len();
    let max_lambda: f64 = u.iter().map(|u| u.to_primitive().max_eigenval()).fold(0.0, |a, b| a.max(b));

    if dx / dt < max_lambda {
        panic!("dx/dt too small!")
    }

    let mut u1 = vec![Conserved(0.0, 0.0, 0.0); n];

    // Update interior cells.
    for i in 1..n-1 {
        let prim_right_iph = u[i+1].to_primitive();
        let prim_left_iph = u[i].to_primitive();
        let prim_right_imh = u[i].to_primitive();
        let prim_left_imh = u[i-1].to_primitive();

        let fr_iph = prim_right_iph.get_fluxes();
        let fl_iph = prim_left_iph.get_fluxes();
        let fr_imh = prim_right_imh.get_fluxes();
        let fl_imh = prim_left_imh.get_fluxes();

        let ur_iph = u[i+1];
        let ul_iph = u[i];
        let ur_imh = u[i];
        let ul_imh = u[i-i];

        // Characteristic speeds for iph.
        let sr_p_iph = prim_right_iph.max_eigenval_signed();
        let sl_p_iph = prim_left_iph.max_eigenval_signed();
        let sr_m_iph = prim_right_iph.min_eigenval_signed();
        let sl_m_iph = prim_left_iph.min_eigenval_signed();

        let alpha_p_iph = sr_p_iph.max(sl_p_iph.max(0.0));
        let alpha_m_iph = sr_m_iph.min(sl_m_iph.min(0.0));

        // Characteristic speeds for imh.
        let sr_p_imh = prim_right_imh.max_eigenval_signed();
        let sl_p_imh = prim_left_imh.max_eigenval_signed();
        let sr_m_imh = prim_right_imh.min_eigenval_signed();
        let sl_m_imh = prim_left_imh.min_eigenval_signed();

        let alpha_p_imh = sr_p_imh.max(sl_p_imh.max(0.0));
        let alpha_m_imh = sr_m_imh.min(sl_m_imh.min(0.0));

        // FHLL.
        let fhll_iph = ( alpha_p_iph * fl_iph - alpha_m_iph * fr_iph + alpha_p_iph * alpha_m_iph * ( ur_iph - ul_iph )) / ( alpha_p_iph - alpha_m_iph );
        let fhll_imh = ( alpha_p_imh * fl_imh - alpha_m_imh * fr_imh + alpha_p_imh * alpha_m_imh * ( ur_imh - ul_imh )) / ( alpha_p_imh - alpha_m_imh );

        u1[i] = u[i] - (dt / dx) * (fhll_iph - fhll_imh);
    }

    u1[0]   = u[0];
    u1[n-1] = u[n-1];
    u1
}

