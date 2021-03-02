



/**
 * Conserved hydrodynamic quantities: density, momentum, entropy density
 */
#[derive(Copy, Clone)]
pub struct Conserved(pub f64, pub f64, pub f64, pub f64);




/**
 * Primitive hydrodynamic quantities: density, velocity, specific internal energy
 */
#[derive(Copy, Clone)]
pub struct Primitive(pub f64, pub f64, pub f64, pub f64);




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

    pub fn passive_scalar_density(self) -> f64 {
        self.3
    }

    pub fn velocity(self) -> f64 {
        self.momentum() / self.density()
    }

    pub fn specific_internal_energy(self) -> f64 {
        self.2 / self.0.powf(2.0 - GAMMA_LAW_INDEX) 
    }

    pub fn passive_scalar(self) -> f64 {
        self.3 / self.0
    }

    pub fn to_primitive(self) -> Primitive {
        Primitive(self.density(), self.velocity(), self.specific_internal_energy(), self.passive_scalar())
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

    pub fn pressure(self) -> f64 {
        self.0 * self.2 * (GAMMA_LAW_INDEX - 1.0)
    }

    pub fn passive_scalar(self) -> f64 {
        self.3
    }

    pub fn to_conserved(self) -> Conserved {
        Conserved(self.density(), self.density() * self.velocity(), self.density() * self.specific_internal_energy() * self.density().powf(1.0 - GAMMA_LAW_INDEX), self.density() * self.passive_scalar())
    }

    pub fn get_fluxes(self) -> Conserved {
        let mass_flux     = self.density() * self.velocity();
        let momentum_flux = self.density() * self.velocity().powi(2) + self.density() * self.specific_internal_energy() * (GAMMA_LAW_INDEX - 1.0);
        let entropy_flux  = self.to_conserved().entropy_density() * self.velocity();
        let scalar_flux   = self.density() * self.passive_scalar() * self.velocity();
        Conserved(mass_flux, momentum_flux, entropy_flux, scalar_flux)
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
impl std::ops::Add<Conserved> for Conserved { type Output = Self; fn add(self, u: Conserved) -> Conserved { Conserved(self.0 + u.0, self.1 + u.1, self.2 + u.2, self.3 + u.3) } }
impl std::ops::Sub<Conserved> for Conserved { type Output = Self; fn sub(self, u: Conserved) -> Conserved { Conserved(self.0 - u.0, self.1 - u.1, self.2 - u.2, self.3 - u.3) } }
impl std::ops::Mul<Conserved> for f64 { type Output = Conserved; fn mul(self, u: Conserved) -> Conserved { Conserved(self * u.0, self * u.1, self * u.2, self * u.3) } }
impl std::ops::Mul<f64> for Conserved { type Output = Conserved; fn mul(self, a: f64) -> Conserved { Conserved(self.0 * a, self.1 * a, self.2 * a, self.3 * a) } }
impl std::ops::Div<f64> for Conserved { type Output = Conserved; fn div(self, a: f64) -> Conserved { Conserved(self.0 / a, self.1 / a, self.2 / a, self.2 / a) } }




// ============================================================================
pub fn flux_split(u: Vec<Conserved>, dx: f64, dt: f64) -> Vec<Conserved> {
    let n = u.len();
    let max_lambda: f64 = u.iter().map(|u| u.to_primitive().max_eigenval()).fold(0.0, |a, b| a.max(b));

    if dx / dt < max_lambda {
        panic!("dx/dt too small!")
    }

    let mut u1 = vec![Conserved(0.0, 0.0, 0.0, 0.0); n];

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

    let mut u1 = vec![Conserved(0.0, 0.0, 0.0, 0.0); n];

    // Update interior cells.
    for i in 1..n-1 {
        let pr_iph = u[i+1].to_primitive();
        let pl_iph = u[i  ].to_primitive();
        let pr_imh = u[i  ].to_primitive();
        let pl_imh = u[i-1].to_primitive();

        let fr_iph = pr_iph.get_fluxes();
        let fl_iph = pl_iph.get_fluxes();
        let fr_imh = pr_imh.get_fluxes();
        let fl_imh = pl_imh.get_fluxes();

        let ur_iph = u[i+1];
        let ul_iph = u[i  ];
        let ur_imh = u[i  ];
        let ul_imh = u[i-1]; //had [i-i] here

        // Characteristic speeds for iph.
        let sr_premax_iph = pr_iph.max_eigenval_signed();
        let sl_premax_iph = pl_iph.max_eigenval_signed();
        let sr_premin_iph = pr_iph.min_eigenval_signed();
        let sl_premin_iph = pl_iph.min_eigenval_signed();

        let sr_iph = sr_premax_iph.max(sl_premax_iph.max(0.0));
        let sl_iph = sr_premin_iph.min(sl_premin_iph.min(0.0));

        // Characteristic speeds for imh.
        let sr_premax_imh = pr_imh.max_eigenval_signed();
        let sl_premax_imh = pl_imh.max_eigenval_signed();
        let sr_premin_imh = pr_imh.min_eigenval_signed();
        let sl_premin_imh = pl_imh.min_eigenval_signed();

        let sr_imh = sr_premax_imh.max(sl_premax_imh.max(0.0));
        let sl_imh = sr_premin_imh.min(sl_premin_imh.min(0.0));

        // FHLL.
        let fhll_iph = ( sr_iph * fl_iph - sl_iph * fr_iph + sr_iph * sl_iph * ( ur_iph - ul_iph )) / ( sr_iph - sl_iph );
        let fhll_imh = ( sr_imh * fl_imh - sl_imh * fr_imh + sr_imh * sl_imh * ( ur_imh - ul_imh )) / ( sr_imh - sl_imh );

        u1[i] = u[i] - (dt / dx) * (fhll_iph - fhll_imh);
    }

    u1[0]   = u[0];
    u1[n-1] = u[n-1];
    u1
}


// ============================================================================
pub fn hllc_solver(u: Vec<Conserved>, dx: f64, dt: f64) -> Vec<Conserved> {
    let n = u.len();
    let max_lambda: f64 = u.iter().map(|u| u.to_primitive().max_eigenval()).fold(0.0, |a, b| a.max(b));
    if dx / dt < max_lambda {
        panic!("dx/dt too small!")
    }

    let mut u1 = vec![Conserved(0.0, 0.0, 0.0, 0.0); n];

    // Update interior cells.
    for i in 1..n-1 {
        let pr_iph = u[i+1].to_primitive();
        let pl_iph = u[i  ].to_primitive();
        let pr_imh = u[i  ].to_primitive();
        let pl_imh = u[i-1].to_primitive();

        let fr_iph = pr_iph.get_fluxes();
        let fl_iph = pl_iph.get_fluxes();
        let fr_imh = pr_imh.get_fluxes();
        let fl_imh = pl_imh.get_fluxes();

        let ur_iph = u[i+1];
        let ul_iph = u[i  ];
        let ur_imh = u[i  ];
        let ul_imh = u[i-1];

        // Characteristic speeds for iph.
        let sr_premax_iph = pr_iph.max_eigenval_signed();
        let sl_premax_iph = pl_iph.max_eigenval_signed();
        let sr_premin_iph = pr_iph.min_eigenval_signed();
        let sl_premin_iph = pl_iph.min_eigenval_signed();

        let sr_iph = sr_premax_iph.max(sl_premax_iph.max(0.0));
        let sl_iph = sr_premin_iph.min(sl_premin_iph.min(0.0));

        // Characteristic speeds for imh.
        let sr_premax_imh = pr_imh.max_eigenval_signed();
        let sl_premax_imh = pl_imh.max_eigenval_signed();
        let sr_premin_imh = pr_imh.min_eigenval_signed();
        let sl_premin_imh = pl_imh.min_eigenval_signed();

        let sr_imh = sr_premax_imh.max(sl_premax_imh.max(0.0));
        let sl_imh = sr_premin_imh.min(sl_premin_imh.min(0.0));

        // HLLC stuff
        let (rhor, vr, presr) = (pr_iph.density(), pr_iph.velocity(), pr_iph.pressure());
        let (rhol, vl, presl) = (pl_iph.density(), pl_iph.velocity(), pl_iph.pressure());
        let smid_iph  = (presr - presl + rhol * vl * (sl_iph - vl) - rhor * vr *(sr_iph - vr)) / (rhol * (sl_iph - vl) - rhor * (sr_iph - vr));
        let rhomidr   = rhor * (sr_iph - vr) / (sr_iph - smid_iph);
        let rhomidl   = rhol * (sl_iph - vl) / (sl_iph - smid_iph);
        let presmidr  = presr + rhor * (sr_iph - vr) *(smid_iph - vr);
        let presmidl  = presl + rhol * (sl_iph - vl) *(smid_iph - vl);
        let epsmidr   = presmidr / rhomidr / (GAMMA_LAW_INDEX - 1.0);
        let epsmidl   = presmidl / rhomidl / (GAMMA_LAW_INDEX - 1.0);
        let umidr_iph = Primitive(rhomidr, smid_iph, epsmidr, pr_iph.passive_scalar()).to_conserved();
        let umidl_iph = Primitive(rhomidl, smid_iph, epsmidl, pl_iph.passive_scalar()).to_conserved();
        let fmidr_iph = fr_iph + (umidr_iph - ur_iph) * sr_iph;
        let fmidl_iph = fl_iph + (umidl_iph - ul_iph) * sl_iph;

        let (rhor, vr, presr) = (pr_imh.density(), pr_imh.velocity(), pr_imh.pressure());
        let (rhol, vl, presl) = (pl_imh.density(), pl_imh.velocity(), pl_imh.pressure());
        let smid_imh  = (presr - presl + rhol * vl * (sl_imh - vl) - rhor * vr *(sr_imh - vr)) / (rhol * (sl_imh - vl) - rhor * (sr_imh - vr));
        let rhomidr   = rhor * (sr_imh - vr) / (sr_imh - smid_imh);
        let rhomidl   = rhol * (sl_imh - vl) / (sl_imh - smid_imh);
        let presmidr  = presr + rhor * (sr_imh - vr) *(smid_imh - vr);
        let presmidl  = presl + rhol * (sl_imh - vl) *(smid_imh - vl);
        let epsmidr   = presmidr / rhomidr / (GAMMA_LAW_INDEX - 1.0);
        let epsmidl   = presmidl / rhomidl / (GAMMA_LAW_INDEX - 1.0);
        let umidr_imh = Primitive(rhomidr, smid_imh, epsmidr, pr_imh.passive_scalar()).to_conserved();
        let umidl_imh = Primitive(rhomidl, smid_imh, epsmidl, pl_imh.passive_scalar()).to_conserved();
        let fmidr_imh = fr_imh + (umidr_imh - ur_imh) * sr_imh;
        let fmidl_imh = fl_imh + (umidl_imh - ul_imh) * sl_imh;

        let fhllc_iph = if 0.0 <= sl_iph {
            fl_iph
        } else if sl_iph < 0.0 && 0.0 <= smid_iph {
            fmidl_iph
        } else if smid_iph < 0.0 && 0.0 <= sr_iph {
            fmidr_iph
        } else {
            fr_iph
        };

        let fhllc_imh = if 0.0 <= sl_imh {
            fl_imh
        } else if sl_imh < 0.0 && 0.0 <= smid_imh {
            fmidl_imh
        } else if smid_imh < 0.0 && 0.0 <= sr_imh {
            fmidr_imh
        } else {
            fr_imh
        };

        u1[i] = u[i] - (dt / dx) * (fhllc_iph - fhllc_imh);
    }

    u1[0]   = u[0];
    u1[n-1] = u[n-1];
    u1
}

