#[derive(Copy, Clone)]
pub struct Conserved (pub f64, pub f64);

#[derive(Copy, Clone)]
pub struct Primitive (pub f64, pub f64);

pub static KAPPA: f64 = 1.0;
pub static GAMMA_LAW_INDEX: f64 = 5.0 / 3.0;


impl Conserved {
    pub fn get_density (self) -> f64 { self.0 }
    pub fn get_momentum (self) -> f64 { self.1 }
    pub fn get_velocity (self) -> f64 { self.get_momentum() / self.get_density() }
    pub fn to_primitive (self) -> Primitive {
        Primitive(self.0, self.get_velocity())
    }
}


impl Primitive {
    pub fn get_density (self) -> f64 { self.0 }
    pub fn get_velocity (self) -> f64 { self.1 }
    pub fn to_conserved (self) -> Conserved {
        Conserved(self.0, self.0 * self.1)
    }

    pub fn get_fluxes (self) -> Conserved {
        Conserved(self.0 * self.1, self.0 * self.1.powi(2) + KAPPA * self.0.powf(GAMMA_LAW_INDEX))
    }

    pub fn sound_speed (self) -> f64 {
        (KAPPA * GAMMA_LAW_INDEX * self.get_density().powf(GAMMA_LAW_INDEX - 1.0)).sqrt()
    }

    pub fn eigen_val_plus (self) -> f64 {
        self.get_velocity() + self.sound_speed()
    }
    pub fn eigen_val_minus (self) -> f64 {
        self.get_velocity() - self.sound_speed()
    }
}


// ============================================================================
impl std::ops::Add<Conserved> for Conserved { type Output = Self; fn add(self, u: Conserved) -> Conserved { Conserved(self.0 + u.0, self.1 + u.1) } }
impl std::ops::Sub<Conserved> for Conserved { type Output = Self; fn sub(self, u: Conserved) -> Conserved { Conserved(self.0 - u.0, self.1 - u.1) } }
impl std::ops::Mul<Conserved> for f64 { type Output = Conserved; fn mul(self, u: Conserved) -> Conserved { Conserved(self * u.0, self * u.1) } }
impl std::ops::Mul<f64> for Conserved { type Output = Conserved; fn mul(self, a: f64) -> Conserved { Conserved(self.0 * a, self.1 * a) } }
impl std::ops::Div<f64> for Conserved { type Output = Conserved; fn div(self, a: f64) -> Conserved { Conserved(self.0 / a, self.1 / a) } }



pub fn flux_split(u: Vec<Conserved>, dx: f64, dt: f64) -> Vec<Conserved> {
    let n = u.len();
    let mut u1 = vec![Conserved(0.0, 0.0); n];

    // Finding max lambda.
    let mut vec_lambdas = vec![];
    for j in 0..n {
        vec_lambdas.push(u[j].to_primitive().eigen_val_minus().abs());
        vec_lambdas.push(u[j].to_primitive().eigen_val_plus().abs());
    }
    let mut max_lambda: f64 = 0.0;
    for k in vec_lambdas {
        if k > max_lambda {
            max_lambda = k;
        } else {}
    }

    if ( dx / dt ) < max_lambda {
        panic!("dx/dt too small!")
    }

    // Update interior.
    for i in 1..n-1 {
        let prim_right = u[i+1].to_primitive();
        let prim_left = u[i-1].to_primitive();
        let fr = prim_right.get_fluxes();
        let fl = prim_left.get_fluxes();

        let flux_cor_term = max_lambda * ( u[i-1] - 2.0 * u[i] + u[i+1] );

        // u1[i] = u[i] - ( dt / dx ) * ( fr - fl );
        u1[i] = u[i] - 0.5 * ( dt / dx ) * ( fr - fl - flux_cor_term );
    }

    // Update boundaries.
    u1[0] = u[0];
    u1[n-1] = u[n-1];

    u1
}

