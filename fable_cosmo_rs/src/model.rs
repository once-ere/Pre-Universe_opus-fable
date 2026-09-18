//! The fableSpinor background model.
//!
//! `fableSpinor` is a real 16-component field carrying an irreducible
//! representation of the universal cover of the real O(4,4). Its only
//! nonvanishing background bilinear is the scalar `S = Psi^T C Psi`, because
//! `C` is symmetric while every `C gamma^a` is antisymmetric, so the vector
//! current of a commuting real spinor vanishes identically.
//!
//! The homogeneous Dirac equation, including the gpt-5.6_bridge scalar
//! connection term, gives the exact dilution law
//!
//! ```text
//!   dS/dN = -nu(s) S,      nu(s) = 3 - xi (1 - h(s)),   h(s) = 3 s^2 - 2 s^3
//! ```
//!
//! where `N = ln a` and `s` is the bridge field interpolating the connection
//! from Levi-Civita (`s = 0`) to Weitzenboeck (`s = 1`). Every formula in this
//! file is proved symbolically in `wolfram/fable_spinor.wls`.

/// Number of integrated state variables: `[ln S, s]`.
pub const N_STATE: usize = 2;
/// The same count as the `sunindextype` the engine's constructors take.
pub const N_EQ: i64 = N_STATE as i64;

/// Index of `ln S` in the state vector.
pub const IDX_LOG_S: usize = 0;
/// Index of the bridge field `s` in the state vector.
pub const IDX_BRIDGE: usize = 1;

/// The integrated state `[ln S, s]`.
pub type State = [f64; N_STATE];

/// Relative tolerance handed to CVODE.
pub const REL_TOL: f64 = 1.0e-12;
/// Absolute tolerances handed to CVODE, one per state variable.
pub const ABS_TOL: State = [1.0e-14, 1.0e-14];
/// Step ceiling for one `CVode` call.
pub const MAX_STEPS_PER_CALL: i64 = 500_000;

/// The Unite-only constant-`w` benchmark this model must reproduce exactly in
/// its massless, teleparallel limit (Camilleri et al., arXiv:2609.05053v2,
/// Section 7.0.3).
pub const BENCHMARK_W: f64 = -0.764;

/// Model parameters. All densities are in units of today's critical density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Params {
    /// Radiation density today.
    pub omega_r0: f64,
    /// Baryon density today.
    pub omega_b0: f64,
    /// Mass coefficient `M`: the dust-like part of the fable sector today.
    pub mass_m: f64,
    /// Potential coefficient `lambda`: the dark-energy-like part today.
    pub lambda: f64,
    /// Potential index `n` in `V(S) = lambda S^n`.
    pub index_n: f64,
    /// Torsion coupling `xi`: the strength of the second mechanism.
    pub xi: f64,
    /// Logistic rate of the bridge field.
    pub gamma_flow: f64,
    /// Bridge field today (`N = 0`).
    pub bridge_s0: f64,
    /// Bilinear today (`N = 0`), the normalisation of `S`.
    pub bilinear_s0: f64,
}

impl Default for Params {
    /// The published scenario: a unified dark sector whose potential index is
    /// the benchmark `n = 1 + BENCHMARK_W`, carrying a dust-like component of
    /// 0.256 and a dark-energy-like component of 0.694910 today, alongside
    /// 0.049 of baryons and 9e-5 of radiation in a spatially flat universe.
    fn default() -> Self {
        let omega_r0 = 9.0e-5;
        let omega_b0 = 0.049;
        let mass_m = 0.256;
        Params {
            omega_r0,
            omega_b0,
            mass_m,
            lambda: 1.0 - omega_r0 - omega_b0 - mass_m,
            index_n: 1.0 + BENCHMARK_W,
            xi: 0.15,
            gamma_flow: 0.35,
            bridge_s0: 0.5,
            bilinear_s0: 1.0,
        }
    }
}

impl Params {
    /// Fail loudly rather than integrate a meaningless model.
    pub fn validate(&self) -> Result<(), String> {
        let all = [
            self.omega_r0,
            self.omega_b0,
            self.mass_m,
            self.lambda,
            self.index_n,
            self.xi,
            self.gamma_flow,
            self.bridge_s0,
            self.bilinear_s0,
        ];
        if all.iter().any(|v| !v.is_finite()) {
            return Err("every parameter must be finite".to_string());
        }
        if self.omega_r0 < 0.0 || self.omega_b0 < 0.0 || self.mass_m < 0.0 {
            return Err("densities must be nonnegative".to_string());
        }
        if self.lambda <= 0.0 {
            return Err("the potential coefficient lambda must be positive".to_string());
        }
        if !(0.0 < self.index_n && self.index_n < 1.0) {
            return Err("the potential index n must lie strictly in (0, 1)".to_string());
        }
        if !(0.0 < self.bridge_s0 && self.bridge_s0 < 1.0) {
            return Err("the bridge field today must lie strictly in (0, 1)".to_string());
        }
        if self.gamma_flow <= 0.0 {
            return Err("the bridge flow rate must be positive".to_string());
        }
        if self.bilinear_s0 <= 0.0 {
            return Err("the bilinear today must be positive".to_string());
        }
        let total = self.omega_r0 + self.omega_b0 + self.mass_m + self.lambda;
        if (total - 1.0).abs() > 1.0e-12 {
            return Err(format!(
                "the universe is not flat: densities today sum to {total}"
            ));
        }
        Ok(())
    }

    /// The same model with the torsion coupling switched off: the limit in
    /// which the dark-energy-like component has exactly `w = n - 1`.
    pub fn torsion_free(&self) -> Params {
        Params { xi: 0.0, ..*self }
    }
}

/// The smoothstep `h(s) = 3 s^2 - 2 s^3` that dials the gpt-5.6_bridge from
/// Levi-Civita (`h = 0`) to Weitzenboeck (`h = 1`).
pub fn bridge_h(s: f64) -> f64 {
    s * s * (3.0 - 2.0 * s)
}

/// The dilution exponent `nu(s) = 3 - xi (1 - h(s))`. It equals 3 exactly at
/// the Weitzenboeck endpoint and whenever the torsion coupling is switched off.
pub fn dilution_exponent(s: f64, xi: f64) -> f64 {
    3.0 - xi * (1.0 - bridge_h(s))
}

/// The logistic flow of the bridge field, `ds/dN = gamma s (1 - s)`.
pub fn bridge_flow(s: f64, gamma_flow: f64) -> f64 {
    gamma_flow * s * (1.0 - s)
}

/// Closed-form bridge field, anchored at `s(0) = s0`.
pub fn bridge_closed_form(n: f64, p: &Params) -> f64 {
    let g = (p.gamma_flow * n).exp();
    p.bridge_s0 * g / (1.0 - p.bridge_s0 + p.bridge_s0 * g)
}

/// The quadrature `(1/gamma)(ln s + s - s^2)` whose derivative along the
/// logistic flow is exactly `1 - h(s)`.
fn bilinear_quadrature(s: f64, gamma_flow: f64) -> f64 {
    (s.ln() + s - s * s) / gamma_flow
}

/// Closed-form `ln S(N)`, the exact solution of the dilution law.
pub fn log_bilinear_closed_form(n: f64, p: &Params) -> f64 {
    let s_n = bridge_closed_form(n, p);
    p.bilinear_s0.ln() - 3.0 * n
        + p.xi
            * (bilinear_quadrature(s_n, p.gamma_flow)
                - bilinear_quadrature(p.bridge_s0, p.gamma_flow))
}

/// The exact closed-form state at e-fold `n`. Used to seed the integrator and
/// to report "today" independently of where the output grid happens to land.
pub fn closed_form_state(n: f64, p: &Params) -> State {
    let mut state = [0.0; N_STATE];
    state[IDX_LOG_S] = log_bilinear_closed_form(n, p);
    state[IDX_BRIDGE] = bridge_closed_form(n, p);
    state
}

/// Every derived background quantity at one sample.
#[derive(Clone, Copy, Debug)]
pub struct Sample {
    pub e_folds: f64,
    pub scale_factor: f64,
    pub redshift: f64,
    pub bridge: f64,
    pub bridge_h: f64,
    pub dilution: f64,
    pub bilinear: f64,
    /// Dust-like part of the fable sector, `M S`.
    pub density_dust: f64,
    /// Dark-energy-like part of the fable sector, `lambda S^n`.
    pub density_potential: f64,
    pub density_fable: f64,
    pub pressure_dust: f64,
    pub pressure_potential: f64,
    pub pressure_fable: f64,
    /// Equation of state of the whole fable sector.
    pub w_fable: f64,
    /// Equation of state of the dust-like part, `nu/3 - 1`.
    pub w_dust: f64,
    /// Equation of state of the dark-energy-like part, `n nu/3 - 1`.
    pub w_potential: f64,
    pub hubble_over_h0: f64,
    pub omega_fable: f64,
    pub deceleration: f64,
    /// Residual of `drho/dN + 3 (rho + p)`, divided by `rho`.
    pub continuity_residual: f64,
    /// Difference between the integrated and closed-form `ln S`.
    pub closed_form_residual: f64,
    /// Difference between the integrated and closed-form bridge field.
    pub bridge_residual: f64,
}

/// Build every derived quantity from the integrated state `[ln S, s]`.
pub fn sample_from_state(e_folds: f64, state: &State, p: &Params) -> Sample {
    let log_s = state[IDX_LOG_S];
    let s_bridge = state[IDX_BRIDGE];
    let bilinear = log_s.exp();
    let nu = dilution_exponent(s_bridge, p.xi);

    let density_dust = p.mass_m * bilinear;
    let density_potential = p.lambda * bilinear.powf(p.index_n);
    let density_fable = density_dust + density_potential;

    // p_i = w_i rho_i with w_dust = nu/3 - 1 and w_potential = n nu/3 - 1: the
    // unique pressures for which the sector is covariantly conserved under the
    // dilution law dS/dN = -nu S.
    let w_dust = nu / 3.0 - 1.0;
    let w_potential = p.index_n * nu / 3.0 - 1.0;
    let pressure_dust = w_dust * density_dust;
    let pressure_potential = w_potential * density_potential;
    let pressure_fable = pressure_dust + pressure_potential;

    let scale_factor = e_folds.exp();
    let density_radiation = p.omega_r0 / scale_factor.powi(4);
    let density_baryon = p.omega_b0 / scale_factor.powi(3);
    let density_total = density_radiation + density_baryon + density_fable;
    let pressure_total = density_radiation / 3.0 + pressure_fable;

    // drho_fable/dN = -nu (M S + n lambda S^n), from the dilution law.
    let d_density_d_efolds = -nu * (density_dust + p.index_n * density_potential);
    let continuity_residual = (d_density_d_efolds + 3.0 * (density_fable + pressure_fable))
        / density_fable;

    Sample {
        e_folds,
        scale_factor,
        redshift: 1.0 / scale_factor - 1.0,
        bridge: s_bridge,
        bridge_h: bridge_h(s_bridge),
        dilution: nu,
        bilinear,
        density_dust,
        density_potential,
        density_fable,
        pressure_dust,
        pressure_potential,
        pressure_fable,
        w_fable: pressure_fable / density_fable,
        w_dust,
        w_potential,
        hubble_over_h0: density_total.sqrt(),
        omega_fable: density_fable / density_total,
        deceleration: 0.5 * (1.0 + 3.0 * pressure_total / density_total),
        continuity_residual,
        closed_form_residual: log_s - log_bilinear_closed_form(e_folds, p),
        bridge_residual: s_bridge - bridge_closed_form(e_folds, p),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_parameters_are_valid_and_flat() {
        let p = Params::default();
        p.validate().unwrap();
        let total = p.omega_r0 + p.omega_b0 + p.mass_m + p.lambda;
        assert!((total - 1.0).abs() <= 1.0e-12);
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        let base = Params::default();
        assert!(Params { index_n: 1.5, ..base }.validate().is_err());
        assert!(Params { index_n: 0.0, ..base }.validate().is_err());
        assert!(Params { bridge_s0: 0.0, ..base }.validate().is_err());
        assert!(Params { gamma_flow: -1.0, ..base }.validate().is_err());
        assert!(Params { mass_m: 0.5, ..base }.validate().is_err());
        assert!(Params { xi: f64::NAN, ..base }.validate().is_err());
    }

    #[test]
    fn smoothstep_endpoints() {
        assert_eq!(bridge_h(0.0), 0.0);
        assert_eq!(bridge_h(1.0), 1.0);
        assert_eq!(bridge_h(0.5), 0.5);
    }

    #[test]
    fn dilution_exponent_is_three_at_the_weitzenboeck_endpoint() {
        assert_eq!(dilution_exponent(1.0, 0.37), 3.0);
        assert_eq!(dilution_exponent(0.25, 0.0), 3.0);
    }

    #[test]
    fn closed_form_solves_the_dilution_law() {
        // d(ln S)/dN must equal -nu(s(N)) along the logistic flow.
        let p = Params::default();
        let step = 1.0e-6;
        let mut worst = 0.0_f64;
        for k in 0..=400 {
            let n = -6.0 + 7.0 * k as f64 / 400.0;
            let numerical = (log_bilinear_closed_form(n + step, &p)
                - log_bilinear_closed_form(n - step, &p))
                / (2.0 * step);
            let analytic = -dilution_exponent(bridge_closed_form(n, &p), p.xi);
            worst = worst.max((numerical - analytic).abs());
        }
        assert!(worst < 1.0e-8, "worst residual {worst}");
    }

    #[test]
    fn closed_form_state_is_anchored_today() {
        let p = Params::default();
        let today = closed_form_state(0.0, &p);
        assert_eq!(today[IDX_BRIDGE], p.bridge_s0);
        assert_eq!(today[IDX_LOG_S], p.bilinear_s0.ln());
    }

    #[test]
    fn sector_is_covariantly_conserved() {
        let p = Params::default();
        for k in 0..=50 {
            let n = -7.0 + 8.0 * k as f64 / 50.0;
            let s = sample_from_state(n, &closed_form_state(n, &p), &p);
            assert!(s.continuity_residual.abs() < 1.0e-13, "N = {n}");
        }
    }

    #[test]
    fn torsion_free_limit_reproduces_the_benchmark_exactly() {
        let p = Params::default().torsion_free();
        let s = sample_from_state(0.0, &closed_form_state(0.0, &p), &p);
        assert_eq!(s.w_potential, BENCHMARK_W);
        assert_eq!(s.w_dust, 0.0);
    }

    #[test]
    fn dark_matter_early_dark_energy_late() {
        let p = Params::default().torsion_free();
        let early = sample_from_state(-12.0, &closed_form_state(-12.0, &p), &p);
        let late = sample_from_state(6.0, &closed_form_state(6.0, &p), &p);
        assert!(early.w_fable.abs() < 1.0e-3, "early w = {}", early.w_fable);
        assert!(
            (late.w_fable - BENCHMARK_W).abs() < 1.0e-3,
            "late w = {}",
            late.w_fable
        );
        assert!(early.w_fable > -1.0 && late.w_fable > -1.0);
    }
}
