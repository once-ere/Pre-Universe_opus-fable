//! The CVODE driver: one BDF + Newton + dense integration in e-folds.
//!
//! The vendored pure-Rust SUNDIALS 7.8.0 is the only integrator; nothing in
//! this crate ever steps the system by hand. Every flag the engine returns is
//! checked and turned into a named error; none is discarded.

use std::any::Any;

use cvode_rs::prelude::*;

use crate::model::{
    closed_form_state, sample_from_state, Params, Sample, State, ABS_TOL, IDX_BRIDGE,
    IDX_LOG_S, MAX_STEPS_PER_CALL, N_EQ, REL_TOL,
};
use crate::rhs;

/// Accumulated CVODE statistics.
#[derive(Clone, Copy, Debug, Default)]
pub struct Stats {
    pub n_steps: i64,
    pub n_rhs: i64,
}

/// The output grid of one run.
pub struct Run {
    pub samples: Vec<Sample>,
    pub stats: Stats,
}

/// Turn a SUNDIALS return flag into a named error.
fn check(flag: i32, what: &str) -> Result<(), String> {
    if flag == CV_SUCCESS {
        Ok(())
    } else {
        Err(format!("{what} failed: {flag}"))
    }
}

/// Copy a slice into a vector under a scoped borrow, so the `RefMut` from
/// `N_VGetArrayPointer` is never held across an engine call.
fn fill(v: &N_Vector, values: &[f64], what: &str) -> Result<(), String> {
    let mut d = N_VGetArrayPointer(v)
        .ok_or_else(|| format!("N_VGetArrayPointer returned None for {what}"))?;
    d.copy_from_slice(values);
    Ok(())
}

fn read_state(v: &N_Vector) -> Result<State, String> {
    let d = N_VGetArrayPointer(v)
        .ok_or_else(|| "N_VGetArrayPointer returned None for y".to_string())?;
    let mut state = [0.0; crate::model::N_STATE];
    state[IDX_LOG_S] = d[IDX_LOG_S];
    state[IDX_BRIDGE] = d[IDX_BRIDGE];
    Ok(state)
}

/// Integrate the background from `n_start` to `n_end` on a uniform grid of
/// `samples` points, starting from the exact closed-form state at `n_start`.
pub fn integrate(
    p: &Params,
    n_start: f64,
    n_end: f64,
    samples: usize,
) -> Result<Run, String> {
    p.validate()?;
    if !(n_end > n_start) {
        return Err(format!(
            "integrate: n_end ({n_end}) must exceed n_start ({n_start})"
        ));
    }
    if samples < 2 {
        return Err("integrate: at least two output samples are required".to_string());
    }

    let mut ctx_out: Option<SUNContext> = None;
    let rc = SUNContext_Create(SUN_COMM_NULL, &mut ctx_out);
    if rc != 0 {
        return Err(format!("SUNContext_Create failed: {rc}"));
    }
    let ctx = ctx_out.ok_or_else(|| "SUNContext_Create returned no context".to_string())?;

    let initial = closed_form_state(n_start, p);
    let y = N_VNew_Serial(N_EQ, &ctx)
        .ok_or_else(|| "N_VNew_Serial(y) returned None".to_string())?;
    fill(&y, &initial, "y")?;

    let abstol = N_VNew_Serial(N_EQ, &ctx)
        .ok_or_else(|| "N_VNew_Serial(abstol) returned None".to_string())?;
    fill(&abstol, &ABS_TOL, "abstol")?;

    let cv = CVodeCreate(CV_BDF, &ctx)
        .ok_or_else(|| "CVodeCreate(CV_BDF) returned None".to_string())?;
    check(CVodeInit(&cv, rhs::rhs, n_start, &y), "CVodeInit")?;
    check(CVodeSVtolerances(&cv, REL_TOL, &abstol), "CVodeSVtolerances")?;

    let a_mat = SUNDenseMatrix(N_EQ, N_EQ, &ctx)
        .ok_or_else(|| "SUNDenseMatrix returned None".to_string())?;
    let ls = SUNLinSol_Dense(&y, &a_mat, &ctx)
        .ok_or_else(|| "SUNLinSol_Dense returned None".to_string())?;
    check(CVodeSetLinearSolver(&cv, &ls, Some(&a_mat)), "CVodeSetLinearSolver")?;

    let user: Box<dyn Any> = Box::new(*p);
    check(CVodeSetUserData(&cv, Some(user)), "CVodeSetUserData")?;
    check(CVodeSetMaxNumSteps(&cv, MAX_STEPS_PER_CALL), "CVodeSetMaxNumSteps")?;
    check(CVodeSetStopTime(&cv, n_end), "CVodeSetStopTime")?;

    let mut out = Vec::with_capacity(samples);
    out.push(sample_from_state(n_start, &initial, p));

    let step = (n_end - n_start) / (samples as f64 - 1.0);
    let mut n_now = n_start;
    for k in 1..samples {
        // The last grid point is n_end up to rounding; the clamp guarantees it.
        let tout = (n_start + step * k as f64).min(n_end);
        let flag = CVode(&cv, tout, &y, &mut n_now, CV_NORMAL);
        if flag < 0 {
            return Err(format!("CVode failed with flag {flag} at N = {n_now}"));
        }
        let state = read_state(&y)?;
        out.push(sample_from_state(n_now, &state, p));
    }

    let mut stats = Stats::default();
    check(CVodeGetNumSteps(&cv, &mut stats.n_steps), "CVodeGetNumSteps")?;
    check(CVodeGetNumRhsEvals(&cv, &mut stats.n_rhs), "CVodeGetNumRhsEvals")?;

    Ok(Run {
        samples: out,
        stats,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_a_degenerate_range_and_too_few_samples() {
        let p = Params::default();
        assert!(integrate(&p, 0.0, 0.0, 10).is_err());
        assert!(integrate(&p, 1.0, 0.0, 10).is_err());
        assert!(integrate(&p, -1.0, 1.0, 1).is_err());
    }

    #[test]
    fn tracks_the_closed_form_and_lands_exactly_on_the_end_point() {
        let p = Params::default();
        let run = integrate(&p, -7.003, 1.0986122886681098, 201).unwrap();
        assert_eq!(run.samples.len(), 201);
        assert_eq!(run.samples.last().unwrap().e_folds, 1.0986122886681098);
        let worst = run
            .samples
            .iter()
            .map(|s| s.closed_form_residual.abs().max(s.bridge_residual.abs()))
            .fold(0.0_f64, f64::max);
        assert!(worst < 1.0e-9, "worst residual {worst}");
        assert!(run.stats.n_steps > 0 && run.stats.n_rhs >= run.stats.n_steps);
    }
}
