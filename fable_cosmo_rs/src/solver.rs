#![forbid(unsafe_code)]
#![deny(warnings)]
#![allow(non_snake_case)]

//! The CVODE driver: one BDF + Newton + dense integration in e-folds.
//!
//! The vendored pure-Rust SUNDIALS 7.8.0 is the only integrator; nothing in
//! this crate ever steps the system by hand.

use std::any::Any;

use cvode_rs::prelude::*;

use crate::model::{
    bridge_closed_form, log_bilinear_closed_form, sample_from_state, Params, Sample, ABS_TOL,
    IDX_BRIDGE, IDX_LOG_S, MAX_STEPS_PER_CALL, N_EQ, REL_TOL,
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

fn with_data_mut<R>(v: &N_Vector, f: impl FnOnce(&mut [f64]) -> R) -> Option<R> {
    let mut d = N_VGetArrayPointer(v)?;
    Some(f(&mut d))
}

fn read_state(v: &N_Vector) -> Result<[f64; 2], String> {
    let d = N_VGetArrayPointer(v)
        .ok_or_else(|| "N_VGetArrayPointer returned None for y".to_string())?;
    Ok([d[IDX_LOG_S], d[IDX_BRIDGE]])
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

    let y = N_VNew_Serial(N_EQ, &ctx)
        .ok_or_else(|| "N_VNew_Serial(y) returned None".to_string())?;
    let initial = [
        log_bilinear_closed_form(n_start, p),
        bridge_closed_form(n_start, p),
    ];
    with_data_mut(&y, |d| d.copy_from_slice(&initial))
        .ok_or_else(|| "N_VGetArrayPointer returned None for y".to_string())?;

    let abstol = N_VNew_Serial(N_EQ, &ctx)
        .ok_or_else(|| "N_VNew_Serial(abstol) returned None".to_string())?;
    with_data_mut(&abstol, |d| d.copy_from_slice(&ABS_TOL))
        .ok_or_else(|| "N_VGetArrayPointer returned None for abstol".to_string())?;

    let cv = CVodeCreate(CV_BDF, &ctx)
        .ok_or_else(|| "CVodeCreate(CV_BDF) returned None".to_string())?;

    let mut f = CVodeInit(&cv, rhs::rhs, n_start, &y);
    if f != CV_SUCCESS {
        return Err(format!("CVodeInit failed: {f}"));
    }
    f = CVodeSVtolerances(&cv, REL_TOL, &abstol);
    if f != CV_SUCCESS {
        return Err(format!("CVodeSVtolerances failed: {f}"));
    }

    let a_mat = SUNDenseMatrix(N_EQ, N_EQ, &ctx)
        .ok_or_else(|| "SUNDenseMatrix returned None".to_string())?;
    let ls = SUNLinSol_Dense(&y, &a_mat, &ctx)
        .ok_or_else(|| "SUNLinSol_Dense returned None".to_string())?;
    f = CVodeSetLinearSolver(&cv, &ls, Some(&a_mat));
    if f != CV_SUCCESS {
        return Err(format!("CVodeSetLinearSolver failed: {f}"));
    }

    let user: Box<dyn Any> = Box::new(*p);
    f = CVodeSetUserData(&cv, Some(user));
    if f != CV_SUCCESS {
        return Err(format!("CVodeSetUserData failed: {f}"));
    }
    f = CVodeSetMaxNumSteps(&cv, MAX_STEPS_PER_CALL);
    if f != CV_SUCCESS {
        return Err(format!("CVodeSetMaxNumSteps failed: {f}"));
    }
    f = CVodeSetStopTime(&cv, n_end);
    if f != CV_SUCCESS {
        return Err(format!("CVodeSetStopTime failed: {f}"));
    }

    let mut out = Vec::with_capacity(samples);
    out.push(sample_from_state(n_start, &initial, p));

    let step = (n_end - n_start) / (samples as f64 - 1.0);
    let mut n_now = n_start;
    for k in 1..samples {
        let mut tout = n_start + step * k as f64;
        if tout > n_end {
            tout = n_end;
        }
        let flag = CVode(&cv, tout, &y, &mut n_now, CV_NORMAL);
        if flag < 0 {
            return Err(format!("CVode failed with flag {flag} at N = {n_now}"));
        }
        let state = read_state(&y)?;
        out.push(sample_from_state(n_now, &state, p));
    }

    let mut stats = Stats::default();
    let flag = CVodeGetNumSteps(&cv, &mut stats.n_steps);
    if flag != CV_SUCCESS {
        return Err(format!("CVodeGetNumSteps failed: {flag}"));
    }
    let flag = CVodeGetNumRhsEvals(&cv, &mut stats.n_rhs);
    if flag != CV_SUCCESS {
        return Err(format!("CVodeGetNumRhsEvals failed: {flag}"));
    }

    Ok(Run {
        samples: out,
        stats,
    })
}
