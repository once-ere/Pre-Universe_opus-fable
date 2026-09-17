#![forbid(unsafe_code)]
#![deny(warnings)]
#![allow(non_snake_case)]

//! CSV and summary writers.
//!
//! Every floating-point field goes through the engine's own `fmt_e`/`fmt_g`,
//! which reproduce C's `%.*e` and `%.*g` exactly, so the output is stable
//! across hosts and byte-comparable between runs.

use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use sundials_core::sundials_utils::{fmt_e, fmt_g};

use crate::model::{Params, Sample};
use crate::solver::Stats;

pub const CSV_HEADER: &str = "e_folds,scale_factor,redshift,bridge,bridge_h,dilution,\
bilinear,density_dust,density_potential,density_fable,pressure_dust,pressure_potential,\
pressure_fable,w_fable,w_dust,w_potential,hubble_over_h0,omega_fable,deceleration,\
continuity_residual,closed_form_residual,bridge_residual";

fn row(s: &Sample) -> String {
    let f = |x: f64| fmt_e(x, 16);
    [
        f(s.e_folds),
        f(s.scale_factor),
        f(s.redshift),
        f(s.bridge),
        f(s.bridge_h),
        f(s.dilution),
        f(s.bilinear),
        f(s.density_dust),
        f(s.density_potential),
        f(s.density_fable),
        f(s.pressure_dust),
        f(s.pressure_potential),
        f(s.pressure_fable),
        f(s.w_fable),
        f(s.w_dust),
        f(s.w_potential),
        f(s.hubble_over_h0),
        f(s.omega_fable),
        f(s.deceleration),
        f(s.continuity_residual),
        f(s.closed_form_residual),
        f(s.bridge_residual),
    ]
    .join(",")
}

/// Write the full background table.
pub fn write_csv(path: &Path, samples: &[Sample]) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        create_dir_all(parent).map_err(|e| format!("create_dir_all({parent:?}) failed: {e}"))?;
    }
    let file = File::create(path).map_err(|e| format!("File::create({path:?}) failed: {e}"))?;
    let mut out = BufWriter::new(file);
    writeln!(out, "{CSV_HEADER}").map_err(|e| format!("write failed: {e}"))?;
    for s in samples {
        writeln!(out, "{}", row(s)).map_err(|e| format!("write failed: {e}"))?;
    }
    out.flush().map_err(|e| format!("flush failed: {e}"))?;
    Ok(())
}

/// Write the machine-readable run summary consumed by the notebook and tests.
#[allow(clippy::too_many_arguments)]
pub fn write_summary(
    path: &Path,
    p: &Params,
    stats: &Stats,
    today: &Sample,
    max_continuity: f64,
    max_closed_form: f64,
    max_bridge: f64,
) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        create_dir_all(parent).map_err(|e| format!("create_dir_all({parent:?}) failed: {e}"))?;
    }
    let file = File::create(path).map_err(|e| format!("File::create({path:?}) failed: {e}"))?;
    let mut out = BufWriter::new(file);
    let g = |x: f64| fmt_g(x, 17);
    let e = |x: f64| fmt_e(x, 16);
    writeln!(out, "{{").map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "  \"engine\": \"sundials_rs 7.8.0 (pure Rust, Linux x86-64)\",")
        .map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "  \"method\": \"CVODE BDF + Newton + dense\",")
        .map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "  \"parameters\": {{").map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"omega_r0\": {},", g(p.omega_r0)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"omega_b0\": {},", g(p.omega_b0)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"mass_m\": {},", g(p.mass_m)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"lambda\": {},", g(p.lambda)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"index_n\": {},", g(p.index_n)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"xi\": {},", g(p.xi)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"gamma_flow\": {},", g(p.gamma_flow)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"bridge_s0\": {},", g(p.bridge_s0)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"bilinear_s0\": {}", g(p.bilinear_s0)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "  }},").map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "  \"today\": {{").map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"w_fable\": {},", g(today.w_fable)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"w_dust\": {},", g(today.w_dust)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"w_potential\": {},", g(today.w_potential)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"omega_fable\": {},", g(today.omega_fable)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"deceleration\": {},", g(today.deceleration)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"dilution\": {}", g(today.dilution)).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "  }},").map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "  \"invariants\": {{").map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"max_continuity_residual\": {},", e(max_continuity))
        .map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"max_closed_form_residual\": {},", e(max_closed_form))
        .map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"max_bridge_residual\": {}", e(max_bridge))
        .map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "  }},").map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "  \"cvode\": {{").map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"n_steps\": {},", stats.n_steps).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "    \"n_rhs_evals\": {}", stats.n_rhs).map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "  }}").map_err(|x| format!("write failed: {x}"))?;
    writeln!(out, "}}").map_err(|x| format!("write failed: {x}"))?;
    out.flush().map_err(|x| format!("flush failed: {x}"))?;
    Ok(())
}
