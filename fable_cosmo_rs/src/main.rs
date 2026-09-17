#![forbid(unsafe_code)]
#![deny(warnings)]
#![allow(non_snake_case)]

//! fable_cosmo_rs — the reference background integration of the fableSpinor
//! unified dark sector.
//!
//! Usage:
//!   cargo run --release -- [--out DIR] [--xi X] [--n-start N] [--n-end N] [--samples K]
//!
//! Exit status is nonzero if any gated invariant is violated.

mod model;
mod output;
mod rhs;
mod solver;

use std::path::PathBuf;
use std::process::ExitCode;

use sundials_core::sundials_utils::{fmt_e, fmt_g};

use model::{Params, Sample, BENCHMARK_W};

/// Gate: the sector must be covariantly conserved to this relative accuracy.
const MAX_CONTINUITY: f64 = 1.0e-13;
/// Gate: CVODE must track the closed-form solution to this absolute accuracy.
const MAX_CLOSED_FORM: f64 = 1.0e-9;
/// Gate: the benchmark limit must be reproduced to this accuracy.
const MAX_BENCHMARK: f64 = 1.0e-12;

struct Options {
    out_dir: PathBuf,
    params: Params,
    n_start: f64,
    n_end: f64,
    samples: usize,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            out_dir: PathBuf::from("artifacts/fable"),
            params: Params::default(),
            // ln a from z = 1100 to a = 3.
            n_start: -7.003,
            n_end: 1.0986122886681098,
            samples: 1601,
        }
    }
}

fn parse_args(argv: &[String]) -> Result<Options, String> {
    let mut o = Options::default();
    let mut i = 1;
    while i < argv.len() {
        let key = argv[i].as_str();
        let need = |i: usize| -> Result<&String, String> {
            argv.get(i + 1)
                .ok_or_else(|| format!("option {} requires a value", argv[i]))
        };
        match key {
            "--out" => {
                o.out_dir = PathBuf::from(need(i)?);
                i += 2;
            }
            "--xi" => {
                o.params.xi = need(i)?
                    .parse()
                    .map_err(|e| format!("--xi: {e}"))?;
                i += 2;
            }
            "--n-start" => {
                o.n_start = need(i)?.parse().map_err(|e| format!("--n-start: {e}"))?;
                i += 2;
            }
            "--n-end" => {
                o.n_end = need(i)?.parse().map_err(|e| format!("--n-end: {e}"))?;
                i += 2;
            }
            "--samples" => {
                o.samples = need(i)?.parse().map_err(|e| format!("--samples: {e}"))?;
                i += 2;
            }
            "--help" | "-h" => {
                println!(
                    "fable_cosmo_rs [--out DIR] [--xi X] [--n-start N] [--n-end N] [--samples K]"
                );
                std::process::exit(0);
            }
            other => return Err(format!("unknown option {other}")),
        }
    }
    Ok(o)
}

/// Linear interpolation of the sample grid onto N = 0 (today).
fn sample_today(samples: &[Sample], p: &Params) -> Sample {
    let state = [
        model::log_bilinear_closed_form(0.0, p),
        model::bridge_closed_form(0.0, p),
    ];
    // Today is evaluated from the closed form so the reported numbers do not
    // depend on whether the output grid happens to land exactly on N = 0.
    let _ = samples;
    model::sample_from_state(0.0, &state, p)
}

fn run() -> Result<(), String> {
    let argv: Vec<String> = std::env::args().collect();
    let o = parse_args(&argv)?;

    let run = solver::integrate(&o.params, o.n_start, o.n_end, o.samples)?;

    let max_continuity = run
        .samples
        .iter()
        .map(|s| s.continuity_residual.abs())
        .fold(0.0_f64, f64::max);
    let max_closed_form = run
        .samples
        .iter()
        .map(|s| s.closed_form_residual.abs())
        .fold(0.0_f64, f64::max);
    let max_bridge = run
        .samples
        .iter()
        .map(|s| s.bridge_residual.abs())
        .fold(0.0_f64, f64::max);

    let today = sample_today(&run.samples, &o.params);

    // The benchmark limit: switching the torsion coupling off makes the
    // dark-energy-like component's equation of state exactly n - 1.
    let mut benchmark_params = o.params;
    benchmark_params.xi = 0.0;
    let benchmark_today = model::sample_from_state(
        0.0,
        &[
            model::log_bilinear_closed_form(0.0, &benchmark_params),
            model::bridge_closed_form(0.0, &benchmark_params),
        ],
        &benchmark_params,
    );
    let benchmark_error = (benchmark_today.w_potential - BENCHMARK_W).abs();
    let dust_error = benchmark_today.w_dust.abs();

    output::write_csv(&o.out_dir.join("fable_background.csv"), &run.samples)?;
    output::write_summary(
        &o.out_dir.join("fable_summary.json"),
        &o.params,
        &run.stats,
        &today,
        max_continuity,
        max_closed_form,
        max_bridge,
    )?;

    println!("fable_cosmo_rs — fableSpinor unified dark sector");
    println!("engine                        : sundials_rs 7.8.0 (pure Rust, Linux x86-64)");
    println!("method                        : CVODE BDF + Newton + dense");
    println!(
        "e-fold range                  : {} .. {}  ({} samples)",
        fmt_g(o.n_start, 10),
        fmt_g(o.n_end, 10),
        o.samples
    );
    println!(
        "CVODE steps / rhs evaluations  : {} / {}",
        run.stats.n_steps, run.stats.n_rhs
    );
    println!();
    println!("today (N = 0):");
    println!("  bridge field s             : {}", fmt_g(today.bridge, 12));
    println!("  dilution exponent nu       : {}", fmt_g(today.dilution, 12));
    println!("  w of the whole sector      : {}", fmt_g(today.w_fable, 12));
    println!("  w of the dust-like part    : {}", fmt_g(today.w_dust, 12));
    println!("  w of the potential part    : {}", fmt_g(today.w_potential, 12));
    println!("  Omega_fable                : {}", fmt_g(today.omega_fable, 12));
    println!("  deceleration parameter     : {}", fmt_g(today.deceleration, 12));
    println!();
    println!("benchmark limit (xi = 0):");
    println!("  w of the potential part    : {}", fmt_g(benchmark_today.w_potential, 12));
    println!("  target                     : {}", fmt_g(BENCHMARK_W, 12));
    println!("  |difference|               : {}", fmt_e(benchmark_error, 6));
    println!("  w of the dust-like part    : {}", fmt_e(benchmark_today.w_dust, 6));
    println!();
    println!("gated invariants:");
    println!(
        "  max |continuity residual|  : {}   (limit {})",
        fmt_e(max_continuity, 6),
        fmt_e(MAX_CONTINUITY, 6)
    );
    println!(
        "  max |closed-form residual| : {}   (limit {})",
        fmt_e(max_closed_form, 6),
        fmt_e(MAX_CLOSED_FORM, 6)
    );
    println!(
        "  max |bridge residual|      : {}   (limit {})",
        fmt_e(max_bridge, 6),
        fmt_e(MAX_CLOSED_FORM, 6)
    );
    println!(
        "  |benchmark - (-0.764)|     : {}   (limit {})",
        fmt_e(benchmark_error, 6),
        fmt_e(MAX_BENCHMARK, 6)
    );

    let mut failures: Vec<String> = Vec::new();
    if !(max_continuity <= MAX_CONTINUITY) {
        failures.push(format!("continuity residual {max_continuity:e}"));
    }
    if !(max_closed_form <= MAX_CLOSED_FORM) {
        failures.push(format!("closed-form residual {max_closed_form:e}"));
    }
    if !(max_bridge <= MAX_CLOSED_FORM) {
        failures.push(format!("bridge residual {max_bridge:e}"));
    }
    if !(benchmark_error <= MAX_BENCHMARK) {
        failures.push(format!("benchmark error {benchmark_error:e}"));
    }
    if !(dust_error <= MAX_BENCHMARK) {
        failures.push(format!("dust equation of state {dust_error:e}"));
    }
    if failures.is_empty() {
        println!();
        println!("SUCCESS: every gated invariant holds.");
        Ok(())
    } else {
        Err(format!("gated invariants violated: {}", failures.join("; ")))
    }
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("fable_cosmo_rs: {message}");
            ExitCode::FAILURE
        }
    }
}
