#![forbid(unsafe_code)]
#![deny(warnings)]

//! fable_cosmo_rs — the reference background integration of the fableSpinor
//! unified dark sector.
//!
//! Usage:
//!   cargo run --release -- [--out DIR] [--xi X] [--n-start N] [--n-end N] [--samples K]
//!
//! Exit status is nonzero if any gated invariant is violated. Every number is
//! printed through the engine's own `fmt_e`/`fmt_g`, never `{:e}`.

mod model;
mod output;
mod rhs;
mod solver;

use std::path::PathBuf;
use std::process::ExitCode;

use sundials_core::sundials_utils::{fmt_e, fmt_g};

use model::{closed_form_state, sample_from_state, Params, Sample, BENCHMARK_W};
use output::Invariants;

/// Gate: the sector must be covariantly conserved to this relative accuracy.
const MAX_CONTINUITY: f64 = 1.0e-13;
/// Gate: CVODE must track the closed-form solution to this absolute accuracy.
const MAX_CLOSED_FORM: f64 = 1.0e-9;
/// Gate: the benchmark limit must be reproduced to this accuracy.
const MAX_BENCHMARK: f64 = 1.0e-12;

const USAGE: &str =
    "fable_cosmo_rs [--out DIR] [--xi X] [--n-start N] [--n-end N] [--samples K]";

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

/// What `parse_args` asks the caller to do.
enum Command {
    Run(Options),
    Help,
}

fn parse_value<T: std::str::FromStr>(argv: &[String], i: usize) -> Result<T, String>
where
    T::Err: std::fmt::Display,
{
    let key = &argv[i];
    let raw = argv
        .get(i + 1)
        .ok_or_else(|| format!("option {key} requires a value"))?;
    raw.parse::<T>().map_err(|e| format!("{key}: {e}"))
}

fn parse_args(argv: &[String]) -> Result<Command, String> {
    let mut o = Options::default();
    let mut i = 1;
    while i < argv.len() {
        match argv[i].as_str() {
            "--out" => o.out_dir = PathBuf::from(parse_value::<String>(argv, i)?),
            "--xi" => o.params.xi = parse_value(argv, i)?,
            "--n-start" => o.n_start = parse_value(argv, i)?,
            "--n-end" => o.n_end = parse_value(argv, i)?,
            "--samples" => o.samples = parse_value(argv, i)?,
            "--help" | "-h" => return Ok(Command::Help),
            other => return Err(format!("unknown option {other}\nusage: {USAGE}")),
        }
        i += 2;
    }
    Ok(Command::Run(o))
}

/// Largest absolute value of one residual field over the whole grid.
fn max_abs(samples: &[Sample], field: fn(&Sample) -> f64) -> f64 {
    samples
        .iter()
        .map(|s| field(s).abs())
        .fold(0.0_f64, f64::max)
}

/// "Today" evaluated from the closed form, so the reported numbers do not
/// depend on whether the output grid happens to land exactly on N = 0.
fn today(p: &Params) -> Sample {
    sample_from_state(0.0, &closed_form_state(0.0, p), p)
}

fn run(o: &Options) -> Result<(), String> {
    let run = solver::integrate(&o.params, o.n_start, o.n_end, o.samples)?;

    let inv = Invariants {
        max_continuity: max_abs(&run.samples, |s| s.continuity_residual),
        max_closed_form: max_abs(&run.samples, |s| s.closed_form_residual),
        max_bridge: max_abs(&run.samples, |s| s.bridge_residual),
    };
    let now = today(&o.params);

    // The benchmark limit: switching the torsion coupling off makes the
    // dark-energy-like component's equation of state exactly n - 1.
    let benchmark = today(&o.params.torsion_free());
    let benchmark_error = (benchmark.w_potential - BENCHMARK_W).abs();
    let dust_error = benchmark.w_dust.abs();

    output::write_csv(&o.out_dir.join("fable_background.csv"), &run.samples)?;
    output::write_summary(
        &o.out_dir.join("fable_summary.json"),
        &o.params,
        &run.stats,
        &now,
        &inv,
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
    println!("  bridge field s             : {}", fmt_g(now.bridge, 12));
    println!("  dilution exponent nu       : {}", fmt_g(now.dilution, 12));
    println!("  w of the whole sector      : {}", fmt_g(now.w_fable, 12));
    println!("  w of the dust-like part    : {}", fmt_g(now.w_dust, 12));
    println!("  w of the potential part    : {}", fmt_g(now.w_potential, 12));
    println!("  Omega_fable                : {}", fmt_g(now.omega_fable, 12));
    println!("  deceleration parameter     : {}", fmt_g(now.deceleration, 12));
    println!();
    println!("benchmark limit (xi = 0):");
    println!("  w of the potential part    : {}", fmt_g(benchmark.w_potential, 12));
    println!("  target                     : {}", fmt_g(BENCHMARK_W, 12));
    println!("  |difference|               : {}", fmt_e(benchmark_error, 6));
    println!("  w of the dust-like part    : {}", fmt_e(benchmark.w_dust, 6));
    println!();
    println!("gated invariants:");
    let gate = |label: &str, value: f64, limit: f64| {
        println!(
            "  {label:<26} : {}   (limit {})",
            fmt_e(value, 6),
            fmt_e(limit, 6)
        );
    };
    gate("max |continuity residual|", inv.max_continuity, MAX_CONTINUITY);
    gate("max |closed-form residual|", inv.max_closed_form, MAX_CLOSED_FORM);
    gate("max |bridge residual|", inv.max_bridge, MAX_CLOSED_FORM);
    gate("|benchmark - (-0.764)|", benchmark_error, MAX_BENCHMARK);

    // `!(x <= limit)` rather than `x > limit` so a NaN residual fails the gate.
    let checks = [
        ("continuity residual", inv.max_continuity, MAX_CONTINUITY),
        ("closed-form residual", inv.max_closed_form, MAX_CLOSED_FORM),
        ("bridge residual", inv.max_bridge, MAX_CLOSED_FORM),
        ("benchmark error", benchmark_error, MAX_BENCHMARK),
        ("dust equation of state", dust_error, MAX_BENCHMARK),
    ];
    let failures: Vec<String> = checks
        .iter()
        .filter(|(_, value, limit)| !(value <= limit))
        .map(|(name, value, _)| format!("{name} {}", fmt_e(*value, 6)))
        .collect();

    if failures.is_empty() {
        println!();
        println!("SUCCESS: every gated invariant holds.");
        Ok(())
    } else {
        Err(format!("gated invariants violated: {}", failures.join("; ")))
    }
}

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().collect();
    let outcome = match parse_args(&argv) {
        Ok(Command::Help) => {
            println!("{USAGE}");
            return ExitCode::SUCCESS;
        }
        Ok(Command::Run(options)) => run(&options),
        Err(message) => Err(message),
    };
    match outcome {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("fable_cosmo_rs: {message}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(args: &[&str]) -> Vec<String> {
        std::iter::once("fable_cosmo_rs")
            .chain(args.iter().copied())
            .map(String::from)
            .collect()
    }

    #[test]
    fn defaults_when_no_arguments_are_given() {
        match parse_args(&argv(&[])).unwrap() {
            Command::Run(o) => {
                assert_eq!(o.samples, 1601);
                assert_eq!(o.params, Params::default());
                assert_eq!(o.out_dir, PathBuf::from("artifacts/fable"));
            }
            Command::Help => panic!("expected a run"),
        }
    }

    #[test]
    fn parses_every_option() {
        let cmd = parse_args(&argv(&[
            "--out", "x/y", "--xi", "0.25", "--n-start", "-12", "--n-end", "6", "--samples", "11",
        ]))
        .unwrap();
        match cmd {
            Command::Run(o) => {
                assert_eq!(o.out_dir, PathBuf::from("x/y"));
                assert_eq!(o.params.xi, 0.25);
                assert_eq!(o.n_start, -12.0);
                assert_eq!(o.n_end, 6.0);
                assert_eq!(o.samples, 11);
            }
            Command::Help => panic!("expected a run"),
        }
    }

    #[test]
    fn help_is_a_command_not_an_exit() {
        assert!(matches!(parse_args(&argv(&["--help"])).unwrap(), Command::Help));
        assert!(matches!(parse_args(&argv(&["-h"])).unwrap(), Command::Help));
    }

    #[test]
    fn rejects_missing_values_and_unknown_options() {
        assert!(parse_args(&argv(&["--xi"])).is_err());
        assert!(parse_args(&argv(&["--xi", "not-a-number"])).is_err());
        assert!(parse_args(&argv(&["--bogus"])).is_err());
    }

    #[test]
    fn nan_fails_a_gate() {
        // The gate must use !(x <= limit), so a NaN residual is a failure.
        let value = f64::NAN;
        assert!(!(value <= MAX_CONTINUITY));
    }
}
