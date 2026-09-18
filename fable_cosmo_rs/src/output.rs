//! CSV and summary writers.
//!
//! Every floating-point field goes through the engine's own `fmt_e`/`fmt_g`,
//! which reproduce C's `%.*e` and `%.*g` exactly, so the output is stable
//! across hosts and byte-comparable between runs. Each document is assembled
//! as a `String` first and written once, so a write failure is reported at a
//! single, named site.

use std::fs::{create_dir_all, write};
use std::path::Path;

use sundials_core::sundials_utils::{fmt_e, fmt_g};

use crate::model::{Params, Sample};
use crate::solver::Stats;

/// One CSV column: its header and the field it reads. Header and row are
/// derived from this single table, so they cannot drift apart.
const COLUMNS: [(&str, fn(&Sample) -> f64); 22] = [
    ("e_folds", |s| s.e_folds),
    ("scale_factor", |s| s.scale_factor),
    ("redshift", |s| s.redshift),
    ("bridge", |s| s.bridge),
    ("bridge_h", |s| s.bridge_h),
    ("dilution", |s| s.dilution),
    ("bilinear", |s| s.bilinear),
    ("density_dust", |s| s.density_dust),
    ("density_potential", |s| s.density_potential),
    ("density_fable", |s| s.density_fable),
    ("pressure_dust", |s| s.pressure_dust),
    ("pressure_potential", |s| s.pressure_potential),
    ("pressure_fable", |s| s.pressure_fable),
    ("w_fable", |s| s.w_fable),
    ("w_dust", |s| s.w_dust),
    ("w_potential", |s| s.w_potential),
    ("hubble_over_h0", |s| s.hubble_over_h0),
    ("omega_fable", |s| s.omega_fable),
    ("deceleration", |s| s.deceleration),
    ("continuity_residual", |s| s.continuity_residual),
    ("closed_form_residual", |s| s.closed_form_residual),
    ("bridge_residual", |s| s.bridge_residual),
];

/// The largest residuals over the whole output grid.
#[derive(Clone, Copy, Debug)]
pub struct Invariants {
    pub max_continuity: f64,
    pub max_closed_form: f64,
    pub max_bridge: f64,
}

/// The CSV header line, derived from `COLUMNS`.
pub fn csv_header() -> String {
    COLUMNS
        .iter()
        .map(|(name, _)| *name)
        .collect::<Vec<_>>()
        .join(",")
}

fn csv_row(s: &Sample) -> String {
    // 16 digits after the point is 17 significant digits: round-trip exact.
    COLUMNS
        .iter()
        .map(|(_, field)| fmt_e(field(s), 16))
        .collect::<Vec<_>>()
        .join(",")
}

/// The full background table as text.
pub fn csv_document(samples: &[Sample]) -> String {
    let mut text = String::with_capacity(samples.len() * 24 * 22 + 512);
    text.push_str(&csv_header());
    text.push('\n');
    for s in samples {
        text.push_str(&csv_row(s));
        text.push('\n');
    }
    text
}

/// The machine-readable run summary as JSON text.
pub fn summary_document(p: &Params, stats: &Stats, today: &Sample, inv: &Invariants) -> String {
    let g = |x: f64| fmt_g(x, 17);
    let e = |x: f64| fmt_e(x, 16);
    let mut t = String::with_capacity(1024);
    t.push_str("{\n");
    t.push_str("  \"engine\": \"sundials_rs 7.8.0 (pure Rust, Linux x86-64)\",\n");
    t.push_str("  \"method\": \"CVODE BDF + Newton + dense\",\n");
    t.push_str("  \"parameters\": {\n");
    t.push_str(&format!("    \"omega_r0\": {},\n", g(p.omega_r0)));
    t.push_str(&format!("    \"omega_b0\": {},\n", g(p.omega_b0)));
    t.push_str(&format!("    \"mass_m\": {},\n", g(p.mass_m)));
    t.push_str(&format!("    \"lambda\": {},\n", g(p.lambda)));
    t.push_str(&format!("    \"index_n\": {},\n", g(p.index_n)));
    t.push_str(&format!("    \"xi\": {},\n", g(p.xi)));
    t.push_str(&format!("    \"gamma_flow\": {},\n", g(p.gamma_flow)));
    t.push_str(&format!("    \"bridge_s0\": {},\n", g(p.bridge_s0)));
    t.push_str(&format!("    \"bilinear_s0\": {}\n", g(p.bilinear_s0)));
    t.push_str("  },\n");
    t.push_str("  \"today\": {\n");
    t.push_str(&format!("    \"w_fable\": {},\n", g(today.w_fable)));
    t.push_str(&format!("    \"w_dust\": {},\n", g(today.w_dust)));
    t.push_str(&format!("    \"w_potential\": {},\n", g(today.w_potential)));
    t.push_str(&format!("    \"omega_fable\": {},\n", g(today.omega_fable)));
    t.push_str(&format!("    \"deceleration\": {},\n", g(today.deceleration)));
    t.push_str(&format!("    \"dilution\": {}\n", g(today.dilution)));
    t.push_str("  },\n");
    t.push_str("  \"invariants\": {\n");
    t.push_str(&format!(
        "    \"max_continuity_residual\": {},\n",
        e(inv.max_continuity)
    ));
    t.push_str(&format!(
        "    \"max_closed_form_residual\": {},\n",
        e(inv.max_closed_form)
    ));
    t.push_str(&format!("    \"max_bridge_residual\": {}\n", e(inv.max_bridge)));
    t.push_str("  },\n");
    t.push_str("  \"cvode\": {\n");
    t.push_str(&format!("    \"n_steps\": {},\n", stats.n_steps));
    t.push_str(&format!("    \"n_rhs_evals\": {}\n", stats.n_rhs));
    t.push_str("  }\n");
    t.push_str("}\n");
    t
}

/// Create the parent directory if needed and write the whole document once.
fn write_document(path: &Path, text: &str) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        create_dir_all(parent)
            .map_err(|err| format!("create_dir_all({}) failed: {err}", parent.display()))?;
    }
    write(path, text).map_err(|err| format!("write({}) failed: {err}", path.display()))
}

/// Write the full background table.
pub fn write_csv(path: &Path, samples: &[Sample]) -> Result<(), String> {
    write_document(path, &csv_document(samples))
}

/// Write the machine-readable run summary consumed by the notebook and tests.
pub fn write_summary(
    path: &Path,
    p: &Params,
    stats: &Stats,
    today: &Sample,
    inv: &Invariants,
) -> Result<(), String> {
    write_document(path, &summary_document(p, stats, today, inv))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{closed_form_state, sample_from_state};

    #[test]
    fn header_and_row_have_the_same_arity() {
        let p = Params::default();
        let s = sample_from_state(0.0, &closed_form_state(0.0, &p), &p);
        assert_eq!(csv_header().split(',').count(), COLUMNS.len());
        assert_eq!(csv_row(&s).split(',').count(), COLUMNS.len());
    }

    #[test]
    fn header_starts_and_ends_as_documented() {
        let header = csv_header();
        assert!(header.starts_with("e_folds,scale_factor,redshift,"));
        assert!(header.ends_with(",closed_form_residual,bridge_residual"));
    }

    #[test]
    fn csv_values_round_trip_exactly() {
        let p = Params::default();
        let s = sample_from_state(-3.0, &closed_form_state(-3.0, &p), &p);
        let parsed: Vec<f64> = csv_row(&s)
            .split(',')
            .map(|v| v.parse::<f64>().unwrap())
            .collect();
        for (value, (_, field)) in parsed.iter().zip(COLUMNS.iter()) {
            assert_eq!(*value, field(&s));
        }
    }

    #[test]
    fn summary_is_well_formed_json() {
        // No serde here by design (zero external crates): check the shape the
        // Python side actually parses.
        let p = Params::default();
        let s = sample_from_state(0.0, &closed_form_state(0.0, &p), &p);
        let inv = Invariants {
            max_continuity: 1.0e-15,
            max_closed_form: 3.0e-11,
            max_bridge: 5.0e-11,
        };
        let text = summary_document(&p, &Stats::default(), &s, &inv);
        assert!(text.starts_with("{\n"));
        assert!(text.ends_with("}\n"));
        assert_eq!(text.matches('{').count(), text.matches('}').count());
        for key in ["\"engine\"", "\"parameters\"", "\"today\"", "\"invariants\"", "\"cvode\""] {
            assert!(text.contains(key), "missing {key}");
        }
        assert!(text.contains("\"index_n\": 0.23599999999999999,"));
    }
}
