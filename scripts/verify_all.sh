#!/usr/bin/env bash

set -euo pipefail

repository_root="$(realpath "$(dirname "${BASH_SOURCE[0]}")/..")"
cd "$repository_root"

python_executable="${PYTHON:-.venv/bin/python}"
if [[ ! -x "$python_executable" ]]; then
  printf 'Missing Python environment: %s\n' "$python_executable" >&2
  printf 'Create it with: python3 -m venv .venv && %s -m pip install -r requirements.txt\n' \
    "$python_executable" >&2
  exit 2
fi

mkdir -p build logs

printf '\n[1/12] Python tests\n'
"$python_executable" -m pytest -q | tee logs/verify-python.log

printf '\n[2/12] Numerical artifacts\n'
"$python_executable" scripts/run_cosmology.py > logs/verify-cosmology.json

printf '\n[3/12] Jupyter notebook build and execution\n'
"$python_executable" scripts/build_cosmology_notebook.py
"$python_executable" -m jupyter nbconvert --to notebook --execute --inplace \
  notebooks/gpt5_6_cosmology.ipynb --ExecutePreprocessor.timeout=180 \
  --ExecutePreprocessor.record_timing=False \
  2>&1 | tee logs/verify-notebook-execution.log
"$python_executable" -m jupyter nbconvert --to html --output-dir build \
  notebooks/gpt5_6_cosmology.ipynb \
  2>&1 | tee logs/verify-notebook-html.log

printf '\n[4/12] Structured cosmology checks\n'
"$python_executable" - <<'PY'
import json
from pathlib import Path

notebook = json.loads(Path("notebooks/gpt5_6_cosmology.ipynb").read_text())
errors = [
    output
    for cell in notebook["cells"]
    for output in cell.get("outputs", [])
    if output.get("output_type") == "error"
]
assert not errors, errors

summary = json.loads(Path("artifacts/gpt5_6_summary.json").read_text())
diagnostics = summary["diagnostics"]
for key in (
    "max_relative_bilinear_error",
    "max_relative_density_error",
    "max_relative_friedmann_error",
):
    assert diagnostics[key] < 2.0e-9, (key, diagnostics[key])
assert summary["observational_benchmarks"]["unite_only_constant_w"]["datasets"] == ["Unite"]
assert summary["observational_benchmarks"]["unite_bao_cmb_cpl"]["datasets"] == ["Unite", "BAO", "CMB"]

html = Path("build/gpt5_6_cosmology.html").read_text()
assert html.count("<img") >= 4
for caption in (
    "CPL, Unite-only constant-w, and Lambda-CDM comparison",
    "Radiation, matter, and reconstructed spinor density",
    "Scale factor versus time and the two q=0 crossings",
    "Reconstructed U(S) and on-shell Dirac kinetic density",
):
    assert f'alt="{caption}"' in html

print("Notebook errors: 0")
print("Invariant thresholds: passed")
print("HTML figure alt text: passed")
PY

printf '\n[5/12] Markdown and Python hygiene\n'
git diff --check
"$python_executable" -m compileall -q src scripts tests

printf '\n[6/12] LaTeX/PDF report\n'
bash scripts/build_documentation.sh > logs/verify-documentation.log 2>&1
for report in gpt5_6_cosmology gpt5_6_dark_sector_relationships; do
  if grep -En 'Warning|Error|Undefined|undefined|Overfull|Underfull' \
    "build/${report}.log"; then
    printf 'The final TeX log for %s contains a diagnostic.\n' "$report" >&2
    exit 1
  fi
  printf '%s:\n' "$report"
  pdfinfo "docs/${report}.pdf" | grep -E '^(Pages|File size|PDF version):'
done

printf '\n[7/12] Wolfram source and generated notebook\n'
"$python_executable" scripts/build_gpt56_notebook.py
wolframscript -file wolfram/gpt56_bridge.wls \
  2>&1 | tee logs/verify-wolfram-source.log
wolframscript -file scripts/run_gpt56_notebook.wls \
  2>&1 | tee logs/verify-wolfram-notebook.log
grep -q 'Tests succeeded: 36' logs/verify-wolfram-source.log
grep -q 'Tests failed: 0' logs/verify-wolfram-source.log
grep -q 'Notebook cells with messages: 0' logs/verify-wolfram-notebook.log
grep -q 'Notebook tests succeeded: 36' logs/verify-wolfram-notebook.log
grep -q 'Notebook tests failed: 0' logs/verify-wolfram-notebook.log

printf '\n[8/12] Required deliverables\n'
for path in \
  README.md \
  PROVENANCE.md \
  docs/gpt5_6_cosmology.md \
  docs/gpt5_6_cosmology.tex \
  docs/gpt5_6_cosmology.pdf \
  docs/gpt5_6_dark_sector_relationships.md \
  docs/gpt5_6_dark_sector_relationships.tex \
  docs/gpt5_6_dark_sector_relationships.pdf \
  notebooks/gpt5_6_cosmology.ipynb \
  notebooks/gpt-5.6_bridge.nb \
  artifacts/SHA256SUMS \
  artifacts/gpt5_6_summary.json \
  artifacts/figures/gpt5_6_equation_of_state.png \
  artifacts/figures/gpt5_6_energy_budget.png \
  artifacts/figures/gpt5_6_expansion.png \
  artifacts/figures/gpt5_6_potential.png; do
  test -s "$path"
done
sha256sum -c artifacts/SHA256SUMS > logs/verify-checksums.log
printf 'Release checksums: %s files passed\n' \
  "$(wc -l < logs/verify-checksums.log)"

printf '\n[9/12] fableSpinor symbolic proofs (Wolfram)\n'
wolframscript -file wolfram/fable_spinor.wls \
  2>&1 | tee logs/verify-fable-spinor.log
wolframscript -file wolfram/gpt56_bridge_dynamic.wls \
  2>&1 | tee logs/verify-gpt56-bridge-dynamic.log
grep -q 'Tests succeeded: 30' logs/verify-fable-spinor.log
grep -q 'Tests failed: 0' logs/verify-fable-spinor.log
grep -q 'Pin(4,4) commutant dimension          : 1' logs/verify-fable-spinor.log
grep -q 'Spin(4,4) commutant dimension         : 2' logs/verify-fable-spinor.log
grep -q 'Succeeded    : 17' logs/verify-gpt56-bridge-dynamic.log
grep -q 'Failed       : 0' logs/verify-gpt56-bridge-dynamic.log

printf '\n[10/12] fable_cosmo_rs: pure-Rust SUNDIALS reference integration\n'
if [[ ! -e vendor/sundials_rs/crates/cvode_rs/Cargo.toml ]]; then
  printf 'The SUNDIALS submodule is missing. Run:\n' >&2
  printf '  git submodule update --init --recursive\n' >&2
  exit 2
fi
( cd fable_cosmo_rs && cargo build --release ) 2>&1 \
  | tee logs/verify-fable-cosmo-build.log
if grep -Eq '^(warning|error)' logs/verify-fable-cosmo-build.log; then
  printf 'The Rust build emitted a diagnostic.\n' >&2
  exit 1
fi
./fable_cosmo_rs/target/release/fable_cosmo_rs --out artifacts/fable \
  2>&1 | tee logs/verify-fable-cosmo-run.log
grep -q 'SUCCESS: every gated invariant holds.' logs/verify-fable-cosmo-run.log

printf '\n[11/12] Three-way cross-check and figures\n'
"$python_executable" scripts/run_fable_cosmology.py \
  2>&1 | tee logs/verify-fable-cross-check.log
grep -q 'SUCCESS: every gate holds.' logs/verify-fable-cross-check.log

printf '\n[12/12] fableSpinor notebook, report and deliverables\n'
"$python_executable" scripts/build_fable_notebook.py
"$python_executable" -m jupyter nbconvert --to notebook --execute --inplace \
  notebooks/fable_spinor_dark_energy.ipynb --ExecutePreprocessor.timeout=900 \
  --ExecutePreprocessor.record_timing=False \
  2>&1 | tee logs/verify-fable-notebook-execution.log
"$python_executable" -m jupyter nbconvert --to html --output-dir build \
  notebooks/fable_spinor_dark_energy.ipynb \
  2>&1 | tee logs/verify-fable-notebook-html.log
"$python_executable" scripts/check_fable_deliverables.py

if grep -En 'Warning|Error|Undefined|undefined|Overfull|Underfull' \
  build/fable_spinor.log; then
  printf 'The final TeX log for fable_spinor contains a diagnostic.\n' >&2
  exit 1
fi
printf 'fable_spinor:\n'
pdfinfo docs/fable_spinor.pdf | grep -E '^(Pages|File size|PDF version):'

printf '\nAll repository verification gates passed.\n'
