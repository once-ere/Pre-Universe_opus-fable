#!/usr/bin/env python3
"""Determinism gate: execute every notebook and artifact twice, compare bytes.

Run A and run B are complete, independent regenerations. Every tracked artifact
must be byte-identical between them, or the artifact is not reproducible and
must not be pinned in artifacts/SHA256SUMS.

Run from the repository root:

    .venv/bin/python scripts/determinism_gate.py 2>&1 | tee logs/determinism_gate.log
"""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
PYTHON = REPOSITORY_ROOT / ".venv/bin/python"

WATCHED = [
    "artifacts/gpt5_6_summary.json",
    "artifacts/figures/gpt5_6_energy_budget.png",
    "artifacts/figures/gpt5_6_energy_budget.pdf",
    "artifacts/figures/gpt5_6_equation_of_state.png",
    "artifacts/figures/gpt5_6_equation_of_state.pdf",
    "artifacts/figures/gpt5_6_expansion.png",
    "artifacts/figures/gpt5_6_expansion.pdf",
    "artifacts/figures/gpt5_6_potential.png",
    "artifacts/figures/gpt5_6_potential.pdf",
    "notebooks/gpt5_6_cosmology.ipynb",
    "notebooks/gpt-5.6_bridge.nb",
    "artifacts/fable/fable_summary.json",
    "artifacts/fable/fable_published_summary.json",
    "artifacts/fable/figures/fable_equation_of_state.png",
    "artifacts/fable/figures/fable_equation_of_state.pdf",
    "artifacts/fable/figures/fable_dark_sector.png",
    "artifacts/fable/figures/fable_dark_sector.pdf",
    "artifacts/fable/figures/fable_second_mechanism.png",
    "artifacts/fable/figures/fable_second_mechanism.pdf",
    "artifacts/fable/figures/fable_expansion.png",
    "artifacts/fable/figures/fable_expansion.pdf",
    "artifacts/fable/figures/fable_cross_check.png",
    "artifacts/fable/figures/fable_cross_check.pdf",
    "notebooks/fable_spinor_dark_energy.ipynb",
    "docs/gpt5_6_cosmology.pdf",
    "docs/gpt5_6_dark_sector_relationships.pdf",
    "docs/fable_spinor.pdf",
]


def run(command: list[str], label: str) -> None:
    print(f"    {label}", flush=True)
    completed = subprocess.run(command, cwd=REPOSITORY_ROOT, capture_output=True, text=True)
    if completed.returncode != 0:
        print(completed.stdout[-4000:])
        print(completed.stderr[-4000:], file=sys.stderr)
        raise SystemExit(f"{label} failed with status {completed.returncode}")


def regenerate() -> None:
    """One complete regeneration, in the same order and with the same tools
    as scripts/verify_all.sh, so what is compared is exactly what is committed."""
    py = str(PYTHON)
    nbconvert = [py, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute",
                 "--inplace", "--ExecutePreprocessor.record_timing=False"]
    run([py, "scripts/run_cosmology.py"], "gpt5_6 numerical artifacts")
    run([py, "scripts/build_cosmology_notebook.py"], "gpt5_6 notebook build")
    run(nbconvert + ["--ExecutePreprocessor.timeout=600", "notebooks/gpt5_6_cosmology.ipynb"],
        "gpt5_6 notebook execution")
    run(["./fable_cosmo_rs/target/release/fable_cosmo_rs", "--out", "artifacts/fable"],
        "fable_cosmo_rs reference integration")
    run([py, "scripts/run_fable_cosmology.py"], "fable three-way cross-check and figures")
    run([py, "scripts/build_fable_notebook.py"], "fable notebook build")
    run(nbconvert + ["--ExecutePreprocessor.timeout=900", "notebooks/fable_spinor_dark_energy.ipynb"],
        "fable notebook execution")
    # The committed notebooks are the normalized ones, so normalize before hashing.
    run([py, "scripts/normalize_notebooks.py"], "notebook normalization")
    run([py, "scripts/build_gpt56_notebook.py"], "gpt-5.6 bridge notebook build")
    # execute_gpt56_notebook.wls writes the outputs back into the .nb; the older
    # run_gpt56_notebook.wls evaluates but discards them and would leave the
    # watched file unexecuted.
    run(["wolframscript", "-file", "scripts/execute_gpt56_notebook.wls"],
        "gpt-5.6 bridge notebook execution in place")
    run(["bash", "scripts/build_documentation.sh"], "LaTeX reports")


def digests() -> dict[str, str]:
    out = {}
    for name in WATCHED:
        path = REPOSITORY_ROOT / name
        if not path.is_file():
            raise SystemExit(f"watched artifact is missing: {name}")
        out[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return out


def main() -> int:
    if not PYTHON.is_file():
        raise SystemExit(
            f"missing {PYTHON}; create it with `uv venv .venv && uv pip install "
            "--python .venv/bin/python -r requirements.txt` (or python3 -m venv .venv)"
        )
    print("Run A: regenerating every artifact")
    regenerate()
    first = digests()

    print("Run B: regenerating every artifact again")
    regenerate()
    second = digests()

    drifted = [name for name in WATCHED if first[name] != second[name]]
    width = max(len(name) for name in WATCHED)
    for name in WATCHED:
        mark = "DRIFT" if name in drifted else "same "
        print(f"  {mark}  {name:<{width}}  {second[name][:16]}")

    print()
    print(f"watched artifacts : {len(WATCHED)}")
    print(f"byte-identical    : {len(WATCHED) - len(drifted)}")
    print(f"drifted           : {len(drifted)}")
    if drifted:
        print("\nFAILED: these artifacts are not byte-reproducible and must not be pinned:")
        for name in drifted:
            print(f"  {name}")
        return 1
    print("\nSUCCESS: every watched artifact is byte-identical across two full executions.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
