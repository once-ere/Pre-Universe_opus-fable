#!/usr/bin/env python3
"""Gate the fableSpinor deliverables: notebook, published gates, HTML, files.

Called by scripts/verify_all.sh stage 12. Exits nonzero on any failure.
"""

from __future__ import annotations

import json
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent

REQUIRED = [
    "docs/PROVENANCE-FABLESPINOR-ALGEBRA.md",
    "docs/PROVENANCE-GPT-5.6_BRIDGE-DYNAMIC.md",
    "docs/PROVENANCE-FABLE_COSMO_RS.md",
    "docs/PROVENANCE-FABLESPINOR-NOTEBOOK.md",
    "docs/STUDENT-GUIDE-FABLESPINOR.md",
    "docs/fable_spinor.md",
    "docs/fable_spinor.tex",
    "docs/fable_spinor.pdf",
    "wolfram/fable_spinor.wls",
    "wolfram/gpt56_bridge_dynamic.wls",
    "notebooks/fable_spinor_dark_energy.ipynb",
    "fable_cosmo_rs/Cargo.toml",
    "fable_cosmo_rs/src/main.rs",
    "fable_cosmo_rs/src/model.rs",
    "fable_cosmo_rs/src/rhs.rs",
    "fable_cosmo_rs/src/solver.rs",
    "fable_cosmo_rs/src/output.rs",
    "src/fable_spinor.py",
    "tests/test_fable_spinor.py",
    "artifacts/fable/fable_published_summary.json",
    "artifacts/fable/figures/fable_equation_of_state.png",
    "artifacts/fable/figures/fable_dark_sector.png",
    "artifacts/fable/figures/fable_second_mechanism.png",
    "artifacts/fable/figures/fable_expansion.png",
    "artifacts/fable/figures/fable_cross_check.png",
]


def main() -> int:
    notebook = json.loads(
        (REPOSITORY_ROOT / "notebooks/fable_spinor_dark_energy.ipynb").read_text(
            encoding="utf-8"
        )
    )
    code = [c for c in notebook["cells"] if c["cell_type"] == "code"]
    errors = [
        o
        for c in notebook["cells"]
        for o in c.get("outputs", [])
        if o.get("output_type") == "error"
    ]
    unrun = [c for c in code if c.get("execution_count") is None]
    images = [
        o
        for c in notebook["cells"]
        for o in c.get("outputs", [])
        if "image/png" in o.get("data", {})
    ]
    assert not errors, errors
    assert not unrun, "some code cells were never executed"
    assert len(images) >= 2, len(images)

    summary = json.loads(
        (REPOSITORY_ROOT / "artifacts/fable/fable_published_summary.json").read_text(
            encoding="utf-8"
        )
    )
    gates, tolerances = summary["gates"], summary["tolerances"]
    for gate, tolerance in (
        ("max_cross_check_difference", "max_cross_check_difference"),
        ("max_closed_form_residual", "max_closed_form_residual"),
        ("max_continuity_residual", "max_continuity_residual"),
        ("benchmark_w_potential_error", "benchmark"),
        ("benchmark_w_dust_error", "benchmark"),
    ):
        assert gates[gate] <= tolerances[tolerance], (gate, gates[gate])

    html = (REPOSITORY_ROOT / "build/fable_spinor_dark_energy.html").read_text(
        encoding="utf-8"
    )
    assert html.count("<img") >= 2, html.count("<img")

    missing = [p for p in REQUIRED if not (REPOSITORY_ROOT / p).is_file()
               or (REPOSITORY_ROOT / p).stat().st_size == 0]
    assert not missing, f"missing or empty deliverables: {missing}"

    print(f"Notebook: {len(code)} code cells, {len(images)} figures, 0 errors")
    print("fableSpinor gates: passed")
    print(f"fableSpinor deliverables: {len(REQUIRED)} present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
