#!/usr/bin/env python3
"""Solve the gpt5_6 background model and write reproducible artifacts."""

from __future__ import annotations

import json
from pathlib import Path
import sys

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from gpt5_6_cosmology import (
    UNITE_BAO_CMB_CPL_BENCHMARK,
    UNITE_ONLY_CONSTANT_W_BENCHMARK,
    CosmologyParameters,
    ObservationalBenchmark,
    ParameterEstimate,
    save_figures,
    solution_diagnostics,
    solve_background,
    spinor_thermodynamics,
)


def benchmark_record(benchmark: ObservationalBenchmark) -> dict[str, object]:
    """Convert an observational benchmark dataclass into JSON-ready metadata."""

    def estimate_record(
        estimate: ParameterEstimate | None,
    ) -> dict[str, float] | None:
        if estimate is None:
            return None
        return {
            "value": estimate.value,
            "lower_error": estimate.lower_error,
            "upper_error": estimate.upper_error,
        }

    return {
        "label": benchmark.label,
        "datasets": list(benchmark.datasets),
        "model": benchmark.model,
        "omega_m": estimate_record(benchmark.omega_m),
        "w0": estimate_record(benchmark.w0),
        "wa": estimate_record(benchmark.wa),
        "source": benchmark.source,
    }


def main() -> None:
    artifact_directory = REPOSITORY_ROOT / "artifacts"
    figure_directory = artifact_directory / "figures"
    artifact_directory.mkdir(parents=True, exist_ok=True)

    parameters = CosmologyParameters()
    solution = solve_background(parameters)
    diagnostics = solution_diagnostics(solution, parameters)
    today = {
        key: float(value)
        for key, value in spinor_thermodynamics(1.0, parameters).items()
    }
    figures = save_figures(solution, parameters, figure_directory)

    summary = {
        "model": "gpt5_6 16-component nonlinear spinor condensate",
        "components": 16,
        "parameters": {
            "omega_m0": parameters.omega_m0,
            "omega_r0": parameters.omega_r0,
            "omega_gpt5_6_0": parameters.omega_gpt5_6_0,
            "w0": parameters.w0,
            "wa": parameters.wa,
            "h0_km_s_mpc": parameters.h0_km_s_mpc,
            "spinor_bilinear0": parameters.spinor_bilinear0,
        },
        "observational_benchmarks": {
            "unite_only_constant_w": benchmark_record(
                UNITE_ONLY_CONSTANT_W_BENCHMARK
            ),
            "unite_bao_cmb_cpl": benchmark_record(
                UNITE_BAO_CMB_CPL_BENCHMARK
            ),
        },
        "today": today,
        "diagnostics": diagnostics,
        "figures": [str(path.relative_to(REPOSITORY_ROOT)) for path in figures],
        "interpretation": {
            "background_fit_only": True,
            "phantom_crossing_requires_nec_violation": True,
            "perturbative_stability_established": False,
            "dark_matter_included_as_separate_fluid": True,
        },
    }
    summary_path = artifact_directory / "gpt5_6_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()