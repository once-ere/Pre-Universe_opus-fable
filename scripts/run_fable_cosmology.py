#!/usr/bin/env python3
"""Run the fableSpinor background two independent ways and publish the result.

The pure-Rust SUNDIALS 7.8.0 CVODE integration in ``fable_cosmo_rs`` is the
reference. SciPy's Radau integration is an independent cross-check. Both are
compared against the closed-form solution, and every comparison is gated: this
script exits nonzero if any tolerance is exceeded.

Run from the repository root:

    python3 scripts/run_fable_cosmology.py 2>&1 | tee logs/run_fable_cosmology.log
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

import fable_spinor as fs

ARTIFACT_DIR = Path("artifacts/fable")
FIGURE_DIR = ARTIFACT_DIR / "figures"

#: Pinned so regenerated PDF figures are byte-identical between runs.
RELEASE_DATETIME = datetime(2026, 9, 17, tzinfo=timezone.utc)

#: scipy and the Rust reference must agree to this mixed absolute/relative
#: tolerance, |a - b| / (1 + |b|).
MAX_CROSS_CHECK = 1.0e-10
#: both integrations must track the closed form to this absolute tolerance.
MAX_CLOSED_FORM = 1.0e-9
#: the sector must be covariantly conserved to this relative tolerance.
MAX_CONTINUITY = 1.0e-13
#: the benchmark limit must be exact to this tolerance.
MAX_BENCHMARK = 1.0e-12


def _save(figure: plt.Figure, stem: str) -> None:
    target = REPOSITORY_ROOT / FIGURE_DIR
    target.mkdir(parents=True, exist_ok=True)
    figure.savefig(target / f"{stem}.png", dpi=160, bbox_inches="tight")
    # A pinned timestamp is what makes the PDF byte-reproducible.
    figure.savefig(
        target / f"{stem}.pdf",
        dpi=160,
        bbox_inches="tight",
        metadata={"CreationDate": RELEASE_DATETIME, "ModDate": RELEASE_DATETIME},
    )
    plt.close(figure)


def figure_equation_of_state(background: fs.FableBackground, p: fs.FableParameters) -> None:
    """w of the whole sector and of its two components, against redshift."""
    mask = background.redshift >= -0.5
    z = background.redshift[mask]
    figure, axes = plt.subplots(figsize=(7.2, 4.4))
    axes.plot(z, background.w_fable[mask], lw=2.0, label=r"$w$ of the whole fable sector")
    axes.plot(z, background.w_potential[mask], lw=1.8, ls="--",
              label=r"$w$ of the potential part ($n\nu/3-1$)")
    axes.plot(z, background.w_dust[mask], lw=1.8, ls=":",
              label=r"$w$ of the dust-like part ($\nu/3-1$)")
    axes.axhline(fs.BENCHMARK_W, color="crimson", lw=1.2, ls="-.",
                 label=rf"Unite-only benchmark $w={fs.BENCHMARK_W}$")
    axes.set_xscale("symlog", linthresh=1.0)
    axes.set_xlabel("redshift $z$")
    axes.set_ylabel("equation of state $w$")
    axes.set_title(
        rf"fableSpinor equation of state ($n={p.index_n:.3f}$, $\xi={p.xi}$)"
    )
    axes.grid(alpha=0.3)
    axes.legend(fontsize=8, loc="best")
    _save(figure, "fable_equation_of_state")


def figure_dark_sector(background: fs.FableBackground, p: fs.FableParameters) -> None:
    """The dark-matter-like to dark-energy-like handover inside one field."""
    a = background.scale_factor
    figure, axes = plt.subplots(figsize=(7.2, 4.4))
    axes.loglog(a, background.density_dust, lw=2.0,
                label=r"dust-like part $M\,S \propto a^{-\nu}$")
    axes.loglog(a, background.density_potential, lw=2.0, ls="--",
                label=r"dark-energy-like part $\lambda S^{n}$")
    axes.loglog(a, background.density_fable, lw=1.4, ls=":", color="black",
                label="total fable sector")
    axes.loglog(a, p.omega_r0 / a**4, lw=1.0, alpha=0.6, label="radiation")
    axes.loglog(a, p.omega_b0 / a**3, lw=1.0, alpha=0.6, label="baryons")
    axes.axvline(1.0, color="grey", lw=0.8)
    axes.set_xlabel("scale factor $a$")
    axes.set_ylabel(r"density / $\rho_{c,0}$")
    axes.set_title("One field, two behaviours: the unified dark sector")
    axes.grid(alpha=0.3, which="both")
    axes.legend(fontsize=8, loc="best")
    _save(figure, "fable_dark_sector")


def figure_second_mechanism(p: fs.FableParameters) -> None:
    """The torsion mechanism moves w even with the potential index fixed."""
    frozen = fs.FableParameters(**{**vars(p), "xi": 0.0})
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for xi, style in ((0.0, "-"), (0.10, "--"), (0.15, "-."), (0.25, ":")):
        variant = fs.FableParameters(**{**vars(p), "xi": xi})
        solved = fs.solve_background(variant)
        mask = solved.redshift >= -0.5
        axes[0].plot(solved.redshift[mask], solved.w_potential[mask], style,
                     label=rf"$\xi={xi}$")
        axes[1].plot(solved.redshift[mask], solved.w_dust[mask], style,
                     label=rf"$\xi={xi}$")
    for panel, title, target in (
        (axes[0], "potential part", fs.BENCHMARK_W),
        (axes[1], "dust-like part", 0.0),
    ):
        panel.axhline(target, color="crimson", lw=1.0, alpha=0.7)
        panel.set_xscale("symlog", linthresh=1.0)
        panel.set_xlabel("redshift $z$")
        panel.set_ylabel("$w$")
        panel.set_title(f"Second mechanism: {title}")
        panel.grid(alpha=0.3)
        panel.legend(fontsize=8)
    figure.suptitle(
        "Torsion of the gpt-5.6_bridge shifts $w$ at fixed potential index "
        rf"$n={frozen.index_n:.3f}$"
    )
    _save(figure, "fable_second_mechanism")


def figure_expansion(background: fs.FableBackground) -> None:
    """Expansion history and the acceleration it implies."""
    mask = background.redshift >= -0.5
    z = background.redshift[mask]
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    axes[0].plot(z, background.hubble_over_h0[mask], lw=2.0)
    axes[0].set_xscale("symlog", linthresh=1.0)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("redshift $z$")
    axes[0].set_ylabel("$H/H_0$")
    axes[0].set_title("Expansion history")
    axes[0].grid(alpha=0.3, which="both")
    axes[1].plot(z, background.deceleration[mask], lw=2.0)
    axes[1].axhline(0.0, color="crimson", lw=1.0)
    axes[1].set_xscale("symlog", linthresh=1.0)
    axes[1].set_xlabel("redshift $z$")
    axes[1].set_ylabel("$q$")
    axes[1].set_title("Deceleration parameter")
    axes[1].grid(alpha=0.3)
    _save(figure, "fable_expansion")


def figure_cross_check(background: fs.FableBackground,
                       reference: dict[str, np.ndarray]) -> None:
    """SciPy against the pure-Rust SUNDIALS reference, and both against the
    closed form."""
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for name, label in (
        ("bilinear", r"$S$"),
        ("bridge", r"$s$"),
        ("w_fable", r"$w$"),
        ("hubble_over_h0", r"$H/H_0$"),
    ):
        axes[0].semilogy(
            background.e_folds,
            fs.mixed_error(getattr(background, name), reference[name]) + 1e-300,
            lw=1.4,
            label=label,
        )
    axes[0].axhline(MAX_CROSS_CHECK, color="crimson", lw=1.0, ls="--",
                    label="gate")
    axes[0].set_xlabel(r"e-folds $N=\ln a$")
    axes[0].set_ylabel(r"$|a-b|/(1+|b|)$, SciPy Radau vs pure-Rust CVODE")
    axes[0].set_title("Two independent integrators")
    axes[0].grid(alpha=0.3, which="both")
    axes[0].legend(fontsize=8)

    axes[1].semilogy(background.e_folds,
                     np.abs(background.closed_form_residual) + 1e-300,
                     lw=1.4, label=r"SciPy vs closed form, $\ln S$")
    axes[1].semilogy(background.e_folds,
                     np.abs(reference["closed_form_residual"]) + 1e-300,
                     lw=1.4, ls="--", label=r"CVODE vs closed form, $\ln S$")
    axes[1].semilogy(background.e_folds,
                     np.abs(reference["continuity_residual"]) + 1e-300,
                     lw=1.2, ls=":", label="CVODE continuity residual")
    axes[1].axhline(MAX_CLOSED_FORM, color="crimson", lw=1.0, ls="--", label="gate")
    axes[1].set_xlabel(r"e-folds $N=\ln a$")
    axes[1].set_ylabel("absolute residual")
    axes[1].set_title("Against the exact solution")
    axes[1].grid(alpha=0.3, which="both")
    axes[1].legend(fontsize=8)
    _save(figure, "fable_cross_check")


def main() -> int:
    parameters = fs.FableParameters()
    parameters.validate()

    print("fableSpinor background — reference and cross-check")
    print("=" * 66)

    summary = fs.run_reference(ARTIFACT_DIR, repository_root=REPOSITORY_ROOT)
    print(summary["stdout"].rstrip())
    print("=" * 66)

    reference = fs.load_reference_csv(
        REPOSITORY_ROOT / ARTIFACT_DIR / "fable_background.csv"
    )
    background = fs.solve_background(parameters)

    differences = fs.compare_to_reference(background, reference)
    max_cross_check = max(differences.values())
    max_closed_form = max(
        float(np.max(np.abs(background.closed_form_residual))),
        float(np.max(np.abs(reference["closed_form_residual"]))),
    )
    max_continuity = max(
        float(np.max(np.abs(background.continuity_residual))),
        float(np.max(np.abs(reference["continuity_residual"]))),
    )

    benchmark = fs.FableParameters(**{**vars(parameters), "xi": 0.0})
    benchmark_solution = fs.solve_background(benchmark)
    benchmark_error = float(
        np.max(np.abs(benchmark_solution.w_potential - fs.BENCHMARK_W))
    )
    dust_error = float(np.max(np.abs(benchmark_solution.w_dust)))

    print("SciPy Radau vs pure-Rust CVODE, largest |a-b|/(1+|b|):")
    for name, value in sorted(differences.items()):
        print(f"  {name:<18}: {value:.6e}")
    print()
    print(f"max |closed-form residual| : {max_closed_form:.6e}  (gate {MAX_CLOSED_FORM:.1e})")
    print(f"max |continuity residual|  : {max_continuity:.6e}  (gate {MAX_CONTINUITY:.1e})")
    print(f"benchmark |w_pot - (-0.764)|: {benchmark_error:.6e}  (gate {MAX_BENCHMARK:.1e})")
    print(f"benchmark |w_dust|          : {dust_error:.6e}  (gate {MAX_BENCHMARK:.1e})")

    figure_equation_of_state(background, parameters)
    figure_dark_sector(background, parameters)
    figure_second_mechanism(parameters)
    figure_expansion(background)
    figure_cross_check(background, reference)
    print()
    print(f"figures written to {FIGURE_DIR}")

    today = summary["today"]
    published = {
        "engine": summary["engine"],
        "method": summary["method"],
        "cross_check_method": "scipy.integrate.solve_ivp Radau, rtol=1e-12, atol=1e-14",
        "parameters": summary["parameters"],
        "today": today,
        "gates": {
            "max_cross_check_difference": max_cross_check,
            "max_closed_form_residual": max_closed_form,
            "max_continuity_residual": max_continuity,
            "benchmark_w_potential_error": benchmark_error,
            "benchmark_w_dust_error": dust_error,
        },
        "tolerances": {
            "max_cross_check_difference": MAX_CROSS_CHECK,
            "max_closed_form_residual": MAX_CLOSED_FORM,
            "max_continuity_residual": MAX_CONTINUITY,
            "benchmark": MAX_BENCHMARK,
        },
        "cvode": summary["cvode"],
    }
    out = REPOSITORY_ROOT / ARTIFACT_DIR / "fable_published_summary.json"
    out.write_text(json.dumps(published, indent=2) + "\n", encoding="utf-8")
    print(f"summary written to {out.relative_to(REPOSITORY_ROOT)}")

    failures = []
    if not max_cross_check <= MAX_CROSS_CHECK:
        failures.append(f"cross-check {max_cross_check:.3e}")
    if not max_closed_form <= MAX_CLOSED_FORM:
        failures.append(f"closed form {max_closed_form:.3e}")
    if not max_continuity <= MAX_CONTINUITY:
        failures.append(f"continuity {max_continuity:.3e}")
    if not benchmark_error <= MAX_BENCHMARK:
        failures.append(f"benchmark {benchmark_error:.3e}")
    if not dust_error <= MAX_BENCHMARK:
        failures.append(f"dust {dust_error:.3e}")
    if failures:
        print("\nFAILED gates: " + "; ".join(failures))
        return 1
    print("\nSUCCESS: every gate holds.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
