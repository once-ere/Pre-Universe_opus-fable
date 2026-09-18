#!/usr/bin/env python3
"""Build the executable teaching notebook for the gpt5_6 cosmology model."""

from __future__ import annotations

from pathlib import Path

import nbformat


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATH = REPOSITORY_ROOT / "notebooks" / "gpt5_6_cosmology.ipynb"


def markdown(cell_id: str, source: str) -> nbformat.NotebookNode:
    return nbformat.v4.new_markdown_cell(source.strip(), id=cell_id)


def code(cell_id: str, source: str) -> nbformat.NotebookNode:
    return nbformat.v4.new_code_cell(source.strip(), id=cell_id)


def build_notebook() -> int:
    cells = [
        markdown(
            "title",
            r"""
# $\mathrm{gpt5\_6}$: a 16-component nonlinear spinor cosmology

This notebook defines a four-flavor Dirac field, derives its homogeneous
stress-energy tensor, reconstructs a time-varying dark-energy background, and
checks the result numerically. It is written for a reader who has not solved a
cosmological ODE before.

**Claim boundary.** This is an effective background reconstruction. It is not a
fit to supernova light curves, a new irreducible Lorentz representation, or a
proof of perturbative or quantum stability.
""",
        ),
        markdown(
            "roadmap",
            r"""
## What will be built

1. Define $\mathrm{gpt5\_6}\in\mathbb C^{16}$ as four ordinary Dirac flavors.
2. Write the covariant action, field equation, and symmetric stress tensor.
3. Derive kinetic density, potential, pressure, density, and $w=p/\rho$.
4. Reconstruct a potential that realizes the CPL history
   $w(a)=w_0+w_a(1-a)$.
5. Integrate the FLRW background in $N=\ln a$, check analytic invariants, and
   locate the phantom/NEC and acceleration transitions.
6. Compare with the separate Unite-only constant-$w$ result and discuss dark
   matter, dark energy, and the model's limitations.
""",
        ),
        markdown(
            "setup",
            r"""
## Setup from an empty Python environment

Open a terminal in the repository root. On Linux or macOS, create an isolated
environment and install the pinned packages:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
```

On Windows PowerShell, replace `.venv/bin/python` with
`.venv\Scripts\python.exe`. Start Jupyter with
`.venv/bin/python -m jupyter lab`, open this file, and choose **Run All**.

The notebook imports the tested implementation from `src/`; equations are not
silently reimplemented in notebook cells.
""",
        ),
        code(
            "imports",
            r"""
from pathlib import Path
import sys

import numpy as np
from IPython.display import Markdown, display

REPOSITORY_ROOT = Path.cwd().resolve()
if not (REPOSITORY_ROOT / "src").is_dir():
    REPOSITORY_ROOT = REPOSITORY_ROOT.parent
assert (REPOSITORY_ROOT / "src" / "gpt5_6_cosmology.py").is_file()
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from gpt5_6_cosmology import (
    GAMMA_16,
    UNITE_BAO_CMB_CPL_BENCHMARK,
    UNITE_ONLY_CONSTANT_W_BENCHMARK,
    CosmologyParameters,
    bilinear_from_scale_factor,
    cpl_w,
    crossing_scale_factors,
    energy_momentum_mixed,
    gpt5_6,
    hubble_squared,
    save_figures,
    scalar_bilinear,
    solution_diagnostics,
    solve_background,
    spinor_potential,
    spinor_thermodynamics,
    vector_current,
)

parameters = UNITE_BAO_CMB_CPL_BENCHMARK.cosmology_parameters()
constant_w_parameters = UNITE_ONLY_CONSTANT_W_BENCHMARK.cosmology_parameters()
parameters
""",
        ),
        markdown(
            "observations",
            r"""
## Observational inputs: do not mix these two fits

Camilleri et al. report two relevant, but statistically different, results:

| data | model | $\Omega_m$ | equation of state |
|---|---|---:|---:|
| Unite supernovae only | flat $w$CDM | $0.197^{+0.056}_{-0.054}$ | $w=-0.764^{+0.078}_{-0.096}$ |
| Unite + BAO + CMB | flat $w_0w_a$CDM | $0.305\pm0.004$ | $w_0=-0.861^{+0.044}_{-0.042}$, $w_a=-0.60^{+0.17}_{-0.19}$ |

The first row assumes constant $w$. The second uses the Chevallier-Polarski-
Linder (CPL) form and supplies the central curve reconstructed here. These
numbers cannot be combined as though they came from one likelihood. The paper
also reports that Bayesian evidence for evolving dark energy is weak even where
frequentist summaries reach roughly $3.1$--$3.3\sigma$; neither statement is a
discovery of this spinor model.
""",
        ),
        code(
            "benchmark-check",
            r"""
assert UNITE_ONLY_CONSTANT_W_BENCHMARK.datasets == ("Unite",)
assert UNITE_ONLY_CONSTANT_W_BENCHMARK.wa is None
assert UNITE_BAO_CMB_CPL_BENCHMARK.datasets == ("Unite", "BAO", "CMB")

sample_a = np.array([0.1, 1.0, 2.0])
display(Markdown(
    "**CPL central values:** "
    + ", ".join(
        f"$w({a:g})={w:.3f}$"
        for a, w in zip(sample_a, cpl_w(sample_a, parameters), strict=True)
    )
))
display(Markdown(
    f"**Separate constant-$w$ benchmark:** $w={constant_w_parameters.w0:.3f}$ "
    "at every scale factor."
))
""",
        ),
        markdown(
            "field-definition",
            r"""
## 1. The 16-component field

In $3+1$ dimensions an ordinary Dirac spinor has four complex components. We
define the new field as a reducible four-flavor multiplet

$$
\mathrm{gpt5\_6}\equiv\Psi=
(\psi_1,\psi_2,\psi_3,\psi_4)^T\in\mathbb C^{16},\qquad
\Gamma^\mu=I_4\otimes\gamma^\mu,
$$

so $\{\Gamma^\mu,\Gamma^\nu\}=2\eta^{\mu\nu}I_{16}$ for signature
$(+---)$. The Dirac adjoint is $\bar\Psi=\Psi^\dagger\Gamma^0$ and the
Lorentz scalar driving the condensate is $S=\bar\Psi\Psi$.

For the background calculation, `gpt5_6(a)` returns one convenient homogeneous
representative with $S=S_0a^{-3}$ and vanishing spatial vector current. The
observable background depends on $S$, not on this representative's arbitrary
flavor orientation.
""",
        ),
        code(
            "field-check",
            r"""
metric = np.diag([1.0, -1.0, -1.0, -1.0])
identity16 = np.eye(16)
for mu in range(4):
    for nu in range(4):
        anticommutator = GAMMA_16[mu] @ GAMMA_16[nu] + GAMMA_16[nu] @ GAMMA_16[mu]
        np.testing.assert_allclose(anticommutator, 2.0 * metric[mu, nu] * identity16)

field_today = gpt5_6(1.0, parameters.spinor_bilinear0)
assert field_today.shape == (16,)
np.testing.assert_allclose(scalar_bilinear(field_today), parameters.spinor_bilinear0)
np.testing.assert_allclose(vector_current(field_today)[1:], 0.0, atol=1.0e-14)

{
    "components": field_today.size,
    "S_today": scalar_bilinear(field_today),
    "spatial_current": vector_current(field_today)[1:].tolist(),
    "Clifford_checks": "16 of 16 passed",
}
""",
        ),
        markdown(
            "action",
            r"""
## 2. Action, field equation, and stress-energy tensor

Use units $c=\hbar=1$ and the torsion-free spinor covariant derivative
$\nabla_\mu$. The complete background action is

$$
\mathcal S=\int d^4x\sqrt{-g}\left[
\frac{R}{2\kappa}+\mathcal L_\Psi\right]+\mathcal S_m+\mathcal S_r,
$$

with the Hermitian first-order spinor Lagrangian

$$
\mathcal L_\Psi=\frac{i}{2}\left[
\bar\Psi\Gamma^\mu\nabla_\mu\Psi-
(\nabla_\mu\bar\Psi)\Gamma^\mu\Psi\right]-U(S).
$$

Independent variation of $\bar\Psi$ gives

$$i\Gamma^\mu\nabla_\mu\Psi-U_{,S}\Psi=0,$$

and variation with respect to the metric/tetrad gives the symmetric tensor

$$
T_{\mu\nu}=\frac{i}{4}\left[
\bar\Psi\Gamma_\mu\nabla_\nu\Psi+
\bar\Psi\Gamma_\nu\nabla_\mu\Psi-
(\nabla_\mu\bar\Psi)\Gamma_\nu\Psi-
(\nabla_\nu\bar\Psi)\Gamma_\mu\Psi\right]
-g_{\mu\nu}\mathcal L_\Psi.
$$

The adjoint variation gives the conjugate Dirac equation. Together they imply
covariant conservation $\nabla_\mu T^{\mu\nu}=0$ on shell.
""",
        ),
        markdown(
            "thermodynamics",
            r"""
## 3. Homogeneous FLRW thermodynamics

For $ds^2=dt^2-a(t)^2d\mathbf x^2$ and a homogeneous isotropic scalar
condensate, the Dirac equations reduce to

$$\dot S+3HS=0,\qquad S(a)=S_0a^{-3}.$$

The mixed stress tensor has perfect-fluid form
$T^\mu{}_{\nu}=\operatorname{diag}(\rho,-p,-p,-p)$, with

$$
K_D\equiv\frac{i}{2}\left[
\bar\Psi\Gamma^\mu\nabla_\mu\Psi-
(\nabla_\mu\bar\Psi)\Gamma^\mu\Psi\right]=S U_{,S},
$$

$$
\rho_\Psi=U(S),\qquad
p_\Psi=K_D-U=SU_{,S}-U,\qquad
w_\Psi=-1+\frac{S U_{,S}}{U}.
$$

Thus the on-shell Lagrangian density equals the pressure. $K_D$ is a
first-order Dirac kinetic density, not the nonnegative
$\dot\phi^2/2$ of canonical scalar quintessence. In particular, this effective
spinor background can cross $w=-1$; canonical single-field quintessence cannot.
""",
        ),
        code(
            "stress-check",
            r"""
today = spinor_thermodynamics(1.0, parameters)
stress_today = energy_momentum_mixed(1.0, parameters)
expected_stress = np.diag([
    today["density"],
    -today["pressure"],
    -today["pressure"],
    -today["pressure"],
])
np.testing.assert_allclose(stress_today, expected_stress)
np.testing.assert_allclose(today["pressure"], today["kinetic"] - today["potential"])
np.testing.assert_allclose(today["w"], today["pressure"] / today["density"])
stress_today
""",
        ),
        markdown(
            "reconstruction",
            r"""
## 4. Reconstructing the CPL potential

The CPL equation of state is

$$w(a)=w_0+w_a(1-a).$$

Energy conservation, $d\rho/d\ln a=-3(1+w)\rho$, integrates to

$$
\rho_\Psi(a)=\rho_{\Psi0}\,
a^{-3(1+w_0+w_a)}\exp[-3w_a(1-a)].
$$

Substitute $a=(S_0/S)^{1/3}$ to obtain the nonlinear potential

$$
U(S)=\rho_{\Psi0}
\left(\frac{S}{S_0}\right)^{1+w_0+w_a}
\exp\left[-3w_a\left(1-
\left(\frac{S_0}{S}\right)^{1/3}\right)\right].
$$

Direct differentiation gives $S U_{,S}/U=1+w(a)$, which proves that
$p/\rho$ reproduces the prescribed CPL curve. This is a reconstruction: the
observations select a background history, and that history is used to define
$U(S)$. They do not independently select this microscopic spinor theory.
""",
        ),
        code(
            "reconstruction-check",
            r"""
scale_factor_check = np.geomspace(0.05, 2.5, 301)
thermodynamics_check = spinor_thermodynamics(scale_factor_check, parameters)
np.testing.assert_allclose(
    thermodynamics_check["w"],
    cpl_w(scale_factor_check, parameters),
    rtol=2.0e-14,
)
np.testing.assert_allclose(
    thermodynamics_check["density"],
    spinor_potential(bilinear_from_scale_factor(scale_factor_check, parameters), parameters),
)
"Exact CPL reconstruction passed on 301 scale factors."
""",
        ),
        markdown(
            "ode",
            r"""
## 5. Numerical initial-value problem

Choose $N=\ln a$ as the independent variable and
$\tau=H_0(t-t_0)$. Flatness fixes
$\Omega_{\Psi0}=1-\Omega_{m0}-\Omega_{r0}$. The three integrated variables
obey

$$
\frac{d\tau}{dN}=\frac{1}{E},\qquad
\frac{dS}{dN}=-3S,\qquad
\frac{d\rho_\Psi}{dN}=-3[1+w(e^N)]\rho_\Psi,
$$

where

$$
E^2(N)\equiv\frac{H^2}{H_0^2}=
\Omega_{r0}e^{-4N}+\Omega_{m0}e^{-3N}+\rho_\Psi/\rho_{c0}.
$$

Initial conditions at today ($N=0$, $a=1$) are
$\tau=0$, $S=S_0=1$, and
$\rho_\Psi/\rho_{c0}=\Omega_{\Psi0}=0.69491$. `solve_ivp` uses the
eighth-order DOP853 method, relative tolerance $10^{-11}$, and absolute
tolerance $10^{-13}$. It integrates from today into the past and future
separately, then joins the arrays in increasing $a$.
""",
        ),
        code(
            "solve",
            r"""
solution = solve_background(parameters, e_folds_min=-4.0, e_folds_max=1.0, points=1201)
diagnostics = solution_diagnostics(solution, parameters)

assert diagnostics["max_relative_bilinear_error"] < 2.0e-9
assert diagnostics["max_relative_density_error"] < 2.0e-9
assert diagnostics["max_relative_friedmann_error"] < 2.0e-9
np.testing.assert_allclose(
    solution.hubble_over_h0**2,
    hubble_squared(solution.scale_factor, parameters),
    rtol=2.0e-9,
)
diagnostics
""",
        ),
        markdown(
            "results-intro",
            r"""
## 6. Numerical results and transitions

The phantom crossing is defined by $w=-1$. Because
$\rho_\Psi+p_\Psi=SU_{,S}=\rho_\Psi(1+w)$, it is exactly the
background null-energy-condition (NEC) boundary. The acceleration boundaries
are roots of

$$q=-\frac{\ddot a}{aH^2}=\frac12\left(1+3\frac{p_{\rm tot}}{\rho_{\rm tot}}\right).$$

The first $q=0$ crossing is in the observed past. The second is a formal future
consequence of extending the empirical CPL ansatz beyond its fitted domain.
""",
        ),
        code(
            "results-table",
            r"""
acceleration_a = diagnostics["acceleration_crossing_a"]
acceleration_z = diagnostics["acceleration_crossing_z"]
rows = [
    ("present spinor density rho_Psi / rho_c0", today["density"]),
    ("present spinor pressure p_Psi / rho_c0", today["pressure"]),
    ("present Dirac kinetic K_D / rho_c0", today["kinetic"]),
    ("w(0.1)", diagnostics["w_at_a_0_1"]),
    ("w(1)", diagnostics["w_today"]),
    ("w(2)", diagnostics["w_at_a_2"]),
    ("phantom/NEC crossing a", diagnostics["phantom_crossing_a"]),
    ("phantom/NEC crossing z", diagnostics["phantom_crossing_z"]),
    ("acceleration begins a", acceleration_a[0]),
    ("acceleration begins z", acceleration_z[0]),
    ("formal future acceleration end a", acceleration_a[1]),
    ("formal future acceleration end z", acceleration_z[1]),
    ("max relative S invariant error", diagnostics["max_relative_bilinear_error"]),
    ("max relative density invariant error", diagnostics["max_relative_density_error"]),
    ("max relative Friedmann closure error", diagnostics["max_relative_friedmann_error"]),
]
display(Markdown("\n".join([
    "| Quantity | Value |",
    "|---|---:|",
    *[f"| {label} | {float(value):.10g} |" for label, value in rows],
])))
""",
        ),
        code(
            "transition-checks",
            r"""
analytic_crossing = 1.0 + (1.0 + parameters.w0) / parameters.wa
np.testing.assert_allclose(diagnostics["phantom_crossing_a"], analytic_crossing, rtol=2.0e-6)
np.testing.assert_allclose(diagnostics["nec_crossing_a"], analytic_crossing, rtol=2.0e-6)

crossing_terms = spinor_thermodynamics(analytic_crossing, parameters)
np.testing.assert_allclose(crossing_terms["density"] + crossing_terms["pressure"], 0.0, atol=1.0e-14)

constant_grid = np.geomspace(0.05, np.exp(1.0), 301)
constant_values = cpl_w(constant_grid, constant_w_parameters)
assert crossing_scale_factors(constant_grid, constant_values, -1.0) == []
"Transition and constant-w no-crossing checks passed."
""",
        ),
        code(
            "figures",
            r"""
figure_directory = REPOSITORY_ROOT / "artifacts" / "figures"
generated_figures = save_figures(solution, parameters, figure_directory)

for name, caption in (
    ("gpt5_6_equation_of_state.png", "CPL, Unite-only constant-w, and Lambda-CDM comparison"),
    ("gpt5_6_energy_budget.png", "Radiation, matter, and reconstructed spinor density"),
    ("gpt5_6_expansion.png", "Scale factor versus time and the two q=0 crossings"),
    ("gpt5_6_potential.png", "Reconstructed U(S) and on-shell Dirac kinetic density"),
):
    display(Markdown(
        f"### {caption}\n\n![{caption}](../artifacts/figures/{name})"
    ))

[str(path.relative_to(REPOSITORY_ROOT)) for path in generated_figures]
""",
        ),
        markdown(
            "interpretation",
            r"""
## 7. What the model says about dark energy and dark matter

**Dark energy.** In this parameter regime, the spinor condensate is a second
field-theoretic mechanism for the requested evolving dark-energy background.
It has negative pressure today, $w(1)=-0.861$, and drives acceleration when
combined with matter and radiation. Its central CPL trajectory crosses the
phantom divide at $a\simeq0.76833$ ($z\simeq0.30152$), simultaneously changing
the sign of $\rho+p$. Calling the whole trajectory canonical quintessence would
therefore be incorrect.

**Dark matter.** A nonlinear spinor has a dust limit: if $U(S)=mS$, then
$SU_{,S}=U$, hence $p=0$, $w=0$, and $\rho\propto a^{-3}$. That is a mathematical
connection to cold dark matter. It is not the regime used in this run: the
Friedmann equation already contains a separate matter fluid
$\Omega_{m0}a^{-3}$. Reinterpreting this same condensate as all of the dark
matter would double-count the matter density unless a new unified model were
fitted and tested.
""",
        ),
        markdown(
            "limitations",
            r"""
## 8. What remains unproved

- The 16-component object is four Dirac flavors, not a new irreducible spinor.
- The homogeneous field is an effective condensate or mean field; this notebook
  does not quantize it or explain its microscopic origin.
- The potential is reconstructed from a chosen CPL central curve. It is not
  unique, and no observational likelihood or covariance is evaluated here.
- A background NEC crossing is not by itself a stability analysis. Ghost,
  gradient, causal, fermionic perturbation, radiative, and strong-coupling
  questions remain open.
- Structure growth, CMB perturbations, supernova calibration systematics, BAO,
  and local gravity constraints are not calculated.
- CPL is an empirical low-redshift parameterization. Its $a>1$ behavior,
  including the future $q=0$ crossing, is an extrapolation rather than a
  prediction supported by the cited data.

Consequently the calculation demonstrates mathematical background
reconstructibility, not that $\mathrm{gpt5\_6}$ is realized in nature.
""",
        ),
        markdown(
            "rustsolveit",
            r"""
## 9. Why rustSolveIt/SUNDIALS is not used in this notebook

The requested Linux rustSolveIt repository was evaluated. It vendors a
pure-Rust SUNDIALS 7.8.0 implementation and provides its own mechanics-oriented
notebook language and Jupyter wrapper, but it does not expose a drop-in Python
ODE API for this FLRW state vector. Using it here would require cloning a large
external workspace, writing and maintaining a separate Rust executable, and
bridging its output back to Python. The present smooth, nonstiff three-state
problem also has closed-form invariants that DOP853 reproduces at about
$10^{-11}$.

For those reasons the reproducible primary path uses SciPy. The decision is an
interface and scope choice, not a claim that SUNDIALS is unsuitable. A future
stiff perturbation system could justify a dedicated CVODE Rust executable.
""",
        ),
        markdown(
            "reproduction",
            r"""
## 10. Reproduce every deliverable

From the repository root:

```bash
.venv/bin/python -m pytest tests/test_gpt5_6_cosmology.py -q
.venv/bin/python scripts/run_cosmology.py
.venv/bin/python scripts/build_cosmology_notebook.py
.venv/bin/python -m jupyter nbconvert --to notebook --execute --inplace \
    notebooks/gpt5_6_cosmology.ipynb --ExecutePreprocessor.timeout=180 \
    --ExecutePreprocessor.record_timing=False \
    --KernelManager.transport_encryption=auto
.venv/bin/python -m jupyter nbconvert --to html --output-dir build \
  notebooks/gpt5_6_cosmology.ipynb
```

Expected result: 10 tests pass; the notebook finishes without an error output;
all three relative invariant errors are below $2\times10^{-9}$; four PNG/PDF
figure pairs and `artifacts/gpt5_6_summary.json` exist. See `PROVENANCE.md` for
source hashes, exact validation commands, and the report build.
""",
        ),
        markdown(
            "references",
            r"""
## References

1. R. Camilleri et al., *Supernovae Unite: Combining Pantheon+ and
   DES-SN5YR*, arXiv:2609.05053v2, especially Sections 7.0.3--7.0.4.
2. Y.-F. Cai and J. Wang, *Dark Energy Model with Spinor Matter and Its
   Quintom Scenario*, Class. Quantum Grav. **25** (2008) 165014,
   arXiv:0806.3890.
3. M. Chevallier and D. Polarski, Int. J. Mod. Phys. D **10** (2001) 213.
4. E. V. Linder, Phys. Rev. D **68** (2003) 083504, arXiv:astro-ph/0208512.
5. once-ere, *rustSolveIt on pure-Rust SUNDIALS 7.8.0 for Ubuntu Linux*,
   <https://github.com/once-ere/rustSolveIt_linux_SUNDIALS_7_8_0>.
""",
        ),
    ]

    notebook = nbformat.v4.new_notebook(
        cells=cells,
        metadata={
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python", "version": "3"},
        },
    )
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(notebook, NOTEBOOK_PATH)
    return len(cells)


def main() -> int:
    cell_count = build_notebook()
    print(f"cells: {cell_count}")
    print(f"notebook: {NOTEBOOK_PATH.relative_to(REPOSITORY_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())