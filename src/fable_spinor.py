"""fableSpinor: the real Cl(4,4) algebra, the Pin(4,4)-irreducible 16-component
field, and its exactly solvable unified dark-sector background cosmology.

Every algebraic claim in this module is proved symbolically in
``wolfram/fable_spinor.wls``; every numerical claim is cross-checked against the
pure-Rust SUNDIALS 7.8.0 reference integrator in ``fable_cosmo_rs``.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

SPINOR_DIMENSION = 16
TANGENT_DIMENSION = 8

#: The tangent-space metric of the 4+4 split signature.
ETA_4488 = np.diag([1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0])

#: The Unite-only constant-w benchmark (Camilleri et al., arXiv:2609.05053v2,
#: Section 7.0.3). The massless, torsion-free limit of fableSpinor reproduces
#: this value exactly with potential index n = 1 + BENCHMARK_W.
BENCHMARK_W = -0.764


# --------------------------------------------------------------------------- #
# 1. The real Clifford algebra Cl(4,4).
# --------------------------------------------------------------------------- #


def _exterior_basis() -> list[int]:
    """Occupation-number basis of the exterior algebra of R^4.

    Even-parity states occupy positions 0..7 (the type-one block) and
    odd-parity states positions 8..15.
    """
    even = [s for s in range(16) if bin(s).count("1") % 2 == 0]
    odd = [s for s in range(16) if bin(s).count("1") % 2 == 1]
    return even + odd


EXTERIOR_BASIS = _exterior_basis()
_POSITION = {state: index for index, state in enumerate(EXTERIOR_BASIS)}


def _orientation_sign(state: int, index: int) -> int:
    """Koszul sign picked up moving past the already-occupied slots."""
    occupied = sum(1 for j in range(index) if state & (1 << j))
    return -1 if occupied % 2 else 1


def _creation(index: int) -> np.ndarray:
    operator = np.zeros((SPINOR_DIMENSION, SPINOR_DIMENSION))
    for column, state in enumerate(EXTERIOR_BASIS):
        if state & (1 << index):
            continue
        row = _POSITION[state | (1 << index)]
        operator[row, column] = _orientation_sign(state, index)
    return operator


def _annihilation(index: int) -> np.ndarray:
    operator = np.zeros((SPINOR_DIMENSION, SPINOR_DIMENSION))
    for column, state in enumerate(EXTERIOR_BASIS):
        if not state & (1 << index):
            continue
        row = _POSITION[state & ~(1 << index)]
        operator[row, column] = _orientation_sign(state, index)
    return operator


def gamma_matrices_4488() -> np.ndarray:
    """Return the eight real 16x16 generators of Cl(4,4).

    The first four square to +1 and are symmetric; the last four square to -1
    and are antisymmetric. All entries are integers, so the representation is
    real in the strongest possible sense.
    """
    creation = [_creation(i) for i in range(4)]
    annihilation = [_annihilation(i) for i in range(4)]
    plus = [creation[i] + annihilation[i] for i in range(4)]
    minus = [creation[i] - annihilation[i] for i in range(4)]
    return np.array(plus + minus)


GAMMA = gamma_matrices_4488()

#: The chirality operator, the product of all eight generators.
CHIRALITY = np.linalg.multi_dot(list(GAMMA))

#: The spinor metric: the product of the four generators of positive square.
#: It is symmetric, while every ``SPINOR_METRIC @ GAMMA[a]`` is antisymmetric.
SPINOR_METRIC = np.linalg.multi_dot([GAMMA[0], GAMMA[1], GAMMA[2], GAMMA[3]])


def lorentz_generators() -> np.ndarray:
    """The so(4,4) generators ``S^{ab} = (1/4)[gamma^a, gamma^b]``."""
    out = np.zeros(
        (TANGENT_DIMENSION, TANGENT_DIMENSION, SPINOR_DIMENSION, SPINOR_DIMENSION)
    )
    for a in range(TANGENT_DIMENSION):
        for b in range(TANGENT_DIMENSION):
            out[a, b] = 0.25 * (GAMMA[a] @ GAMMA[b] - GAMMA[b] @ GAMMA[a])
    return out


LORENTZ = lorentz_generators()


def commutant_dimension(generators: np.ndarray) -> int:
    """Dimension of the real commutant of a family of 16x16 matrices.

    The commutant is the kernel of the stacked map ``X -> g X - X g``. By real
    Schur / Burnside, dimension 1 means the representation is absolutely
    irreducible over the reals.
    """
    identity = np.eye(SPINOR_DIMENSION)
    blocks = [
        np.kron(g, identity) - np.kron(identity, g.T) for g in generators
    ]
    stacked = np.vstack(blocks)
    rank = np.linalg.matrix_rank(stacked, tol=1e-9)
    return SPINOR_DIMENSION * SPINOR_DIMENSION - rank


def clifford_words(even_only: bool = False) -> np.ndarray:
    """All 256 (or the 128 even) ordered products of distinct generators."""
    words = []
    for mask in range(1 << TANGENT_DIMENSION):
        indices = [i for i in range(TANGENT_DIMENSION) if mask & (1 << i)]
        if even_only and len(indices) % 2:
            continue
        if not indices:
            words.append(np.eye(SPINOR_DIMENSION))
        else:
            words.append(np.linalg.multi_dot([GAMMA[i] for i in indices])
                         if len(indices) > 1 else GAMMA[indices[0]])
    return np.array(words)


def word_span_rank(even_only: bool = False) -> int:
    """Rank of the span of the Clifford words inside the 256-dimensional
    real matrix algebra."""
    words = clifford_words(even_only=even_only)
    return int(np.linalg.matrix_rank(words.reshape(len(words), -1), tol=1e-9))


# --------------------------------------------------------------------------- #
# 2. fableSpinor bilinears.
# --------------------------------------------------------------------------- #


def fable_bilinear(psi: np.ndarray) -> float:
    """The scalar ``S = Psi^T C Psi``, the only nonvanishing background
    bilinear of a commuting real spinor."""
    psi = np.asarray(psi, dtype=float)
    if psi.shape != (SPINOR_DIMENSION,):
        raise ValueError("fableSpinor must have exactly 16 real components.")
    return float(psi @ SPINOR_METRIC @ psi)


def fable_pseudoscalar(psi: np.ndarray) -> float:
    """The pseudoscalar ``P = Psi^T C gamma_9 Psi``."""
    psi = np.asarray(psi, dtype=float)
    if psi.shape != (SPINOR_DIMENSION,):
        raise ValueError("fableSpinor must have exactly 16 real components.")
    return float(psi @ SPINOR_METRIC @ CHIRALITY @ psi)


def fable_current(psi: np.ndarray) -> np.ndarray:
    """The vector current, which vanishes identically for a commuting real
    spinor because every ``C gamma^a`` is antisymmetric."""
    psi = np.asarray(psi, dtype=float)
    if psi.shape != (SPINOR_DIMENSION,):
        raise ValueError("fableSpinor must have exactly 16 real components.")
    return np.array([psi @ SPINOR_METRIC @ GAMMA[a] @ psi for a in range(TANGENT_DIMENSION)])


# --------------------------------------------------------------------------- #
# 3. The background model, mirroring fable_cosmo_rs exactly.
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class FableParameters:
    """Dimensionless flat-FLRW parameters in units of today's critical density."""

    omega_r0: float = 9.0e-5
    omega_b0: float = 0.049
    mass_m: float = 0.256
    lambda_v: float = 1.0 - 9.0e-5 - 0.049 - 0.256
    index_n: float = 1.0 + BENCHMARK_W
    xi: float = 0.15
    gamma_flow: float = 0.35
    bridge_s0: float = 0.5
    bilinear_s0: float = 1.0

    def validate(self) -> None:
        values = np.array(list(asdict(self).values()), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("Every parameter must be finite.")
        if self.lambda_v <= 0:
            raise ValueError("The potential coefficient lambda must be positive.")
        if not 0.0 < self.index_n < 1.0:
            raise ValueError("The potential index n must lie strictly in (0, 1).")
        if not 0.0 < self.bridge_s0 < 1.0:
            raise ValueError("The bridge field today must lie strictly in (0, 1).")
        if self.gamma_flow <= 0:
            raise ValueError("The bridge flow rate must be positive.")
        total = self.omega_r0 + self.omega_b0 + self.mass_m + self.lambda_v
        if abs(total - 1.0) > 1e-12:
            raise ValueError(f"The universe is not flat: densities sum to {total}.")


def bridge_h(s):
    """The smoothstep dialling the gpt-5.6_bridge from Levi-Civita to
    Weitzenboeck."""
    s = np.asarray(s, dtype=float)
    return s * s * (3.0 - 2.0 * s)


def dilution_exponent(s, xi: float):
    """``nu(s) = 3 - xi (1 - h(s))``; exactly 3 at the Weitzenboeck endpoint."""
    return 3.0 - xi * (1.0 - bridge_h(s))


def bridge_closed_form(efolds, p: FableParameters):
    """Closed-form logistic bridge field anchored at ``s(0) = s0``."""
    efolds = np.asarray(efolds, dtype=float)
    growth = np.exp(p.gamma_flow * efolds)
    return p.bridge_s0 * growth / (1.0 - p.bridge_s0 + p.bridge_s0 * growth)


def _quadrature(s, gamma_flow: float):
    return (np.log(s) + s - s * s) / gamma_flow


def log_bilinear_closed_form(efolds, p: FableParameters):
    """Closed-form ``ln S(N)``: the exact solution of the dilution law."""
    efolds = np.asarray(efolds, dtype=float)
    s_now = bridge_closed_form(efolds, p)
    return (
        np.log(p.bilinear_s0)
        - 3.0 * efolds
        + p.xi * (_quadrature(s_now, p.gamma_flow) - _quadrature(p.bridge_s0, p.gamma_flow))
    )


def background_rhs(efolds: float, state: np.ndarray, p: FableParameters) -> np.ndarray:
    """``d/dN [ln S, s]`` — the same two equations the Rust crate integrates."""
    log_s, s = state
    del log_s
    return np.array(
        [-dilution_exponent(s, p.xi), p.gamma_flow * s * (1.0 - s)]
    )


@dataclass(frozen=True)
class FableBackground:
    """A solved background, with every derived quantity on the output grid."""

    e_folds: np.ndarray
    scale_factor: np.ndarray
    redshift: np.ndarray
    bridge: np.ndarray
    dilution: np.ndarray
    bilinear: np.ndarray
    density_dust: np.ndarray
    density_potential: np.ndarray
    density_fable: np.ndarray
    pressure_fable: np.ndarray
    w_fable: np.ndarray
    w_dust: np.ndarray
    w_potential: np.ndarray
    hubble_over_h0: np.ndarray
    omega_fable: np.ndarray
    deceleration: np.ndarray
    continuity_residual: np.ndarray
    closed_form_residual: np.ndarray
    bridge_residual: np.ndarray


def derive(efolds: np.ndarray, log_s: np.ndarray, s: np.ndarray,
           p: FableParameters) -> FableBackground:
    """Build every derived background quantity from the integrated state."""
    bilinear = np.exp(log_s)
    nu = dilution_exponent(s, p.xi)

    density_dust = p.mass_m * bilinear
    density_potential = p.lambda_v * bilinear**p.index_n
    density_fable = density_dust + density_potential

    pressure_dust = (nu / 3.0 - 1.0) * density_dust
    pressure_potential = (p.index_n * nu / 3.0 - 1.0) * density_potential
    pressure_fable = pressure_dust + pressure_potential

    scale_factor = np.exp(efolds)
    density_radiation = p.omega_r0 / scale_factor**4
    density_baryon = p.omega_b0 / scale_factor**3
    density_total = density_radiation + density_baryon + density_fable
    pressure_total = density_radiation / 3.0 + pressure_fable

    d_density = -nu * (density_dust + p.index_n * density_potential)
    continuity = (d_density + 3.0 * (density_fable + pressure_fable)) / density_fable

    return FableBackground(
        e_folds=efolds,
        scale_factor=scale_factor,
        redshift=1.0 / scale_factor - 1.0,
        bridge=s,
        dilution=nu,
        bilinear=bilinear,
        density_dust=density_dust,
        density_potential=density_potential,
        density_fable=density_fable,
        pressure_fable=pressure_fable,
        w_fable=pressure_fable / density_fable,
        w_dust=nu / 3.0 - 1.0,
        w_potential=p.index_n * nu / 3.0 - 1.0,
        hubble_over_h0=np.sqrt(density_total),
        omega_fable=density_fable / density_total,
        deceleration=0.5 * (1.0 + 3.0 * pressure_total / density_total),
        continuity_residual=continuity,
        closed_form_residual=log_s - log_bilinear_closed_form(efolds, p),
        bridge_residual=s - bridge_closed_form(efolds, p),
    )


def solve_background(
    p: FableParameters | None = None,
    n_start: float = -7.003,
    n_end: float = 1.0986122886681098,
    samples: int = 1601,
    rtol: float = 1e-12,
    atol: float = 1e-14,
    method: str = "Radau",
) -> FableBackground:
    """Independent scipy cross-check of the Rust reference integration."""
    p = FableParameters() if p is None else p
    p.validate()
    grid = np.linspace(n_start, n_end, samples)
    initial = np.array(
        [float(log_bilinear_closed_form(n_start, p)), float(bridge_closed_form(n_start, p))]
    )
    solution = solve_ivp(
        background_rhs,
        (n_start, n_end),
        initial,
        t_eval=grid,
        args=(p,),
        method=method,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )
    if not solution.success:
        raise RuntimeError(f"scipy failed to integrate the background: {solution.message}")
    return derive(solution.t, solution.y[0], solution.y[1], p)


# --------------------------------------------------------------------------- #
# 4. The pure-Rust SUNDIALS reference integration.
# --------------------------------------------------------------------------- #

REFERENCE_BINARY = Path("fable_cosmo_rs/target/release/fable_cosmo_rs")


def run_reference(
    out_dir: Path,
    repository_root: Path | None = None,
    binary: Path | None = None,
) -> dict:
    """Run the pure-Rust SUNDIALS CVODE reference integrator.

    Raises if the binary is missing, so a stale or absent build can never be
    silently replaced by the scipy result.
    """
    root = Path.cwd() if repository_root is None else Path(repository_root)
    exe = (root / REFERENCE_BINARY) if binary is None else Path(binary)
    if not exe.is_file():
        raise FileNotFoundError(
            f"The reference binary {exe} does not exist. Build it first with:\n"
            f"    cd {root}/fable_cosmo_rs && cargo build --release"
        )
    target = root / out_dir
    completed = subprocess.run(
        [str(exe), "--out", str(target)],
        cwd=str(root),
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"fable_cosmo_rs exited {completed.returncode}:\n{completed.stdout}\n{completed.stderr}"
        )
    summary = json.loads((target / "fable_summary.json").read_text(encoding="utf-8"))
    summary["stdout"] = completed.stdout
    return summary


def load_reference_csv(path: Path) -> dict[str, np.ndarray]:
    """Load the reference background table written by the Rust binary."""
    raw = np.genfromtxt(path, delimiter=",", names=True)
    return {name: np.asarray(raw[name], dtype=float) for name in raw.dtype.names}


def mixed_error(computed: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Elementwise ``|a - b| / (1 + |b|)``.

    The bilinear ``S`` spans nine decades across the integration range, so a
    bare absolute difference would be meaningless there; for the bounded
    quantities (``w``, the bridge field) the denominator is close to one and
    this reduces to the absolute difference.
    """
    return np.abs(computed - reference) / (1.0 + np.abs(reference))


def compare_to_reference(
    background: FableBackground, reference: dict[str, np.ndarray]
) -> dict[str, float]:
    """Largest mixed absolute/relative disagreement between scipy and Rust."""
    pairs = {
        "e_folds": background.e_folds,
        "bilinear": background.bilinear,
        "bridge": background.bridge,
        "w_fable": background.w_fable,
        "w_potential": background.w_potential,
        "hubble_over_h0": background.hubble_over_h0,
    }
    out = {}
    for name, values in pairs.items():
        if name not in reference:
            raise KeyError(f"The reference table has no column {name!r}.")
        if len(reference[name]) != len(values):
            raise ValueError(
                f"Grid length mismatch for {name!r}: "
                f"{len(reference[name])} reference vs {len(values)} scipy."
            )
        out[name] = float(np.max(mixed_error(values, reference[name])))
    return out
