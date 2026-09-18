"""fableSpinor: the real Cl(4,4) algebra, the Pin(4,4)-irreducible 16-component
field, and its exactly solvable unified dark-sector background cosmology.

Every algebraic claim in this module is proved symbolically in
``wolfram/fable_spinor.wls``; every numerical claim is cross-checked against the
pure-Rust SUNDIALS 7.8.0 reference integrator in ``fable_cosmo_rs``. The
background functions here mirror ``fable_cosmo_rs/src/model.rs`` operation for
operation, including the order of floating-point arithmetic, so the two
implementations can be compared to the last digit.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import asdict, dataclass, replace
from functools import reduce
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
    even = [s for s in range(1 << 4) if s.bit_count() % 2 == 0]
    odd = [s for s in range(1 << 4) if s.bit_count() % 2 == 1]
    return even + odd


EXTERIOR_BASIS = _exterior_basis()
_POSITION = {state: index for index, state in enumerate(EXTERIOR_BASIS)}


def _orientation_sign(state: int, index: int) -> int:
    """Koszul sign picked up moving past the slots below ``index`` that are
    already occupied."""
    occupied_below = (state & ((1 << index) - 1)).bit_count()
    return -1 if occupied_below % 2 else 1


def _ladder(index: int, create: bool) -> np.ndarray:
    """Creation (``create=True``) or annihilation operator for slot ``index``.

    A creation operator acts on states with the slot empty and fills it; an
    annihilation operator acts on states with the slot filled and empties it.
    Both carry the same Koszul sign.
    """
    bit = 1 << index
    operator = np.zeros((SPINOR_DIMENSION, SPINOR_DIMENSION))
    for column, state in enumerate(EXTERIOR_BASIS):
        occupied = bool(state & bit)
        if occupied == create:
            continue
        operator[_POSITION[state ^ bit], column] = _orientation_sign(state, index)
    return operator


def gamma_matrices_4488() -> np.ndarray:
    """Return the eight real 16x16 generators of Cl(4,4).

    ``gamma^i = c_i + a_i`` for the first four (square +1, symmetric) and
    ``gamma^{4+i} = c_i - a_i`` for the last four (square -1, antisymmetric).
    Every entry is an integer, so the representation is real in the strongest
    possible sense; the arrays are float64 only so that downstream linear
    algebra never has to think about integer overflow or dtype promotion.
    """
    creation = [_ladder(i, create=True) for i in range(4)]
    annihilation = [_ladder(i, create=False) for i in range(4)]
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
    products = np.einsum("aij,bjk->abik", GAMMA, GAMMA)
    return 0.25 * (products - products.transpose(1, 0, 2, 3))


LORENTZ = lorentz_generators()


def independent_lorentz_generators() -> np.ndarray:
    """The 28 generators ``S^{ab}`` with ``a < b``: a basis of so(4,4)."""
    upper = np.triu_indices(TANGENT_DIMENSION, k=1)
    return LORENTZ[upper]


def commutant_dimension(generators: np.ndarray) -> int:
    """Dimension of the real commutant of a family of 16x16 matrices.

    The commutant is the kernel of the stacked map ``X -> g X - X g``. By real
    Schur / Burnside, dimension 1 means the representation is absolutely
    irreducible over the reals.
    """
    identity = np.eye(SPINOR_DIMENSION)
    stacked = np.vstack([np.kron(g, identity) - np.kron(identity, g.T) for g in generators])
    rank = np.linalg.matrix_rank(stacked, tol=1e-9)
    return SPINOR_DIMENSION * SPINOR_DIMENSION - rank


def clifford_words(even_only: bool = False) -> np.ndarray:
    """All 256 (or the 128 even) ordered products of distinct generators."""
    identity = np.eye(SPINOR_DIMENSION)
    words = []
    for mask in range(1 << TANGENT_DIMENSION):
        indices = [i for i in range(TANGENT_DIMENSION) if mask & (1 << i)]
        if even_only and len(indices) % 2:
            continue
        words.append(reduce(np.matmul, (GAMMA[i] for i in indices), identity))
    return np.array(words)


def word_span_rank(even_only: bool = False) -> int:
    """Rank of the span of the Clifford words inside the 256-dimensional
    real matrix algebra."""
    words = clifford_words(even_only=even_only)
    return int(np.linalg.matrix_rank(words.reshape(len(words), -1), tol=1e-9))


# --------------------------------------------------------------------------- #
# 2. fableSpinor bilinears.
# --------------------------------------------------------------------------- #


def _as_spinor(psi: np.ndarray) -> np.ndarray:
    psi = np.asarray(psi, dtype=float)
    if psi.shape != (SPINOR_DIMENSION,):
        raise ValueError("fableSpinor must have exactly 16 real components.")
    return psi


def fable_bilinear(psi: np.ndarray) -> float:
    """The scalar ``S = Psi^T C Psi``, the only nonvanishing background
    bilinear of a commuting real spinor."""
    psi = _as_spinor(psi)
    return float(psi @ SPINOR_METRIC @ psi)


def fable_pseudoscalar(psi: np.ndarray) -> float:
    """The pseudoscalar ``P = Psi^T C gamma_9 Psi``."""
    psi = _as_spinor(psi)
    return float(psi @ SPINOR_METRIC @ CHIRALITY @ psi)


def fable_current(psi: np.ndarray) -> np.ndarray:
    """The vector current, which vanishes identically for a commuting real
    spinor because every ``C gamma^a`` is antisymmetric.

    Evaluated as ``((psi @ C) @ gamma^a) @ psi`` per component, in that order:
    the antisymmetric pairs then cancel exactly in floating point and the
    result is a true ``0.0``, which is what the documentation quotes. A single
    ``einsum`` sums in a different order and leaves round-off of order 1e-16.
    """
    psi = _as_spinor(psi)
    return np.array([psi @ SPINOR_METRIC @ GAMMA[a] @ psi for a in range(TANGENT_DIMENSION)])


# --------------------------------------------------------------------------- #
# 3. The background model, mirroring fable_cosmo_rs exactly.
# --------------------------------------------------------------------------- #

# The published scenario's densities today, in units of the critical density.
# Spelled out once so the flat-universe default for lambda is derived rather
# than duplicated; the subtraction order matches fable_cosmo_rs/src/model.rs.
_OMEGA_R0 = 9.0e-5
_OMEGA_B0 = 0.049
_MASS_M = 0.256


@dataclass(frozen=True)
class FableParameters:
    """Dimensionless flat-FLRW parameters in units of today's critical density."""

    omega_r0: float = _OMEGA_R0
    omega_b0: float = _OMEGA_B0
    mass_m: float = _MASS_M
    lambda_v: float = 1.0 - _OMEGA_R0 - _OMEGA_B0 - _MASS_M
    index_n: float = 1.0 + BENCHMARK_W
    xi: float = 0.15
    gamma_flow: float = 0.35
    bridge_s0: float = 0.5
    bilinear_s0: float = 1.0

    def validate(self) -> None:
        """Raise ``ValueError`` on any parameter the model cannot integrate."""
        values = np.array(list(asdict(self).values()), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("Every parameter must be finite.")
        if self.omega_r0 < 0 or self.omega_b0 < 0 or self.mass_m < 0:
            raise ValueError("Densities must be nonnegative.")
        if self.lambda_v <= 0:
            raise ValueError("The potential coefficient lambda must be positive.")
        if not 0.0 < self.index_n < 1.0:
            raise ValueError("The potential index n must lie strictly in (0, 1).")
        if not 0.0 < self.bridge_s0 < 1.0:
            raise ValueError("The bridge field today must lie strictly in (0, 1).")
        if self.gamma_flow <= 0:
            raise ValueError("The bridge flow rate must be positive.")
        if self.bilinear_s0 <= 0:
            raise ValueError("The bilinear today must be positive.")
        total = self.omega_r0 + self.omega_b0 + self.mass_m + self.lambda_v
        if abs(total - 1.0) > 1e-12:
            raise ValueError(f"The universe is not flat: densities sum to {total}.")

    def with_xi(self, xi: float) -> FableParameters:
        """The same model with a different torsion coupling."""
        return replace(self, xi=xi)

    def torsion_free(self) -> FableParameters:
        """The limit in which the dark-energy-like component has exactly
        ``w = n - 1`` and the dust-like component exactly ``w = 0``."""
        return self.with_xi(0.0)


def bridge_h(s: float | np.ndarray) -> np.ndarray:
    """The smoothstep ``h(s) = 3 s^2 - 2 s^3`` dialling the gpt-5.6_bridge from
    Levi-Civita (``h = 0``) to Weitzenboeck (``h = 1``)."""
    s = np.asarray(s, dtype=float)
    return s * s * (3.0 - 2.0 * s)


def dilution_exponent(s: float | np.ndarray, xi: float) -> np.ndarray:
    """``nu(s) = 3 - xi (1 - h(s))``; exactly 3 at the Weitzenboeck endpoint
    and whenever the torsion coupling is switched off."""
    return 3.0 - xi * (1.0 - bridge_h(s))


def bridge_flow(s: float | np.ndarray, gamma_flow: float) -> np.ndarray:
    """The logistic flow of the bridge field, ``ds/dN = gamma s (1 - s)``."""
    s = np.asarray(s, dtype=float)
    return gamma_flow * s * (1.0 - s)


def bridge_closed_form(efolds: float | np.ndarray, p: FableParameters) -> np.ndarray:
    """Closed-form logistic bridge field anchored at ``s(0) = s0``."""
    efolds = np.asarray(efolds, dtype=float)
    growth = np.exp(p.gamma_flow * efolds)
    return p.bridge_s0 * growth / (1.0 - p.bridge_s0 + p.bridge_s0 * growth)


def _quadrature(s: np.ndarray, gamma_flow: float) -> np.ndarray:
    """``(1/gamma)(ln s + s - s^2)``, whose derivative along the logistic flow
    is exactly ``1 - h(s)``."""
    return (np.log(s) + s - s * s) / gamma_flow


def log_bilinear_closed_form(efolds: float | np.ndarray, p: FableParameters) -> np.ndarray:
    """Closed-form ``ln S(N)``: the exact solution of the dilution law."""
    efolds = np.asarray(efolds, dtype=float)
    s_now = bridge_closed_form(efolds, p)
    return (
        np.log(p.bilinear_s0)
        - 3.0 * efolds
        + p.xi * (_quadrature(s_now, p.gamma_flow) - _quadrature(p.bridge_s0, p.gamma_flow))
    )


def closed_form_state(efolds: float, p: FableParameters) -> np.ndarray:
    """The exact state ``[ln S, s]`` at one e-fold, used to seed the integrator."""
    return np.array(
        [float(log_bilinear_closed_form(efolds, p)), float(bridge_closed_form(efolds, p))]
    )


def background_rhs(efolds: float, state: np.ndarray, p: FableParameters) -> np.ndarray:
    """``d/dN [ln S, s]`` — the same two equations the Rust crate integrates.

    Neither equation depends on ``ln S`` itself, so only the bridge field is
    read from the state.
    """
    del efolds
    s = state[1]
    return np.array([-dilution_exponent(s, p.xi), bridge_flow(s, p.gamma_flow)])


@dataclass(frozen=True)
class FableBackground:
    """A solved background, with every derived quantity on the output grid.

    The field names match the columns of the CSV written by ``fable_cosmo_rs``.
    """

    e_folds: np.ndarray
    scale_factor: np.ndarray
    redshift: np.ndarray
    bridge: np.ndarray
    bridge_h: np.ndarray
    dilution: np.ndarray
    bilinear: np.ndarray
    density_dust: np.ndarray
    density_potential: np.ndarray
    density_fable: np.ndarray
    pressure_dust: np.ndarray
    pressure_potential: np.ndarray
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


#: Columns that are physical quantities, compared between the two integrators.
#: The three residual columns are excluded: each is itself round-off noise, so
#: comparing them would measure nothing.
PHYSICAL_COLUMNS: tuple[str, ...] = (
    "e_folds",
    "scale_factor",
    "redshift",
    "bridge",
    "bridge_h",
    "dilution",
    "bilinear",
    "density_dust",
    "density_potential",
    "density_fable",
    "pressure_dust",
    "pressure_potential",
    "pressure_fable",
    "w_fable",
    "w_dust",
    "w_potential",
    "hubble_over_h0",
    "omega_fable",
    "deceleration",
)


def derive(
    efolds: np.ndarray, log_s: np.ndarray, s: np.ndarray, p: FableParameters
) -> FableBackground:
    """Build every derived background quantity from the integrated state.

    Mirrors ``sample_from_state`` in ``fable_cosmo_rs/src/model.rs``.
    """
    bilinear = np.exp(log_s)
    nu = dilution_exponent(s, p.xi)

    density_dust = p.mass_m * bilinear
    density_potential = p.lambda_v * bilinear**p.index_n
    density_fable = density_dust + density_potential

    # p_i = w_i rho_i: the unique pressures for which the sector is covariantly
    # conserved under the dilution law dS/dN = -nu S.
    w_dust = nu / 3.0 - 1.0
    w_potential = p.index_n * nu / 3.0 - 1.0
    pressure_dust = w_dust * density_dust
    pressure_potential = w_potential * density_potential
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
        bridge_h=bridge_h(s),
        dilution=nu,
        bilinear=bilinear,
        density_dust=density_dust,
        density_potential=density_potential,
        density_fable=density_fable,
        pressure_dust=pressure_dust,
        pressure_potential=pressure_potential,
        pressure_fable=pressure_fable,
        w_fable=pressure_fable / density_fable,
        w_dust=w_dust,
        w_potential=w_potential,
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
    """Independent scipy cross-check of the Rust reference integration.

    The defaults match the Rust binary's: e-folds from z = 1100 to a = 3 on
    1601 points, seeded from the exact closed form at ``n_start``.
    """
    p = FableParameters() if p is None else p
    p.validate()
    if not n_end > n_start:
        raise ValueError(f"n_end ({n_end}) must exceed n_start ({n_start}).")
    if samples < 2:
        raise ValueError("At least two output samples are required.")
    grid = np.linspace(n_start, n_end, samples)
    solution = solve_ivp(
        background_rhs,
        (n_start, n_end),
        closed_form_state(n_start, p),
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
REFERENCE_SUMMARY = "fable_summary.json"
REFERENCE_TABLE = "fable_background.csv"


@dataclass(frozen=True)
class ReferenceRun:
    """One completed run of the Rust reference integrator."""

    summary: dict
    stdout: str
    out_dir: Path

    @property
    def table_path(self) -> Path:
        return self.out_dir / REFERENCE_TABLE

    def table(self) -> dict[str, np.ndarray]:
        """The background table the run wrote, one array per column."""
        return load_reference_csv(self.table_path)


def run_reference(
    out_dir: str | Path,
    repository_root: str | Path | None = None,
    binary: str | Path | None = None,
) -> ReferenceRun:
    """Run the pure-Rust SUNDIALS CVODE reference integrator.

    Raises ``FileNotFoundError`` if the binary is missing, so a stale or absent
    build can never be silently replaced by the scipy result, and
    ``RuntimeError`` if the binary's own gated invariants fail.
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
            f"fable_cosmo_rs exited {completed.returncode}:\n"
            f"{completed.stdout}\n{completed.stderr}"
        )
    summary = json.loads((target / REFERENCE_SUMMARY).read_text(encoding="utf-8"))
    return ReferenceRun(summary=summary, stdout=completed.stdout, out_dir=target)


def load_reference_csv(path: str | Path) -> dict[str, np.ndarray]:
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
    background: FableBackground,
    reference: dict[str, np.ndarray],
    columns: tuple[str, ...] = PHYSICAL_COLUMNS,
) -> dict[str, float]:
    """Largest mixed absolute/relative disagreement, per physical column,
    between the scipy background and the Rust reference table."""
    out: dict[str, float] = {}
    for name in columns:
        if name not in reference:
            raise KeyError(f"The reference table has no column {name!r}.")
        values = getattr(background, name)
        if len(reference[name]) != len(values):
            raise ValueError(
                f"Grid length mismatch for {name!r}: "
                f"{len(reference[name])} reference vs {len(values)} scipy."
            )
        out[name] = float(np.max(mixed_error(values, reference[name])))
    return out
