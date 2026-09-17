"""Background cosmology for the 16-component gpt5_6 spinor condensate."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


SPINOR_COMPONENTS = 16
RELEASE_DATETIME = datetime(2026, 9, 17, tzinfo=timezone.utc)


@dataclass(frozen=True)
class CosmologyParameters:
    """Dimensionless flat-FLRW parameters normalized to today's critical density."""

    omega_m0: float = 0.305
    omega_r0: float = 9.0e-5
    w0: float = -0.861
    wa: float = -0.60
    h0_km_s_mpc: float = 70.0
    spinor_bilinear0: float = 1.0

    @property
    def omega_gpt5_6_0(self) -> float:
        return 1.0 - self.omega_m0 - self.omega_r0

    def validate(self) -> None:
        values = np.array(
            [
                self.omega_m0,
                self.omega_r0,
                self.w0,
                self.wa,
                self.h0_km_s_mpc,
                self.spinor_bilinear0,
            ]
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("All cosmological parameters must be finite.")
        if self.omega_m0 < 0 or self.omega_r0 < 0:
            raise ValueError("Matter and radiation densities must be nonnegative.")
        if self.omega_gpt5_6_0 <= 0:
            raise ValueError("The flat-universe gpt5_6 density must be positive.")
        if self.h0_km_s_mpc <= 0 or self.spinor_bilinear0 <= 0:
            raise ValueError("H0 and the present spinor bilinear must be positive.")


@dataclass(frozen=True)
class ParameterEstimate:
    """A reported central value with asymmetric 68.27% credible errors."""

    value: float
    lower_error: float
    upper_error: float


@dataclass(frozen=True)
class ObservationalBenchmark:
    """Published fit metadata kept separate from the spinor model parameters."""

    label: str
    datasets: tuple[str, ...]
    model: str
    omega_m: ParameterEstimate
    w0: ParameterEstimate
    wa: ParameterEstimate | None
    source: str

    def cosmology_parameters(
        self,
        *,
        omega_r0: float = 9.0e-5,
        h0_km_s_mpc: float = 70.0,
        spinor_bilinear0: float = 1.0,
    ) -> CosmologyParameters:
        """Return central values in the numerical model's parameter format."""

        return CosmologyParameters(
            omega_m0=self.omega_m.value,
            omega_r0=omega_r0,
            w0=self.w0.value,
            wa=0.0 if self.wa is None else self.wa.value,
            h0_km_s_mpc=h0_km_s_mpc,
            spinor_bilinear0=spinor_bilinear0,
        )


UNITE_ONLY_CONSTANT_W_BENCHMARK = ObservationalBenchmark(
    label="Unite-only constant-w benchmark",
    datasets=("Unite",),
    model="Flat-wCDM",
    omega_m=ParameterEstimate(0.197, 0.054, 0.056),
    w0=ParameterEstimate(-0.764, 0.096, 0.078),
    wa=None,
    source="Camilleri et al., arXiv:2609.05053v2, Section 7.0.3",
)

UNITE_BAO_CMB_CPL_BENCHMARK = ObservationalBenchmark(
    label="Unite+BAO+CMB evolving-w benchmark",
    datasets=("Unite", "BAO", "CMB"),
    model="Flat-w0waCDM",
    omega_m=ParameterEstimate(0.305, 0.004, 0.004),
    w0=ParameterEstimate(-0.861, 0.042, 0.044),
    wa=ParameterEstimate(-0.60, 0.19, 0.17),
    source="Camilleri et al., arXiv:2609.05053v2, Section 7.0.4",
)


@dataclass(frozen=True)
class BackgroundSolution:
    e_folds: np.ndarray
    scale_factor: np.ndarray
    redshift: np.ndarray
    hubble_time: np.ndarray
    spinor_bilinear: np.ndarray
    spinor_density: np.ndarray
    spinor_pressure: np.ndarray
    spinor_w: np.ndarray
    hubble_over_h0: np.ndarray
    deceleration: np.ndarray


def dirac_gamma_matrices_16() -> tuple[np.ndarray, ...]:
    """Return Gamma^mu = I_4 (flavor) tensor gamma^mu in signature (+---)."""

    identity2 = np.eye(2, dtype=np.complex128)
    zero2 = np.zeros((2, 2), dtype=np.complex128)
    sigma1 = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    sigma2 = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    sigma3 = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    gamma0 = np.block([[identity2, zero2], [zero2, -identity2]])
    gamma_spatial = tuple(
        np.block([[zero2, sigma], [-sigma, zero2]])
        for sigma in (sigma1, sigma2, sigma3)
    )
    flavor_identity = np.eye(4, dtype=np.complex128)
    return tuple(
        np.kron(flavor_identity, gamma)
        for gamma in (gamma0, *gamma_spatial)
    )


GAMMA_16 = dirac_gamma_matrices_16()


def gpt5_6(scale_factor: float, spinor_bilinear0: float = 1.0) -> np.ndarray:
    """Return an isotropic four-flavor representative of the 16-component field."""

    if scale_factor <= 0 or spinor_bilinear0 <= 0:
        raise ValueError("Scale factor and present bilinear must be positive.")
    scalar_bilinear = spinor_bilinear0 / scale_factor**3
    amplitude = np.sqrt(scalar_bilinear / 4.0)
    field = np.zeros(SPINOR_COMPONENTS, dtype=np.complex128)
    field[[0, 5, 8, 13]] = amplitude
    return field


def spinor_adjoint(field: np.ndarray) -> np.ndarray:
    if field.shape != (SPINOR_COMPONENTS,):
        raise ValueError("gpt5_6 must have exactly 16 components.")
    return field.conjugate() @ GAMMA_16[0]


def scalar_bilinear(field: np.ndarray) -> float:
    value = spinor_adjoint(field) @ field
    return float(np.real_if_close(value))


def vector_current(field: np.ndarray) -> np.ndarray:
    adjoint = spinor_adjoint(field)
    return np.real_if_close(np.array([adjoint @ gamma @ field for gamma in GAMMA_16]))


def cpl_w(scale_factor: np.ndarray | float, parameters: CosmologyParameters) -> np.ndarray:
    scale_factor_array = np.asarray(scale_factor, dtype=float)
    return parameters.w0 + parameters.wa * (1.0 - scale_factor_array)


def bilinear_from_scale_factor(
    scale_factor: np.ndarray | float, parameters: CosmologyParameters
) -> np.ndarray:
    scale_factor_array = np.asarray(scale_factor, dtype=float)
    if np.any(scale_factor_array <= 0):
        raise ValueError("Scale factor must be positive.")
    return parameters.spinor_bilinear0 / scale_factor_array**3


def scale_factor_from_bilinear(
    bilinear: np.ndarray | float, parameters: CosmologyParameters
) -> np.ndarray:
    bilinear_array = np.asarray(bilinear, dtype=float)
    if np.any(bilinear_array <= 0):
        raise ValueError("Spinor bilinear must be positive.")
    return (parameters.spinor_bilinear0 / bilinear_array) ** (1.0 / 3.0)


def spinor_potential(
    bilinear: np.ndarray | float, parameters: CosmologyParameters
) -> np.ndarray:
    """Reconstructed U(S) whose perfect-fluid background exactly follows CPL."""

    parameters.validate()
    bilinear_array = np.asarray(bilinear, dtype=float)
    scale_factor = scale_factor_from_bilinear(bilinear_array, parameters)
    density_ratio = scale_factor ** (
        -3.0 * (1.0 + parameters.w0 + parameters.wa)
    ) * np.exp(-3.0 * parameters.wa * (1.0 - scale_factor))
    return parameters.omega_gpt5_6_0 * density_ratio


def spinor_potential_derivative(
    bilinear: np.ndarray | float, parameters: CosmologyParameters
) -> np.ndarray:
    bilinear_array = np.asarray(bilinear, dtype=float)
    scale_factor = scale_factor_from_bilinear(bilinear_array, parameters)
    return spinor_potential(bilinear_array, parameters) * (
        1.0 + cpl_w(scale_factor, parameters)
    ) / bilinear_array


def spinor_thermodynamics(
    scale_factor: np.ndarray | float, parameters: CosmologyParameters
) -> dict[str, np.ndarray]:
    """Return U, first-order Dirac kinetic density, rho, p, and w."""

    bilinear = bilinear_from_scale_factor(scale_factor, parameters)
    potential = spinor_potential(bilinear, parameters)
    kinetic = bilinear * spinor_potential_derivative(bilinear, parameters)
    pressure = kinetic - potential
    density = potential
    return {
        "bilinear": bilinear,
        "kinetic": kinetic,
        "potential": potential,
        "lagrangian": kinetic - potential,
        "density": density,
        "pressure": pressure,
        "w": pressure / density,
    }


def energy_momentum_mixed(
    scale_factor: float, parameters: CosmologyParameters
) -> np.ndarray:
    """Return T^mu_nu = diag(rho, -p, -p, -p) for homogeneous FLRW."""

    terms = spinor_thermodynamics(scale_factor, parameters)
    density = float(terms["density"])
    pressure = float(terms["pressure"])
    return np.diag([density, -pressure, -pressure, -pressure])


def hubble_squared(
    scale_factor: np.ndarray | float, parameters: CosmologyParameters
) -> np.ndarray:
    scale_factor_array = np.asarray(scale_factor, dtype=float)
    spinor_density = spinor_thermodynamics(scale_factor_array, parameters)["density"]
    return (
        parameters.omega_r0 / scale_factor_array**4
        + parameters.omega_m0 / scale_factor_array**3
        + spinor_density
    )


def solve_background(
    parameters: CosmologyParameters,
    e_folds_min: float = -4.0,
    e_folds_max: float = 1.0,
    points: int = 1201,
) -> BackgroundSolution:
    """Integrate cosmic time, S, and rho as functions of N = ln(a)."""

    parameters.validate()
    if not e_folds_min < 0 < e_folds_max:
        raise ValueError("The integration interval must bracket the present N=0.")
    if points < 101:
        raise ValueError("Use at least 101 points for stable diagnostics and plots.")

    initial = np.array(
        [0.0, parameters.spinor_bilinear0, parameters.omega_gpt5_6_0]
    )

    def right_hand_side(e_folds: float, state: np.ndarray) -> np.ndarray:
        scale_factor = np.exp(e_folds)
        density = state[2]
        expansion_squared = (
            parameters.omega_r0 * np.exp(-4.0 * e_folds)
            + parameters.omega_m0 * np.exp(-3.0 * e_folds)
            + density
        )
        if expansion_squared <= 0 or density <= 0:
            raise RuntimeError("The numerical trajectory left the physical domain.")
        return np.array(
            [
                1.0 / np.sqrt(expansion_squared),
                -3.0 * state[1],
                -3.0 * (1.0 + float(cpl_w(scale_factor, parameters))) * density,
            ]
        )

    past_points = points // 2 + 1
    future_points = points - past_points + 1
    past = solve_ivp(
        right_hand_side,
        (0.0, e_folds_min),
        initial,
        method="DOP853",
        t_eval=np.linspace(0.0, e_folds_min, past_points),
        rtol=1.0e-11,
        atol=1.0e-13,
    )
    future = solve_ivp(
        right_hand_side,
        (0.0, e_folds_max),
        initial,
        method="DOP853",
        t_eval=np.linspace(0.0, e_folds_max, future_points),
        rtol=1.0e-11,
        atol=1.0e-13,
    )
    if not past.success:
        raise RuntimeError(f"Past integration failed: {past.message}")
    if not future.success:
        raise RuntimeError(f"Future integration failed: {future.message}")

    e_folds = np.concatenate((past.t[::-1], future.t[1:]))
    state = np.concatenate((past.y[:, ::-1], future.y[:, 1:]), axis=1)
    scale_factor = np.exp(e_folds)
    spinor_w = cpl_w(scale_factor, parameters)
    spinor_density = state[2]
    spinor_pressure = spinor_w * spinor_density
    radiation_density = parameters.omega_r0 / scale_factor**4
    matter_density = parameters.omega_m0 / scale_factor**3
    total_density = radiation_density + matter_density + spinor_density
    total_pressure = radiation_density / 3.0 + spinor_pressure

    return BackgroundSolution(
        e_folds=e_folds,
        scale_factor=scale_factor,
        redshift=1.0 / scale_factor - 1.0,
        hubble_time=state[0],
        spinor_bilinear=state[1],
        spinor_density=spinor_density,
        spinor_pressure=spinor_pressure,
        spinor_w=spinor_w,
        hubble_over_h0=np.sqrt(total_density),
        deceleration=0.5 * (1.0 + 3.0 * total_pressure / total_density),
    )


def solution_diagnostics(
    solution: BackgroundSolution, parameters: CosmologyParameters
) -> dict[str, float | list[float]]:
    expected_bilinear = bilinear_from_scale_factor(solution.scale_factor, parameters)
    expected_density = spinor_potential(solution.spinor_bilinear, parameters)
    expected_hubble_squared = hubble_squared(solution.scale_factor, parameters)
    bilinear_error = np.max(
        np.abs(solution.spinor_bilinear - expected_bilinear) / expected_bilinear
    )
    density_error = np.max(
        np.abs(solution.spinor_density - expected_density) / expected_density
    )
    friedmann_error = np.max(
        np.abs(solution.hubble_over_h0**2 - expected_hubble_squared)
        / expected_hubble_squared
    )
    phantom_crossing = crossing_scale_factor(solution.scale_factor, solution.spinor_w, -1.0)
    nec_density = solution.spinor_density + solution.spinor_pressure
    nec_crossing = crossing_scale_factor(solution.scale_factor, nec_density, 0.0)
    acceleration_crossings = crossing_scale_factors(
        solution.scale_factor, solution.deceleration, 0.0
    )
    acceleration_redshifts = [1.0 / value - 1.0 for value in acceleration_crossings]
    return {
        "max_relative_bilinear_error": float(bilinear_error),
        "max_relative_density_error": float(density_error),
        "max_relative_friedmann_error": float(friedmann_error),
        "phantom_crossing_a": phantom_crossing,
        "phantom_crossing_z": 1.0 / phantom_crossing - 1.0,
        "nec_crossing_a": nec_crossing,
        "nec_crossing_z": 1.0 / nec_crossing - 1.0,
        "acceleration_crossing_a": acceleration_crossings,
        "acceleration_crossing_z": acceleration_redshifts,
        "w_today": float(cpl_w(1.0, parameters)),
        "w_at_a_0_1": float(cpl_w(0.1, parameters)),
        "w_at_a_2": float(cpl_w(2.0, parameters)),
    }


def crossing_scale_factor(
    scale_factor: np.ndarray, values: np.ndarray, target: float
) -> float:
    crossings = crossing_scale_factors(scale_factor, values, target)
    return crossings[-1] if crossings else float("nan")


def crossing_scale_factors(
    scale_factor: np.ndarray, values: np.ndarray, target: float
) -> list[float]:
    shifted = values - target
    crossing_indices = np.flatnonzero(shifted[:-1] * shifted[1:] <= 0)
    crossings: list[float] = []
    for raw_index in crossing_indices:
        index = int(raw_index)
        x0, x1 = scale_factor[index : index + 2]
        y0, y1 = shifted[index : index + 2]
        crossing = x0 if y1 == y0 else x0 - y0 * (x1 - x0) / (y1 - y0)
        if not crossings or not np.isclose(crossing, crossings[-1]):
            crossings.append(float(crossing))
    return crossings


def save_figures(
    solution: BackgroundSolution,
    parameters: CosmologyParameters,
    output_directory: str | Path,
) -> list[Path]:
    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.22,
            "figure.dpi": 140,
        }
    )
    colors = {"spinor": "#b53a2f", "matter": "#16697a", "radiation": "#d49b26"}
    generated: list[Path] = []

    figure, axis = plt.subplots(figsize=(7.2, 4.4))
    axis.plot(
        solution.scale_factor,
        solution.spinor_w,
        color=colors["spinor"],
        lw=2.4,
        label="Unite+BAO+CMB CPL reconstruction",
    )
    axis.axhline(
        UNITE_ONLY_CONSTANT_W_BENCHMARK.w0.value,
        color=colors["matter"],
        ls="-.",
        lw=1.8,
        label=r"Unite-only constant $w=-0.764$",
    )
    axis.axhline(-1.0, color="#222222", ls="--", lw=1.0, label=r"$\Lambda$CDM / phantom divide")
    axis.axhline(-1.0 / 3.0, color="#666666", ls=":", lw=1.2, label="acceleration threshold")
    axis.axvline(1.0, color="#999999", lw=0.8)
    axis.set_xscale("log")
    axis.set_xlim(0.05, np.exp(1.0))
    axis.set_ylim(-1.6, 0.35)
    axis.set_xlabel("scale factor a")
    axis.set_ylabel("spinor equation of state w")
    axis.legend(frameon=False)
    figure.tight_layout()
    generated.extend(_save_figure_pair(figure, output_path / "gpt5_6_equation_of_state"))

    figure, axis = plt.subplots(figsize=(7.2, 4.4))
    scale_factor = solution.scale_factor
    axis.loglog(scale_factor, solution.spinor_density, color=colors["spinor"], lw=2.4, label="gpt5_6")
    axis.loglog(scale_factor, parameters.omega_m0 / scale_factor**3, color=colors["matter"], lw=1.8, label="matter")
    axis.loglog(scale_factor, parameters.omega_r0 / scale_factor**4, color=colors["radiation"], lw=1.6, label="radiation")
    axis.set_xlabel("scale factor a")
    axis.set_ylabel("density / present critical density")
    axis.legend(frameon=False)
    figure.tight_layout()
    generated.extend(_save_figure_pair(figure, output_path / "gpt5_6_energy_budget"))

    figure, axes = plt.subplots(1, 2, figsize=(9.0, 3.9))
    axes[0].plot(solution.hubble_time, scale_factor, color=colors["matter"], lw=2.2)
    axes[0].axvline(0.0, color="#999999", lw=0.8)
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"dimensionless time $H_0(t-t_0)$")
    axes[0].set_ylabel("scale factor a")
    axes[1].plot(scale_factor, solution.deceleration, color=colors["spinor"], lw=2.2)
    axes[1].axhline(0.0, color="#222222", ls="--", lw=1.0)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("scale factor a")
    axes[1].set_ylabel("deceleration q")
    figure.tight_layout()
    generated.extend(_save_figure_pair(figure, output_path / "gpt5_6_expansion"))

    thermodynamics = spinor_thermodynamics(scale_factor, parameters)
    figure, axis = plt.subplots(figsize=(7.2, 4.4))
    axis.plot(
        thermodynamics["bilinear"],
        thermodynamics["potential"],
        color=colors["spinor"],
        lw=2.4,
        label="potential U(S) = density",
    )
    axis.plot(
        thermodynamics["bilinear"],
        thermodynamics["kinetic"],
        color=colors["matter"],
        lw=1.8,
        label="on-shell Dirac kinetic S U'(S)",
    )
    axis.set_xscale("log")
    axis.set_xlabel(r"bilinear $S=\bar{\Psi}\Psi$ (S0 = 1)")
    axis.set_ylabel("density / present critical density")
    axis.legend(frameon=False)
    figure.tight_layout()
    generated.extend(_save_figure_pair(figure, output_path / "gpt5_6_potential"))
    return generated


def _save_figure_pair(figure: plt.Figure, stem: Path) -> list[Path]:
    png_path = stem.with_suffix(".png")
    pdf_path = stem.with_suffix(".pdf")
    figure.savefig(png_path, bbox_inches="tight")
    figure.savefig(
        pdf_path,
        bbox_inches="tight",
        metadata={
            "CreationDate": RELEASE_DATETIME,
            "ModDate": RELEASE_DATETIME,
        },
    )
    plt.close(figure)
    return [png_path, pdf_path]