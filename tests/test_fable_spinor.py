"""Tests for the fableSpinor algebra and its background cosmology.

Run from the repository root (pytest.ini puts src/ on the import path):

    python3 -m pytest tests/test_fable_spinor.py -q
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import fable_spinor as fs

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent


# --------------------------------------------------------------------------- #
# The real Clifford algebra.
# --------------------------------------------------------------------------- #


def test_gammas_are_real_integer_matrices():
    assert fs.GAMMA.shape == (8, 16, 16)
    assert np.array_equal(fs.GAMMA, np.round(fs.GAMMA))
    assert np.isrealobj(fs.GAMMA)


def test_clifford_relation():
    identity = np.eye(fs.SPINOR_DIMENSION)
    for a in range(fs.TANGENT_DIMENSION):
        for b in range(fs.TANGENT_DIMENSION):
            anticommutator = fs.GAMMA[a] @ fs.GAMMA[b] + fs.GAMMA[b] @ fs.GAMMA[a]
            assert np.array_equal(anticommutator, 2 * fs.ETA_4488[a, b] * identity)


def test_signature_is_four_plus_four():
    squares = [float(np.trace(g @ g)) / fs.SPINOR_DIMENSION for g in fs.GAMMA]
    assert squares[:4] == [1.0, 1.0, 1.0, 1.0]
    assert squares[4:] == [-1.0, -1.0, -1.0, -1.0]


def test_spinor_metric_is_symmetric_and_gamma_products_are_antisymmetric():
    assert np.array_equal(fs.SPINOR_METRIC, fs.SPINOR_METRIC.T)
    for g in fs.GAMMA:
        product = fs.SPINOR_METRIC @ g
        assert np.array_equal(product, -product.T)


def test_bilinear_is_invariant_under_the_spin_generators():
    for a in range(fs.TANGENT_DIMENSION):
        for b in range(fs.TANGENT_DIMENSION):
            generator = fs.LORENTZ[a, b]
            residual = generator.T @ fs.SPINOR_METRIC + fs.SPINOR_METRIC @ generator
            assert np.allclose(residual, 0.0, atol=1e-12)


def test_chirality_splits_the_spin_representation_but_not_the_pin_one():
    for a in range(fs.TANGENT_DIMENSION):
        assert np.allclose(fs.CHIRALITY @ fs.GAMMA[a] + fs.GAMMA[a] @ fs.CHIRALITY, 0.0)
        for b in range(fs.TANGENT_DIMENSION):
            generator = fs.LORENTZ[a, b]
            assert np.allclose(
                fs.CHIRALITY @ generator - generator @ fs.CHIRALITY, 0.0, atol=1e-12
            )


# --------------------------------------------------------------------------- #
# Irreducibility under the universal cover of the real O(4,4).
# --------------------------------------------------------------------------- #


def test_pin_representation_is_absolutely_irreducible():
    """Commutant dimension 1 is exactly real Schur's criterion."""
    assert fs.commutant_dimension(fs.GAMMA) == 1


def test_spin_representation_is_reducible():
    """Under the connected group the 16 splits as 8 + 8, so the commutant is
    spanned by the identity and the chirality operator."""
    generators = fs.independent_lorentz_generators()
    assert generators.shape == (28, 16, 16)
    assert fs.commutant_dimension(generators) == 2


def test_clifford_words_span_the_full_real_matrix_algebra():
    assert fs.word_span_rank() == 256
    assert fs.word_span_rank(even_only=True) == 128


# --------------------------------------------------------------------------- #
# Bilinears of a commuting real spinor.
# --------------------------------------------------------------------------- #


def test_scalar_and_pseudoscalar_are_generically_nonzero():
    rng = np.random.default_rng(20260917)
    psi = rng.standard_normal(fs.SPINOR_DIMENSION)
    assert abs(fs.fable_bilinear(psi)) > 1e-6
    assert abs(fs.fable_pseudoscalar(psi)) > 1e-6


def test_vector_current_vanishes_identically():
    rng = np.random.default_rng(1)
    for _ in range(16):
        psi = rng.standard_normal(fs.SPINOR_DIMENSION)
        assert np.allclose(fs.fable_current(psi), 0.0, atol=1e-12)


def test_wrong_component_count_is_rejected():
    with pytest.raises(ValueError):
        fs.fable_bilinear(np.zeros(15))


# --------------------------------------------------------------------------- #
# The background model.
# --------------------------------------------------------------------------- #


def test_bridge_endpoints():
    assert fs.bridge_h(0.0) == 0.0
    assert fs.bridge_h(1.0) == 1.0
    assert fs.dilution_exponent(1.0, xi=0.37) == pytest.approx(3.0)
    assert fs.dilution_exponent(0.0, xi=0.37) == pytest.approx(3.0 - 0.37)


def test_default_parameters_describe_a_flat_universe():
    p = fs.FableParameters()
    p.validate()
    total = p.omega_r0 + p.omega_b0 + p.mass_m + p.lambda_v
    assert total == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize(
    "override",
    [
        {"index_n": 1.5},
        {"index_n": 0.0},
        {"bridge_s0": 0.0},
        {"gamma_flow": -1.0},
        {"mass_m": 0.5},
    ],
)
def test_invalid_parameters_are_rejected(override):
    with pytest.raises(ValueError):
        replace(fs.FableParameters(), **override).validate()


@pytest.mark.parametrize("override", [{"omega_b0": -0.01}, {"bilinear_s0": 0.0}, {"xi": float("nan")}])
def test_negative_density_zero_bilinear_and_nan_are_rejected(override):
    with pytest.raises(ValueError):
        replace(fs.FableParameters(), **override).validate()


def test_parameter_helpers_only_change_xi():
    p = fs.FableParameters()
    assert p.torsion_free() == replace(p, xi=0.0)
    assert p.with_xi(0.25) == replace(p, xi=0.25)
    assert p.torsion_free().validate() is None


def test_closed_form_solves_the_dilution_law():
    """d(ln S)/dN must equal -nu(s(N)) along the logistic bridge flow."""
    p = fs.FableParameters()
    grid = np.linspace(-6.0, 1.0, 400)
    step = 1e-6
    numerical = (
        fs.log_bilinear_closed_form(grid + step, p)
        - fs.log_bilinear_closed_form(grid - step, p)
    ) / (2 * step)
    analytic = -fs.dilution_exponent(fs.bridge_closed_form(grid, p), p.xi)
    assert np.max(np.abs(numerical - analytic)) < 1e-8


def test_scipy_background_tracks_the_closed_form():
    solution = fs.solve_background()
    assert np.max(np.abs(solution.closed_form_residual)) < 1e-9
    assert np.max(np.abs(solution.bridge_residual)) < 1e-9


def test_sector_is_covariantly_conserved():
    solution = fs.solve_background()
    assert np.max(np.abs(solution.continuity_residual)) < 1e-13


def test_massless_torsion_free_limit_reproduces_the_benchmark():
    """The published claim: w = n - 1 exactly, so n = 1 + (-0.764) gives
    a flat, non-evolving w = -0.764."""
    solution = fs.solve_background(fs.FableParameters().torsion_free())
    assert np.allclose(solution.w_potential, fs.BENCHMARK_W, atol=1e-12)
    assert np.allclose(solution.w_dust, 0.0, atol=1e-12)


def test_unified_dark_sector_interpolates_dust_and_dark_energy():
    p = fs.FableParameters().torsion_free()
    solution = fs.solve_background(p, n_start=-12.0, n_end=6.0, samples=2001)
    early = solution.w_fable[0]
    late = solution.w_fable[-1]
    assert early == pytest.approx(0.0, abs=1e-3)
    assert late == pytest.approx(fs.BENCHMARK_W, abs=1e-3)
    assert np.all(np.diff(solution.w_fable) <= 1e-12)


def test_second_mechanism_needs_torsion_and_survives_without_a_potential():
    """w moves with the bridge field only when xi is nonzero, and it moves even
    when the potential index is held fixed."""
    p = fs.FableParameters()
    frozen = fs.solve_background(p.torsion_free())
    active = fs.solve_background(p.with_xi(0.25))
    assert np.ptp(frozen.w_potential) == pytest.approx(0.0, abs=1e-14)
    assert np.ptp(active.w_potential) > 1e-3
    assert np.ptp(frozen.w_dust) == pytest.approx(0.0, abs=1e-14)
    assert np.ptp(active.w_dust) > 1e-3


def test_no_phantom_crossing():
    """The sector never crosses w = -1, so the model cannot fake a phantom."""
    solution = fs.solve_background(n_start=-12.0, n_end=6.0, samples=2001)
    assert np.all(solution.w_fable > -1.0)
    assert np.all(solution.w_potential > -1.0)


def test_reference_binary_is_required_not_optional():
    """A missing Rust build must raise, never silently fall back to scipy."""
    with pytest.raises(FileNotFoundError):
        fs.run_reference(
            Path("artifacts/does-not-exist"),
            repository_root=REPOSITORY_ROOT,
            binary=REPOSITORY_ROOT / "fable_cosmo_rs/target/release/not-a-binary",
        )


def test_solve_background_rejects_a_degenerate_grid():
    with pytest.raises(ValueError):
        fs.solve_background(n_start=1.0, n_end=0.0)
    with pytest.raises(ValueError):
        fs.solve_background(samples=1)


def test_background_columns_match_the_reference_csv_schema():
    """Every column the Rust binary writes has a same-named field here, so
    the two integrators can be compared column for column."""
    fields = set(fs.FableBackground.__dataclass_fields__)
    assert set(fs.PHYSICAL_COLUMNS) <= fields
    residuals = {"continuity_residual", "closed_form_residual", "bridge_residual"}
    assert fields == set(fs.PHYSICAL_COLUMNS) | residuals


def test_compare_to_reference_rejects_a_missing_column_and_a_grid_mismatch():
    background = fs.solve_background(samples=11)
    with pytest.raises(KeyError):
        fs.compare_to_reference(background, {}, columns=("e_folds",))
    with pytest.raises(ValueError):
        fs.compare_to_reference(
            background, {"e_folds": np.zeros(10)}, columns=("e_folds",)
        )


@pytest.mark.skipif(
    not (REPOSITORY_ROOT / fs.REFERENCE_BINARY).is_file(),
    reason="fable_cosmo_rs is not built; run `cargo build --release` in fable_cosmo_rs/",
)
def test_scipy_agrees_with_the_pure_rust_reference(tmp_path):
    """The three-way check, as a test: SciPy Radau against pure-Rust CVODE on
    every physical column, to the mixed tolerance the published gates use."""
    run = fs.run_reference(tmp_path, repository_root=REPOSITORY_ROOT)
    assert "SUCCESS: every gated invariant holds." in run.stdout
    differences = fs.compare_to_reference(fs.solve_background(), run.table())
    assert set(differences) == set(fs.PHYSICAL_COLUMNS)
    worst = max(differences.values())
    assert worst < 1e-10, differences
