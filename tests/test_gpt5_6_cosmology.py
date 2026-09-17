import numpy as np

from src.gpt5_6_cosmology import (
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
    scalar_bilinear,
    solution_diagnostics,
    solve_background,
    spinor_potential,
    spinor_potential_derivative,
    spinor_thermodynamics,
    vector_current,
)


PARAMETERS = CosmologyParameters()


def test_observational_benchmarks_are_distinct_and_normalized() -> None:
    constant_w = UNITE_ONLY_CONSTANT_W_BENCHMARK.cosmology_parameters()
    evolving_w = UNITE_BAO_CMB_CPL_BENCHMARK.cosmology_parameters()

    assert UNITE_ONLY_CONSTANT_W_BENCHMARK.datasets == ("Unite",)
    assert UNITE_ONLY_CONSTANT_W_BENCHMARK.wa is None
    assert UNITE_BAO_CMB_CPL_BENCHMARK.datasets == ("Unite", "BAO", "CMB")
    assert constant_w.omega_m0 == 0.197
    assert constant_w.w0 == -0.764
    assert constant_w.wa == 0.0
    assert evolving_w.omega_m0 == 0.305
    assert evolving_w.w0 == -0.861
    assert evolving_w.wa == -0.60

    scale_factor = np.geomspace(0.1, 2.0, 31)
    np.testing.assert_allclose(cpl_w(scale_factor, constant_w), -0.764)
    assert crossing_scale_factors(scale_factor, cpl_w(scale_factor, constant_w), -1.0) == []
    density = spinor_potential(
        bilinear_from_scale_factor(scale_factor, constant_w), constant_w
    )
    expected_density = constant_w.omega_gpt5_6_0 * scale_factor ** (
        -3.0 * (1.0 + constant_w.w0)
    )
    np.testing.assert_allclose(density, expected_density)
    density_today = spinor_potential(constant_w.spinor_bilinear0, constant_w)
    np.testing.assert_allclose(
        constant_w.omega_m0 + constant_w.omega_r0 + density_today,
        1.0,
    )


def test_sixteen_component_clifford_representation() -> None:
    metric = np.diag([1.0, -1.0, -1.0, -1.0])
    identity = np.eye(16)
    assert len(GAMMA_16) == 4
    assert all(gamma.shape == (16, 16) for gamma in GAMMA_16)
    for mu in range(4):
        for nu in range(4):
            anticommutator = GAMMA_16[mu] @ GAMMA_16[nu] + GAMMA_16[nu] @ GAMMA_16[mu]
            np.testing.assert_allclose(anticommutator, 2.0 * metric[mu, nu] * identity)


def test_gpt5_6_condensate_has_expected_bilinear_and_no_spatial_current() -> None:
    scale_factor = 0.73
    field = gpt5_6(scale_factor)
    assert field.shape == (16,)
    np.testing.assert_allclose(scalar_bilinear(field), scale_factor**-3)
    np.testing.assert_allclose(vector_current(field)[1:], 0.0, atol=1.0e-14)


def test_reconstructed_potential_derivative() -> None:
    bilinear = np.geomspace(0.1, 100.0, 101)
    step = 1.0e-5 * bilinear
    finite_difference = (
        spinor_potential(bilinear + step, PARAMETERS)
        - spinor_potential(bilinear - step, PARAMETERS)
    ) / (2.0 * step)
    analytic = spinor_potential_derivative(bilinear, PARAMETERS)
    np.testing.assert_allclose(analytic, finite_difference, rtol=2.0e-8, atol=1.0e-10)


def test_pressure_and_equation_of_state_exactly_reproduce_cpl() -> None:
    scale_factor = np.geomspace(0.05, 2.5, 301)
    terms = spinor_thermodynamics(scale_factor, PARAMETERS)
    np.testing.assert_allclose(terms["pressure"], terms["kinetic"] - terms["potential"])
    np.testing.assert_allclose(terms["density"], terms["potential"])
    np.testing.assert_allclose(terms["w"], cpl_w(scale_factor, PARAMETERS), rtol=2.0e-14)


def test_homogeneous_energy_momentum_tensor_is_perfect_fluid() -> None:
    tensor = energy_momentum_mixed(1.0, PARAMETERS)
    terms = spinor_thermodynamics(1.0, PARAMETERS)
    expected = np.diag(
        [terms["density"], -terms["pressure"], -terms["pressure"], -terms["pressure"]]
    )
    np.testing.assert_allclose(tensor, expected)


def test_cosmological_constant_and_dust_limits() -> None:
    cosmological_constant = CosmologyParameters(w0=-1.0, wa=0.0)
    dust = CosmologyParameters(w0=0.0, wa=0.0)
    scale_factor = np.geomspace(0.1, 2.0, 31)
    np.testing.assert_allclose(spinor_thermodynamics(scale_factor, cosmological_constant)["w"], -1.0)
    np.testing.assert_allclose(spinor_thermodynamics(scale_factor, dust)["w"], 0.0, atol=1.0e-15)
    dust_bilinear = bilinear_from_scale_factor(scale_factor, dust)
    np.testing.assert_allclose(
        spinor_potential(dust_bilinear, dust) / dust_bilinear,
        np.full_like(scale_factor, dust.omega_gpt5_6_0 / dust.spinor_bilinear0),
    )


def test_numerical_background_matches_bilinear_and_potential_invariants() -> None:
    solution = solve_background(PARAMETERS)
    diagnostics = solution_diagnostics(solution, PARAMETERS)
    assert diagnostics["max_relative_bilinear_error"] < 2.0e-9
    assert diagnostics["max_relative_density_error"] < 2.0e-9
    assert diagnostics["max_relative_friedmann_error"] < 2.0e-9
    np.testing.assert_allclose(diagnostics["w_today"], PARAMETERS.w0)
    expected_crossing = 1.0 + (1.0 + PARAMETERS.w0) / PARAMETERS.wa
    np.testing.assert_allclose(diagnostics["phantom_crossing_a"], expected_crossing, rtol=2.0e-6)
    np.testing.assert_allclose(diagnostics["nec_crossing_a"], expected_crossing, rtol=2.0e-6)
    acceleration_crossings = diagnostics["acceleration_crossing_a"]
    assert isinstance(acceleration_crossings, list)
    assert len(acceleration_crossings) == 2
    assert acceleration_crossings[0] < 1.0 < acceleration_crossings[1]


def test_cpl_values_and_null_energy_condition_change_sign_at_phantom_crossing() -> None:
    scale_factor = np.array([0.1, 1.0, 2.0])
    np.testing.assert_allclose(cpl_w(scale_factor, PARAMETERS), [-1.401, -0.861, -0.261])

    crossing = 1.0 + (1.0 + PARAMETERS.w0) / PARAMETERS.wa
    crossing_terms = spinor_thermodynamics(crossing, PARAMETERS)
    np.testing.assert_allclose(crossing_terms["w"], -1.0, atol=1.0e-14)
    np.testing.assert_allclose(
        crossing_terms["density"] + crossing_terms["pressure"],
        0.0,
        atol=1.0e-14,
    )

    before = spinor_thermodynamics(0.5, PARAMETERS)
    after = spinor_thermodynamics(1.0, PARAMETERS)
    assert before["density"] + before["pressure"] < 0.0
    assert after["density"] + after["pressure"] > 0.0


def test_friedmann_closure_and_acceleration_transitions() -> None:
    solution = solve_background(PARAMETERS)
    np.testing.assert_allclose(
        solution.hubble_over_h0**2,
        hubble_squared(solution.scale_factor, PARAMETERS),
        rtol=2.0e-9,
    )

    crossings = solution_diagnostics(solution, PARAMETERS)["acceleration_crossing_a"]
    assert isinstance(crossings, list)
    np.testing.assert_allclose(crossings, [0.57194538, 1.80120338], rtol=2.0e-6)
    assert solution.deceleration[solution.scale_factor < crossings[0]][-1] > 0.0
    assert solution.deceleration[np.argmin(np.abs(solution.scale_factor - 1.0))] < 0.0
    assert solution.deceleration[solution.scale_factor > crossings[1]][0] > 0.0