"""Tests for surfactantkit.thermodynamics."""

import math
import pytest

from surfactantkit.thermodynamics import (
    cmc_to_mole_fraction,
    counterion_binding_degree,
    gibbs_free_energy_micellization,
    mass_action_free_energy_of_micellization,
    vant_hoff_enthalpy,
    vant_hoff_multi_point_fit,
    entropy_micellization,
)

R_GAS = 8.314462618


def test_cmc_to_mole_fraction_dilute_approximation():
    # for a dilute CMC, X_cmc should be very close to CMC/55.5
    x = cmc_to_mole_fraction(0.008)  # 8 mM
    assert x == pytest.approx(0.008 / 55.5, rel=1e-3)


def test_cmc_to_mole_fraction_rejects_nonpositive():
    with pytest.raises(ValueError):
        cmc_to_mole_fraction(0.0)


def test_counterion_binding_degree_basic():
    # slope drops by half above the CMC -> beta = 0.5
    assert counterion_binding_degree(slope_below_cmc=2.0, slope_above_cmc=1.0) == pytest.approx(0.5)


def test_counterion_binding_degree_rejects_wrong_slope_order():
    with pytest.raises(ValueError):
        counterion_binding_degree(slope_below_cmc=1.0, slope_above_cmc=2.0)


def test_counterion_binding_degree_rejects_nonpositive():
    with pytest.raises(ValueError):
        counterion_binding_degree(0.0, 1.0)


def test_gibbs_free_energy_is_negative_for_spontaneous_micellization():
    """Sanity check: for a realistic ionic-surfactant CMC and counterion
    factor, deltaG_mic should be negative (spontaneous) and land in the
    commonly-reported -20 to -50 kJ/mol range for small-molecule
    surfactants."""
    x_cmc = cmc_to_mole_fraction(0.008)  # ~8 mM, SDS-like
    dg = gibbs_free_energy_micellization(x_cmc, 298.15, counterion_factor=1.6)
    assert dg < 0
    assert -60.0 < dg < -10.0


def test_gibbs_free_energy_nonionic_vs_ionic_factor():
    # same CMC, but the ionic (counterion_factor > 1) case should be
    # more negative than the nonionic (factor=1) case, since 2-beta > 1
    x_cmc = cmc_to_mole_fraction(0.008)
    dg_nonionic = gibbs_free_energy_micellization(x_cmc, 298.15, counterion_factor=1.0)
    dg_ionic = gibbs_free_energy_micellization(x_cmc, 298.15, counterion_factor=1.6)
    assert dg_ionic < dg_nonionic


def test_gibbs_free_energy_rejects_bad_mole_fraction():
    with pytest.raises(ValueError):
        gibbs_free_energy_micellization(0.0, 298.15)
    with pytest.raises(ValueError):
        gibbs_free_energy_micellization(1.0, 298.15)


def test_vant_hoff_enthalpy_round_trip():
    """Construct two (CMC, T) points consistent with a chosen deltaH,
    then verify the function recovers it -- mathematically guaranteed.

    Real bug found 2026-09-04: this test used to construct its synthetic
    data with the SAME (1/T1 - 1/T2) convention the function itself used,
    so it stayed green through a real sign bug in vant_hoff_enthalpy() --
    a round-trip test only proves internal self-consistency, not that the
    convention matches the true physics. Caught via an unaugmented Claude
    pilot answer disagreeing with the tool, then confirmed independently
    from the Gibbs-Helmholtz equation and against this project's own
    existing Fu et al. 2019 literature-validation entry (whose two-point
    estimates had the opposite sign from the paper's own local-derivative
    value, previously misread as pure two-point-vs-polynomial magnitude
    noise). Fixed here to use the corrected (1/T2 - 1/T1) convention,
    matching the now-fixed function; see test below for an independent,
    non-round-trip sign check against real physical reasoning.
    """
    R = 8.314462618
    T1, T2 = 293.15, 313.15
    x1 = 1.5e-4
    delta_h_true_kJ = -15.0  # typical small negative micellization enthalpy

    delta_inv_T = (1.0 / T2) - (1.0 / T1)
    ln_x2 = math.log(x1) + (delta_h_true_kJ * 1000.0 * delta_inv_T) / R
    x2 = math.exp(ln_x2)

    delta_h_recovered = vant_hoff_enthalpy(x1, T1, x2, T2)
    assert delta_h_recovered == pytest.approx(delta_h_true_kJ, abs=1e-6)


def test_vant_hoff_enthalpy_sign_matches_physical_reasoning():
    """Independent, non-round-trip sign check: if CMC DECREASES as
    temperature RISES, micellization is more favourable at higher T,
    which by the van't Hoff/Le Chatelier argument must be an ENDOTHERMIC
    (positive deltaH) process -- this is real physical reasoning, not
    derived from the function under test, so it cannot share a bug with
    it the way the round-trip test above could (and did)."""
    T1, T2 = 293.15, 313.15
    cmc1, cmc2 = 4.0e-5, 3.8e-5  # CMC2 < CMC1: CMC falls as T rises
    delta_h = vant_hoff_enthalpy(cmc1, T1, cmc2, T2)
    assert delta_h > 0, "CMC falling with rising T must give a positive (endothermic) deltaH"


def test_vant_hoff_enthalpy_rejects_equal_temperatures():
    with pytest.raises(ValueError):
        vant_hoff_enthalpy(1e-4, 298.15, 2e-4, 298.15)


# --- vant_hoff_multi_point_fit ----------------------------------------------
# Source: Kantonen, Henriksen & Gilson, BBA Gen. Subj. 1862 (2018) 692-704,
# Eq. 4 (verified via WebFetch before implementing, not guessed) --
# deltaG(Tj) = deltaH(Tr) - Tj*deltaS(Tr) + deltaCp(Tr)*[(Tj-Tr) - Tj*ln(Tj/Tr)].
# Validated via mathematically-guaranteed round trips: construct X_cmc(Tj)
# for chosen true deltaH/deltaS/deltaCp via this exact equation, verify the
# fit recovers them.


def _synthetic_series(true_dH_kJ, true_dS_J, true_dCp_J, tr, temps, counterion_factor=1.0):
    true_dS_kJ = true_dS_J / 1000.0
    true_dCp_kJ = true_dCp_J / 1000.0
    x_values = []
    for t in temps:
        g = (t - tr) - t * math.log(t / tr)
        dg = true_dH_kJ - t * true_dS_kJ + true_dCp_kJ * g
        x_values.append(math.exp(dg * 1000.0 / (counterion_factor * R_GAS * t)))
    return x_values


def test_vant_hoff_multi_point_fit_recovers_true_parameters_round_trip():
    true_dH_kJ, true_dS_J, true_dCp_J = -20.0, -50.0, -300.0
    tr = 298.15
    temps = [278.15, 288.15, 298.15, 308.15, 318.15]
    x_values = _synthetic_series(true_dH_kJ, true_dS_J, true_dCp_J, tr, temps)

    result = vant_hoff_multi_point_fit(x_values, temps, reference_temperature_K=tr)
    assert result.delta_H_kJ_per_mol == pytest.approx(true_dH_kJ, abs=1e-6)
    assert result.delta_S_J_per_mol_K == pytest.approx(true_dS_J, abs=1e-4)
    assert result.delta_Cp_J_per_mol_K == pytest.approx(true_dCp_J, abs=1e-3)
    assert result.r_squared == pytest.approx(1.0, abs=1e-9)
    assert result.n_points == 5
    assert result.reference_temperature_K == pytest.approx(tr)


def test_vant_hoff_multi_point_fit_default_reference_is_mean_temperature():
    true_dH_kJ, true_dS_J, true_dCp_J = -18.0, -40.0, -250.0
    tr = 300.0
    temps = [283.15, 293.15, 303.15, 313.15]
    x_values = _synthetic_series(true_dH_kJ, true_dS_J, true_dCp_J, tr, temps)

    result = vant_hoff_multi_point_fit(x_values, temps)  # no explicit reference
    assert result.reference_temperature_K == pytest.approx(sum(temps) / len(temps))
    # deltaCp and the fitted deltaG curve are reference-independent by construction:
    assert result.delta_Cp_J_per_mol_K == pytest.approx(true_dCp_J, abs=1e-3)
    assert result.r_squared == pytest.approx(1.0, abs=1e-9)


def test_vant_hoff_multi_point_fit_with_ionic_counterion_factor():
    true_dH_kJ, true_dS_J, true_dCp_J = -22.0, -55.0, -320.0
    tr = 298.15
    temps = [288.15, 298.15, 308.15, 318.15]
    factor = 1.6
    x_values = _synthetic_series(true_dH_kJ, true_dS_J, true_dCp_J, tr, temps, counterion_factor=factor)

    result = vant_hoff_multi_point_fit(x_values, temps, counterion_factor=factor, reference_temperature_K=tr)
    assert result.delta_H_kJ_per_mol == pytest.approx(true_dH_kJ, abs=1e-6)
    assert result.delta_Cp_J_per_mol_K == pytest.approx(true_dCp_J, abs=1e-3)


def test_vant_hoff_multi_point_fit_rejects_bad_inputs():
    with pytest.raises(ValueError):
        vant_hoff_multi_point_fit([1e-4, 2e-4], [293.15, 303.15])  # too few points
    with pytest.raises(ValueError):
        vant_hoff_multi_point_fit([1e-4, 2e-4, 3e-4], [293.15, 303.15])  # mismatched lengths
    with pytest.raises(ValueError):
        vant_hoff_multi_point_fit([0.0, 2e-4, 3e-4], [293.15, 303.15, 313.15])  # bad mole fraction
    with pytest.raises(ValueError):
        vant_hoff_multi_point_fit([1e-4, 2e-4, 3e-4], [293.15, 0.0, 313.15])  # bad temperature
    with pytest.raises(ValueError):
        vant_hoff_multi_point_fit([1e-4, 2e-4, 3e-4], [298.15, 298.15, 298.15])  # degenerate (all same T)


# --- mass_action_free_energy_of_micellization (alternative-methods sweep,
# 2026-09-10) -- eq. 21 of Boneva-Aroca et al., RSC Advances 13 (2023) 9387
# (open access, quoted verbatim via WebFetch, not guessed), the same
# well-known classical single-equilibrium mass-action result this project's
# own literature review flagged via Rusanov 2014 (paywalled primary source).


def test_mass_action_matches_cited_formula_directly():
    """Independent recomputation in the test itself, not copy-pasted
    from the implementation."""
    x_cmc, n, T = 1.5e-4, 60, 298.15
    R = 8.314462618
    expected_kJ = (R * T * ((1.0 + 1.0 / n) * math.log(x_cmc) - (1.0 / n) * math.log(n))) / 1000.0
    assert mass_action_free_energy_of_micellization(x_cmc, n, T) == pytest.approx(expected_kJ)


def test_mass_action_converges_to_pseudo_phase_separation_at_large_n():
    """The real, verified property tying the two models together: as n
    grows, mass_action_free_energy_of_micellization must converge to
    gibbs_free_energy_micellization(counterion_factor=1) -- the paper's
    own stated large-n limit, not just an incidental numerical match."""
    x_cmc, T = 2.0e-4, 298.15
    pseudo_phase = gibbs_free_energy_micellization(x_cmc, T, counterion_factor=1.0)
    for n in (10, 100, 1000, 100000):
        mass_action = mass_action_free_energy_of_micellization(x_cmc, n, T)
        diff = abs(mass_action - pseudo_phase)
        if n == 10:
            prev_diff = diff
        else:
            assert diff < prev_diff
            prev_diff = diff
    assert mass_action_free_energy_of_micellization(x_cmc, 100000, T) == pytest.approx(pseudo_phase, abs=1e-3)


def test_mass_action_diverges_meaningfully_from_pseudo_phase_at_small_n():
    """For a small aggregation number, the correction should be real and
    non-negligible (not silently ~0), matching the source paper's own
    finding that the correction remains ~0.1 (in RT units) even at
    n~100 -- here at a much smaller, more divergent n."""
    x_cmc, n, T = 2.0e-4, 15, 298.15
    pseudo_phase = gibbs_free_energy_micellization(x_cmc, T, counterion_factor=1.0)
    mass_action = mass_action_free_energy_of_micellization(x_cmc, n, T)
    assert abs(mass_action - pseudo_phase) > 0.5  # kJ/mol, a real, non-trivial difference


def test_mass_action_rejects_bad_inputs():
    with pytest.raises(ValueError):
        mass_action_free_energy_of_micellization(0.0, 50, 298.15)
    with pytest.raises(ValueError):
        mass_action_free_energy_of_micellization(1.5e-4, 1, 298.15)  # n must be > 1
    with pytest.raises(ValueError):
        mass_action_free_energy_of_micellization(1.5e-4, 50, 0.0)


def test_entropy_micellization_completes_triad_by_construction():
    # deltaS = (deltaH - deltaG) / T -- purely definitional, must hold exactly
    dh, dg, T = -15.0, -35.0, 298.15
    ds = entropy_micellization(dg, dh, T)
    assert ds == pytest.approx(((dh - dg) / T) * 1000.0)


def test_entropy_micellization_rejects_bad_temperature():
    with pytest.raises(ValueError):
        entropy_micellization(-35.0, -15.0, 0.0)
