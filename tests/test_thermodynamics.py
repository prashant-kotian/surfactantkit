"""Tests for surfactantkit.thermodynamics."""

import math
import pytest

from surfactantkit.thermodynamics import (
    cmc_to_mole_fraction,
    counterion_binding_degree,
    gibbs_free_energy_micellization,
    vant_hoff_enthalpy,
    entropy_micellization,
)


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


def test_entropy_micellization_completes_triad_by_construction():
    # deltaS = (deltaH - deltaG) / T -- purely definitional, must hold exactly
    dh, dg, T = -15.0, -35.0, 298.15
    ds = entropy_micellization(dg, dh, T)
    assert ds == pytest.approx(((dh - dg) / T) * 1000.0)


def test_entropy_micellization_rejects_bad_temperature():
    with pytest.raises(ValueError):
        entropy_micellization(-35.0, -15.0, 0.0)
