"""Tests for the Szyszkowski surface-tension equation (adsorption.py),
wetting.py, and solubilization.py."""

import math
import pytest

from surfactantkit.adsorption import szyszkowski_surface_tension, gibbs_gamma_max
from surfactantkit.wetting import work_of_adhesion, spreading_coefficient, capillary_number
from surfactantkit.solubilization import molar_solubilization_ratio


def test_szyszkowski_returns_gamma0_at_zero_concentration():
    # boundary condition, follows directly from ln(1+0)=0
    g = szyszkowski_surface_tension(0.0, gamma0_mN_m=72.0, gamma_max_mol_per_m2=3.0e-6, K=50.0,
                                     system_type="nonionic")
    assert g == pytest.approx(72.0)


def test_szyszkowski_decreases_monotonically_with_concentration():
    g_low = szyszkowski_surface_tension(1.0, 72.0, 3.0e-6, 50.0, "nonionic")
    g_high = szyszkowski_surface_tension(10.0, 72.0, 3.0e-6, 50.0, "nonionic")
    assert g_high < g_low < 72.0


def test_szyszkowski_high_concentration_limit_matches_gibbs_slope_method():
    """At high concentration (K*C >> 1), the numerical slope
    d(gamma)/d(ln C) of the Szyszkowski curve should recover the same
    Gamma_max that was used to construct it, when fed through
    gibbs_gamma_max -- these are two independently-implemented physics
    relations that must agree exactly in this limit; a strong internal
    cross-check, not just a unit test."""
    gamma0, gamma_max_true, K, T = 72.0, 3.0e-6, 50.0, 298.15
    C1, C2 = 100.0, 105.0  # K*C ~ 5000-5250, well into the saturation limit
    g1 = szyszkowski_surface_tension(C1, gamma0, gamma_max_true, K, "nonionic", T)
    g2 = szyszkowski_surface_tension(C2, gamma0, gamma_max_true, K, "nonionic", T)
    slope = (g2 - g1) / (math.log(C2) - math.log(C1))
    recovered_gamma_max = gibbs_gamma_max(slope, system_type="nonionic", temperature_k=T)
    assert recovered_gamma_max == pytest.approx(gamma_max_true, rel=0.005)


def test_szyszkowski_rejects_bad_input():
    with pytest.raises(ValueError):
        szyszkowski_surface_tension(-1.0, 72.0, 3.0e-6, 50.0, "nonionic")
    with pytest.raises(ValueError):
        szyszkowski_surface_tension(1.0, 72.0, 0.0, 50.0, "nonionic")
    with pytest.raises(ValueError):
        szyszkowski_surface_tension(1.0, 72.0, 3.0e-6, -1.0, "nonionic")
    with pytest.raises(ValueError):
        szyszkowski_surface_tension(1.0, 72.0, 3.0e-6, 50.0, "not_a_real_system_type")


def test_gibbs_gamma_max_rejects_invalid_system_type():
    """Regression test for the 2026-09-04 hardening: system_type replaced
    a raw n_factor float precisely because a caller (including an LLM
    using this as a tool) must not be able to silently pass an
    unjustified numeric guess. No default is provided and any string
    outside the three real system types must raise."""
    with pytest.raises(ValueError):
        gibbs_gamma_max(-5.0, "made_up_system_type", 298.15)
    with pytest.raises(TypeError):
        gibbs_gamma_max(-5.0, temperature_k=298.15)  # system_type omitted entirely -- no default


def test_work_of_adhesion_complete_wetting_limit():
    # theta = 0 (complete wetting) -> W_a = 2*gamma_LV exactly
    assert work_of_adhesion(72.0, 0.0) == pytest.approx(144.0)


def test_work_of_adhesion_no_wetting_limit():
    # theta = 180 -> W_a = 0 exactly
    assert work_of_adhesion(72.0, 180.0) == pytest.approx(0.0, abs=1e-9)


def test_work_of_adhesion_rejects_bad_input():
    with pytest.raises(ValueError):
        work_of_adhesion(-1.0, 90.0)
    with pytest.raises(ValueError):
        work_of_adhesion(72.0, 200.0)


def test_spreading_coefficient_complete_wetting_boundary():
    # theta = 0 -> S = 0 exactly (the only case where S >= 0 is achieved)
    assert spreading_coefficient(72.0, 0.0) == pytest.approx(0.0, abs=1e-9)


def test_spreading_coefficient_is_never_positive():
    for theta in [0.0, 30.0, 90.0, 150.0, 180.0]:
        assert spreading_coefficient(72.0, theta) <= 1e-9


def test_capillary_number_matches_direct_formula():
    ca = capillary_number(viscosity_mPas=1.0, velocity_m_per_s=1e-5, interfacial_tension_mN_m=0.01)
    expected = (1.0e-3 * 1e-5) / (0.01e-3)
    assert ca == pytest.approx(expected)


def test_capillary_number_increases_as_interfacial_tension_drops():
    # this IS the whole point of surfactant EOR -- lowering IFT raises Ca
    ca_normal_ift = capillary_number(1.0, 1e-5, 20.0)
    ca_ultralow_ift = capillary_number(1.0, 1e-5, 0.001)
    assert ca_ultralow_ift > ca_normal_ift


def test_capillary_number_rejects_nonpositive():
    with pytest.raises(ValueError):
        capillary_number(0.0, 1e-5, 20.0)


def test_molar_solubilization_ratio_basic():
    # 0.5 mM extra solubilized (beyond intrinsic 0.01 mM) per (10-2) mM micellized surfactant
    msr = molar_solubilization_ratio(
        total_solubilized_M=0.51e-3, intrinsic_water_solubility_M=0.01e-3,
        surfactant_concentration_M=10e-3, cmc_M=2e-3,
    )
    assert msr == pytest.approx((0.51e-3 - 0.01e-3) / (10e-3 - 2e-3))


def test_molar_solubilization_ratio_rejects_below_cmc():
    with pytest.raises(ValueError):
        molar_solubilization_ratio(0.1e-3, 0.01e-3, 1e-3, 2e-3)


def test_molar_solubilization_ratio_rejects_inconsistent_solubility():
    with pytest.raises(ValueError):
        molar_solubilization_ratio(0.005e-3, 0.01e-3, 10e-3, 2e-3)
