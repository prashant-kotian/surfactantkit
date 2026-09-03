"""Tests for surfactantkit.electrostatics and surfactantkit.dynamics."""

import math
import pytest

from surfactantkit.electrostatics import (
    ionic_strength,
    debye_length,
    henry_function,
    zeta_potential_henry,
)
from surfactantkit.dynamics import hydrodynamic_radius_stokes_einstein


def test_ionic_strength_symmetric_1_1_salt():
    # 0.1 M NaCl: I = 0.5*(0.1*1^2 + 0.1*1^2) = 0.1 M
    assert ionic_strength({1: 0.1, -1: 0.1}) == pytest.approx(0.1)


def test_ionic_strength_2_1_salt():
    # 0.1 M CaCl2: Ca2+ at 0.1 M (z=2), Cl- at 0.2 M (z=-1)
    # I = 0.5*(0.1*4 + 0.2*1) = 0.5*(0.4+0.2) = 0.3 M
    assert ionic_strength({2: 0.1, -1: 0.2}) == pytest.approx(0.3)


def test_debye_length_matches_classic_textbook_value():
    """0.1 M 1:1 electrolyte (e.g. NaCl) at 25C has a Debye length of
    ~0.96 nm -- one of the most commonly cited reference values in
    colloid science (e.g. Israelachvili). Confirmed here by direct
    calculation, not just asserted from memory."""
    ld = debye_length(ionic_strength_M=0.1, temperature_K=298.15)
    assert ld == pytest.approx(0.961, abs=0.01)


def test_debye_length_decreases_with_ionic_strength():
    # higher salt -> shorter screening length (physically required)
    ld_low = debye_length(0.001, 298.15)
    ld_high = debye_length(0.1, 298.15)
    assert ld_high < ld_low


def test_debye_length_rejects_bad_input():
    with pytest.raises(ValueError):
        debye_length(0.0, 298.15)
    with pytest.raises(ValueError):
        debye_length(0.1, -1.0)


def test_henry_function_limits():
    assert henry_function("huckel") == pytest.approx(1.0)
    assert henry_function("smoluchowski") == pytest.approx(1.5)
    assert henry_function("HUCKEL") == pytest.approx(1.0)  # case-insensitive


def test_henry_function_no_default_for_unknown_regime():
    with pytest.raises(ValueError):
        henry_function("something_else")


def test_zeta_potential_henry_round_trip():
    """Construct a mobility that produces a chosen zeta (Smoluchowski
    regime), then verify the function recovers that same zeta --
    mathematically guaranteed, no external reference needed."""
    eps0 = 8.8541878128e-12
    eps_r = 78.4
    eta_mPas = 0.89
    zeta_true_mV = 45.0

    f_ka = henry_function("smoluchowski")
    zeta_V = zeta_true_mV / 1000.0
    eta_SI = eta_mPas * 1e-3
    mobility_SI = (2.0 * eps_r * eps0 * zeta_V * f_ka) / (3.0 * eta_SI)
    mobility_um_cm_per_Vs = mobility_SI / (1e-6 * 1e-2)

    zeta_recovered = zeta_potential_henry(mobility_um_cm_per_Vs, eta_mPas, "smoluchowski")
    assert zeta_recovered == pytest.approx(zeta_true_mV, abs=1e-6)


def test_zeta_potential_rejects_bad_viscosity():
    with pytest.raises(ValueError):
        zeta_potential_henry(1.0, 0.0, "huckel")


def test_hydrodynamic_radius_round_trip():
    """Construct a diffusion coefficient that produces a chosen R_h,
    then verify the function recovers it -- mathematically guaranteed."""
    k_B = 1.380649e-23
    T = 298.15
    eta_mPas = 0.89
    r_h_true_nm = 3.5

    eta_SI = eta_mPas * 1e-3
    r_h_m = r_h_true_nm * 1e-9
    D_SI = (k_B * T) / (6.0 * math.pi * eta_SI * r_h_m)
    D_cm2_per_s = D_SI * 1e4

    r_h_recovered = hydrodynamic_radius_stokes_einstein(D_cm2_per_s, eta_mPas, T)
    assert r_h_recovered == pytest.approx(r_h_true_nm, abs=1e-6)


def test_hydrodynamic_radius_rejects_bad_input():
    with pytest.raises(ValueError):
        hydrodynamic_radius_stokes_einstein(0.0, 0.89)


def test_hydrodynamic_radius_matches_independent_hand_calculation():
    """Real-number check, not a round-trip: found during the 2026-09-04
    library-wide audit (prompted by the vant_hoff_enthalpy sign bug) that
    this function had only a round-trip test plus a documented NEGATIVE
    control (literature_validation_notes.md's Milone et al. entry, where
    the mismatch was expected and explained away, not a positive match).
    This is an independent hand-derivation using the exact D value from
    that same negative-control paper (D=1.56e-6 cm^2/s, water at 298.15K,
    eta=0.89 mPa.s) worked out by hand outside this function's own code
    path, giving 1.573 nm -- now locked in as a real regression check."""
    r_h = hydrodynamic_radius_stokes_einstein(
        diffusion_coefficient_cm2_per_s=1.56e-6, viscosity_mPas=0.89, temperature_K=298.15)
    assert r_h == pytest.approx(1.573, abs=0.001)
    with pytest.raises(ValueError):
        hydrodynamic_radius_stokes_einstein(1e-6, -1.0)
