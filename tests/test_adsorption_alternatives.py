"""Tests for surfactantkit.adsorption's Frumkin isotherm, the first
alternative-method addition from the 2026-09-07/08 method-alternatives
literature review (see benchmark/METHOD_ALTERNATIVES_LITERATURE_REVIEW.md).

Frumkin generalizes the existing Szyszkowski/Langmuir isotherm with a
dimensionless lateral-interaction parameter `a` (a > 0 = attraction,
a < 0 = repulsion, a == 0 = plain Langmuir). Sign convention and closed
form verified against two independent real sources before implementing
(see frumkin_theta's docstring) -- no numeric literature validation case
was found for this specific isotherm during that pass, so validation
here leans on the one thing that IS guaranteed: a=0 must reduce EXACTLY
to the already-literature-validated Szyszkowski code, plus qualitative
physical-direction checks (attraction steepens the isotherm, repulsion
flattens it) that don't depend on any external numeric source.
"""

import math

import pytest

from surfactantkit.adsorption import (
    frumkin_surface_tension,
    frumkin_theta,
    szyszkowski_surface_tension,
)


def test_frumkin_theta_zero_a_matches_langmuir_closed_form():
    for K, C in [(0.5, 2.0), (1.0, 10.0), (0.05, 100.0)]:
        theta = frumkin_theta(C, K, a=0.0)
        langmuir_theta = (K * C) / (1.0 + K * C)
        assert theta == pytest.approx(langmuir_theta, abs=1e-6)


def test_frumkin_theta_zero_concentration_is_zero_coverage():
    assert frumkin_theta(0.0, K=1.0, a=0.5) == 0.0


def test_frumkin_theta_rejects_negative_inputs():
    with pytest.raises(ValueError):
        frumkin_theta(-1.0, K=1.0, a=0.0)
    with pytest.raises(ValueError):
        frumkin_theta(1.0, K=-1.0, a=0.0)


def test_frumkin_surface_tension_zero_a_matches_szyszkowski_exactly():
    """a=0 must reduce EXACTLY to the existing, literature-validated
    Szyszkowski code -- this is the real regression guard for the
    Frumkin implementation, not a numeric literature match (none was
    found for this specific isotherm)."""
    cases = [
        (2.0, 72.0, 4.0e-6, 0.5, "nonionic"),
        (5.0, 72.0, 3.0e-6, 0.2, "ionic_no_added_salt"),
        (0.1, 72.0, 5.0e-6, 1.0, "ionic_excess_electrolyte"),
    ]
    for concentration, gamma0, gamma_max, K, system_type in cases:
        frumkin = frumkin_surface_tension(concentration, gamma0, gamma_max, K, a=0.0, system_type=system_type)
        szyszkowski = szyszkowski_surface_tension(concentration, gamma0, gamma_max, K, system_type=system_type)
        assert frumkin == pytest.approx(szyszkowski, abs=1e-9)


def test_frumkin_surface_tension_at_zero_concentration_is_gamma0():
    gamma = frumkin_surface_tension(0.0, 72.0, 4.0e-6, K=0.5, a=0.3, system_type="nonionic")
    assert gamma == pytest.approx(72.0)


def test_frumkin_attraction_gives_higher_coverage_than_langmuir():
    """a > 0 (net attraction) must give a MORE saturated monolayer at
    the same bulk concentration than the a=0 Langmuir case -- attractive
    lateral interactions favor closer packing, steepening the isotherm.
    This is the qualitative physical direction stated by both sources
    checked before implementing (see frumkin_theta's docstring)."""
    K, C = 0.3, 5.0
    theta_ideal = frumkin_theta(C, K, a=0.0)
    theta_attractive = frumkin_theta(C, K, a=1.0)
    assert theta_attractive > theta_ideal


def test_frumkin_repulsion_gives_lower_coverage_than_langmuir():
    K, C = 0.3, 5.0
    theta_ideal = frumkin_theta(C, K, a=0.0)
    theta_repulsive = frumkin_theta(C, K, a=-1.0)
    assert theta_repulsive < theta_ideal


def test_frumkin_theta_is_a_real_root_of_the_isotherm_equation():
    """Round-trip check independent of any external source: whatever
    theta is returned must actually satisfy the defining implicit
    equation K*C = theta/(1-theta) * exp(-2*a*theta)."""
    for K, a, C in [(0.4, 0.8, 3.0), (0.4, -0.8, 3.0), (2.0, 1.5, 0.5)]:
        theta = frumkin_theta(C, K, a)
        lhs = K * C
        rhs = (theta / (1.0 - theta)) * math.exp(-2.0 * a * theta)
        assert lhs == pytest.approx(rhs, rel=1e-6)


def test_frumkin_surface_tension_decreases_monotonically_with_concentration():
    """Below the condensation regime (|a| < 2), surface tension must
    decrease monotonically as concentration increases -- basic physical
    sanity, same property already asserted implicitly by the existing
    Szyszkowski behavior."""
    concentrations = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    gammas = [
        frumkin_surface_tension(c, 72.0, 4.0e-6, K=0.5, a=0.5, system_type="nonionic")
        for c in concentrations
    ]
    assert gammas == sorted(gammas, reverse=True)


def test_frumkin_surface_tension_rejects_bad_system_type():
    with pytest.raises(ValueError):
        frumkin_surface_tension(1.0, 72.0, 4.0e-6, K=0.5, a=0.0, system_type="not_a_real_type")
