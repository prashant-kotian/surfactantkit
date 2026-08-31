"""Tests for surfactantkit.hlb and surfactantkit.cpp.

Constants (Davies group numbers, Tanford formula coefficients, CPP
morphology thresholds) were checked against cited sources before being
hardcoded -- see the module docstrings in hlb.py and cpp.py. Tests here
verify: (1) exact boundary conditions that follow directly from the
formula definitions, (2) a worked example (SDS via Davies' method) built
purely from the verified group numbers, and (3) the CPP classification
thresholds at and around each boundary.
"""

import pytest

from surfactantkit.hlb import hlb_griffin, hlb_davies
from surfactantkit.cpp import (
    tanford_tail_volume,
    tanford_critical_length,
    aggregation_number_spherical,
    critical_packing_parameter,
    classify_aggregate_morphology,
)


def test_hlb_griffin_boundary_conditions():
    # fully hydrophilic -> HLB = 20; fully lipophilic -> HLB = 0
    assert hlb_griffin(100.0, 100.0) == pytest.approx(20.0)
    assert hlb_griffin(0.0, 100.0) == pytest.approx(0.0)
    assert hlb_griffin(50.0, 100.0) == pytest.approx(10.0)


def test_hlb_griffin_rejects_bad_input():
    with pytest.raises(ValueError):
        hlb_griffin(150.0, 100.0)  # hydrophilic > total
    with pytest.raises(ValueError):
        hlb_griffin(10.0, 0.0)


def test_hlb_davies_sds_worked_example():
    """SDS = CH3-(CH2)11-SO4Na: 1 CH3 + 11 CH2 + 1 SO4Na headgroup.
    HLB = 7 + 38.6 - 12*0.475 = 39.9, consistent with the commonly-cited
    ballpark of ~40 for SDS (not asserted as an independently-verified
    literature value here -- this checks the summation logic against
    the group numbers already verified in hlb.py)."""
    hlb = hlb_davies({"SO4Na": 1, "CH2": 11, "CH3": 1})
    assert hlb == pytest.approx(39.9, abs=0.01)


def test_hlb_davies_unknown_group_raises_instead_of_guessing():
    """Quaternary ammonium has no verified Davies number in this table
    (see hlb.py docstring) -- must raise, not silently return a guess."""
    with pytest.raises(KeyError):
        hlb_davies({"quaternary_ammonium": 1, "CH2": 11})


def test_tanford_tail_volume_matches_additive_formula():
    """Cross-check the shorthand (27.4 + 26.9*nc) against the more
    detailed per-group formula (v_CH3 + (nc-1)*v_CH2) for a C12 chain --
    they must agree exactly, since the shorthand is algebraically
    derived from the detailed formula, not an independent claim."""
    n_carbons = 12
    shorthand = 27.4 + 26.9 * n_carbons
    detailed = 54.3 + (n_carbons - 1) * 26.9
    assert shorthand == pytest.approx(detailed, abs=1e-9)
    assert tanford_tail_volume(n_carbons) == pytest.approx(detailed, abs=1e-9)


def test_tanford_critical_length_c12():
    # lc = 1.5 + 1.265*12 = 16.68 Angstrom -- in the commonly-cited
    # ballpark (~16.7 A) for a fully-extended dodecyl chain.
    assert tanford_critical_length(12) == pytest.approx(16.68, abs=0.01)


def test_aggregation_number_spherical_c12_lands_in_literature_range():
    """C12 chain (SDS/DTAB scale) geometric aggregation number should
    land near the commonly-reported literature range (~55-70) for these
    surfactants -- a sanity/order-of-magnitude check, not an exact
    literature match (this is a geometric estimate, not a substitute for
    a directly-measured aggregation number; see cpp.py docstring)."""
    v = tanford_tail_volume(12)
    lc = tanford_critical_length(12)
    nagg = aggregation_number_spherical(v, lc)
    assert 40.0 < nagg < 80.0


def test_aggregation_number_spherical_matches_direct_geometry_formula():
    import math

    v, r = 350.2, 16.68
    expected = ((4.0 / 3.0) * math.pi * r ** 3) / v
    assert aggregation_number_spherical(v, r) == pytest.approx(expected)


def test_aggregation_number_spherical_rejects_nonpositive():
    with pytest.raises(ValueError):
        aggregation_number_spherical(0, 16.68)
    with pytest.raises(ValueError):
        aggregation_number_spherical(350.2, -1.0)


def test_tanford_rejects_bad_chain_length():
    with pytest.raises(ValueError):
        tanford_tail_volume(0)
    with pytest.raises(ValueError):
        tanford_critical_length(-1)


def test_critical_packing_parameter_basic():
    # v=350.2 (C12), a0=50 A^2 (typical ionic headgroup area), lc=16.68 A
    cpp = critical_packing_parameter(350.2, 50.0, 16.68)
    assert cpp == pytest.approx(350.2 / (50.0 * 16.68))


def test_critical_packing_parameter_rejects_nonpositive():
    with pytest.raises(ValueError):
        critical_packing_parameter(0, 50.0, 16.68)
    with pytest.raises(ValueError):
        critical_packing_parameter(350.2, -1.0, 16.68)


@pytest.mark.parametrize(
    "cpp,expected",
    [
        (0.20, "spherical micelle"),
        (1.0 / 3.0, "spherical micelle"),          # boundary, inclusive
        (0.34, "cylindrical/rodlike micelle"),
        (0.5, "cylindrical/rodlike micelle"),       # boundary, inclusive
        (0.51, "vesicle/bilayer"),
        (1.0, "vesicle/bilayer"),                   # boundary, inclusive
        (1.01, "inverted structure"),
        (2.0, "inverted structure"),
    ],
)
def test_classify_aggregate_morphology_thresholds(cpp, expected):
    assert classify_aggregate_morphology(cpp) == expected


def test_classify_aggregate_morphology_rejects_nonpositive():
    with pytest.raises(ValueError):
        classify_aggregate_morphology(0)
    with pytest.raises(ValueError):
        classify_aggregate_morphology(-0.5)
