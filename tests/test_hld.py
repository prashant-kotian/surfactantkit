"""Tests for surfactantkit.hld -- the Hydrophilic-Lipophilic Difference
framework, added specifically to fill the gap left by Davies (1957) never
publishing a quaternary ammonium group number for the classic HLB scale
(see hlb.py and benchmark/METHOD_ALTERNATIVES_LITERATURE_REVIEW.md).

Source: Schirone, Tartaro, Gentile & Palazzo, JCIS Open 4 (2021) 100033,
doi:10.1016/j.jciso.2021.100033, read in full (PDF provided by the user
after automated fetch attempts were blocked -- see hld.py's module
docstring).
"""

import math

import pytest

from surfactantkit.hld import (
    CATIONIC_QUAT_CC,
    CATIONIC_QUAT_HLB_DAVIES,
    CATIONIC_QUAT_K,
    cc_mixing_rule,
    fit_k_and_cc_from_salinity_scan,
    hld_ionic,
    hld_nonionic,
    hlb_from_cc_cationic_quat,
    optimal_salinity_ionic,
)


def test_hld_ionic_zero_at_optimal_salinity():
    """Mathematically guaranteed round-trip: optimal_salinity_ionic is
    the direct algebraic inverse of hld_ionic at HLD=0."""
    k, eacn, cc = CATIONIC_QUAT_K, 8.0, CATIONIC_QUAT_CC["CTAB"][0]
    s_star = optimal_salinity_ionic(k, eacn, cc)
    assert hld_ionic(s_star, k, eacn, cc) == pytest.approx(0.0, abs=1e-9)


def test_hld_ionic_round_trip_with_nonzero_temperature_and_alpha():
    k, eacn, cc, alpha, delta_t = 0.7, 10.0, -5.0, 0.02, 5.0
    s_star = optimal_salinity_ionic(k, eacn, cc, alpha=alpha, delta_T_K=delta_t)
    assert hld_ionic(s_star, k, eacn, cc, alpha=alpha, delta_T_K=delta_t) == pytest.approx(0.0, abs=1e-9)


def test_hld_ionic_rejects_nonpositive_salinity():
    with pytest.raises(ValueError):
        hld_ionic(0.0, CATIONIC_QUAT_K, 8.0, -5.0)
    with pytest.raises(ValueError):
        hld_ionic(-1.0, CATIONIC_QUAT_K, 8.0, -5.0)


def test_hld_ionic_matches_manual_formula():
    """Independent recomputation in the test itself, not copy-pasted
    from the implementation."""
    salinity, k, eacn, cc, alpha, delta_t = 10.0, 0.7, 8.0, -5.7, 0.01, 3.0
    expected = math.log(salinity) - k * eacn - alpha * delta_t + cc
    assert hld_ionic(salinity, k, eacn, cc, alpha=alpha, delta_T_K=delta_t) == pytest.approx(expected)


def test_hld_nonionic_matches_manual_formula():
    salinity, b, k, eacn, cc, alpha, delta_t = 5.0, 0.13, 0.15, 6.0, 2.0, -0.06, -2.0
    expected = b * salinity - k * eacn - alpha * delta_t + cc
    assert hld_nonionic(salinity, b, k, eacn, cc, alpha=alpha, delta_T_K=delta_t) == pytest.approx(expected)


def test_cc_mixing_rule_reduces_to_pure_component_at_boundaries():
    """Mathematically guaranteed: at chi=[1,0] or [0,1] the mixing rule
    must exactly recover the pure-component Cc."""
    cc_ddab, cc_ctab = 8.3, -6.5
    assert cc_mixing_rule([cc_ddab, cc_ctab], [1.0, 0.0]) == pytest.approx(cc_ddab)
    assert cc_mixing_rule([cc_ddab, cc_ctab], [0.0, 1.0]) == pytest.approx(cc_ctab)


def test_cc_mixing_rule_matches_manual_weighted_average():
    cc_values = [8.3, -6.5, 0.1]
    chis = [0.5, 0.3, 0.2]
    expected = sum(c * x for c, x in zip(cc_values, chis))
    assert cc_mixing_rule(cc_values, chis) == pytest.approx(expected)


def test_cc_mixing_rule_rejects_bad_inputs():
    with pytest.raises(ValueError):
        cc_mixing_rule([1.0, 2.0], [0.5])  # mismatched lengths
    with pytest.raises(ValueError):
        cc_mixing_rule([1.0, 2.0], [0.3, 0.3])  # doesn't sum to 1
    with pytest.raises(ValueError):
        cc_mixing_rule([], [])


# --- Real numeric validation against the source paper's own Table 1 / Fig. 8 -


def test_hlb_from_cc_matches_paper_own_reported_pairs():
    """Real cross-check: for the surfactants where the source paper
    reports BOTH an independently-measured Cc (HLD-titration) and a
    Davies-scale HLB (from Proverbio et al. 2003), the fitted
    hlb_from_cc_cationic_quat must reproduce the real reported HLB
    reasonably closely -- this is the actual real-data regression the
    source paper fit, not an independent claim."""
    for name, real_hlb in CATIONIC_QUAT_HLB_DAVIES.items():
        cc, _ = CATIONIC_QUAT_CC[name]
        predicted_hlb = hlb_from_cc_cationic_quat(cc)
        assert predicted_hlb == pytest.approx(real_hlb, abs=1.5)


def test_cc_over_k_ratios_match_paper_table_1():
    """Internal-consistency check: Table 1's own reported Cc/k column
    should match Cc/CATIONIC_QUAT_K for each surfactant (both taken
    directly from the same real table) -- catches a transcription
    error in either the Cc or k constants."""
    reported_cc_over_k = {
        "LTAB": -14.0, "MTAB": -12.0, "CTAB": -8.14,
        "CMIC": -10.0, "BDHC": 0.14, "DDAB": 12.0,
    }
    for name, expected_ratio in reported_cc_over_k.items():
        cc, _ = CATIONIC_QUAT_CC[name]
        assert cc / CATIONIC_QUAT_K == pytest.approx(expected_ratio, abs=0.5)


def test_cationic_quat_hlb_davies_matches_direct_group_contribution_computation():
    """Cross-module consistency check: hld.py's CATIONIC_QUAT_HLB_DAVIES
    (sourced via Schirone et al. 2021's reproduction of Proverbio et al.
    2003) must match hlb.py's hlb_davies() computed DIRECTLY from the
    real quaternary-ammonium group numbers resolved 2026-09-08 (also
    sourced from Proverbio et al. 2003, independently). Two different
    routes to the same real data landing on the same numbers is a real,
    strong consistency check, not a tautology -- these were implemented
    in separate modules from separate readings of separate sources."""
    from surfactantkit.hlb import hlb_davies

    single_tail_trimethyl_carbons = {"LTAB": 12, "MTAB": 14, "CTAB": 16}
    for name, n_carbons in single_tail_trimethyl_carbons.items():
        computed = hlb_davies({"N_quat_trimethyl": 1, "CH2": n_carbons - 1, "CH3": 1})
        assert computed == pytest.approx(CATIONIC_QUAT_HLB_DAVIES[name], abs=0.01)

    computed_ddab = hlb_davies({"N_quat_dimethyl_dialkyl": 1, "CH2": 22, "CH3": 2})
    assert computed_ddab == pytest.approx(CATIONIC_QUAT_HLB_DAVIES["DDAB"], abs=0.01)


# --- fit_k_and_cc_from_salinity_scan ----------------------------------------
# Source paper's own Eq. 5 (this module's optimal_salinity_ionic), inverted:
# ln(S*) = k*EACN + alpha*deltaT - Cc is linear in EACN -- the real,
# standard "salinity scan across several oils" method that produces k/Cc
# from raw data. Validated here via a mathematically-guaranteed round trip
# (construct S* from a chosen true k/Cc via optimal_salinity_ionic itself,
# fit them back) -- no external multi-EACN raw dataset was sourced this
# session, same disclosed pattern as the other curve/regression tools built
# this session with no clean external numeric example.


def test_salinity_scan_fit_recovers_true_k_and_cc_round_trip():
    true_k, true_cc = CATIONIC_QUAT_K, CATIONIC_QUAT_CC["CTAB"][0]
    eacn_series = [6.0, 8.0, 10.0, 12.0, 14.0]
    s_star_series = [optimal_salinity_ionic(true_k, eacn, true_cc) for eacn in eacn_series]

    result = fit_k_and_cc_from_salinity_scan(eacn_series, s_star_series)
    assert result.k == pytest.approx(true_k, abs=1e-6)
    assert result.cc == pytest.approx(true_cc, abs=1e-6)
    assert result.r_squared == pytest.approx(1.0, abs=1e-6)
    assert result.n_points == 5


def test_salinity_scan_fit_round_trip_with_nonzero_temperature():
    true_k, true_cc, alpha, delta_t = 0.16, -3.0, 0.02, 4.0
    eacn_series = [4.0, 7.0, 9.0, 11.0]
    s_star_series = [optimal_salinity_ionic(true_k, eacn, true_cc, alpha=alpha, delta_T_K=delta_t) for eacn in eacn_series]

    result = fit_k_and_cc_from_salinity_scan(eacn_series, s_star_series, alpha=alpha, delta_T_K=delta_t)
    assert result.k == pytest.approx(true_k, abs=1e-6)
    assert result.cc == pytest.approx(true_cc, abs=1e-6)


def test_salinity_scan_fit_matches_manual_ols():
    """Independent recomputation of the OLS slope/intercept in the test
    itself, not copy-pasted from the implementation."""
    eacn_series = [5.0, 8.0, 11.0, 13.0]
    s_star_series = [3.0, 5.5, 9.0, 14.0]
    ln_s = [math.log(s) for s in s_star_series]

    n = len(eacn_series)
    mean_x = sum(eacn_series) / n
    mean_y = sum(ln_s) / n
    sxx = sum((x - mean_x) ** 2 for x in eacn_series)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(eacn_series, ln_s))
    expected_k = sxy / sxx
    expected_intercept = mean_y - expected_k * mean_x
    expected_cc = -expected_intercept  # alpha*deltaT=0 default

    result = fit_k_and_cc_from_salinity_scan(eacn_series, s_star_series)
    assert result.k == pytest.approx(expected_k)
    assert result.cc == pytest.approx(expected_cc)


def test_salinity_scan_fit_rejects_bad_inputs():
    with pytest.raises(ValueError):
        fit_k_and_cc_from_salinity_scan([8.0], [5.0])  # too few points
    with pytest.raises(ValueError):
        fit_k_and_cc_from_salinity_scan([6.0, 8.0], [5.0])  # mismatched lengths
    with pytest.raises(ValueError):
        fit_k_and_cc_from_salinity_scan([6.0, 8.0], [5.0, -1.0])  # non-positive salinity
    with pytest.raises(ValueError):
        fit_k_and_cc_from_salinity_scan([8.0, 8.0, 8.0], [5.0, 6.0, 7.0])  # all-identical EACN


def test_cationic_quat_cc_values_are_real_and_ordered_by_hydrophilicity():
    """Real physical/chemical sanity check from the source paper's own
    discussion: the single-tailed quats' Cc should become less negative
    (more hydrophobic-leaning) roughly with chain length increasing from
    C12 (LTAB) to C14 (MTAB) to C16 (CTAB), and DDAB (double-tailed,
    genuinely lipophilic/reverse-micelle-forming) should be positive."""
    assert CATIONIC_QUAT_CC["LTAB"][0] < CATIONIC_QUAT_CC["MTAB"][0] < CATIONIC_QUAT_CC["CTAB"][0]
    assert CATIONIC_QUAT_CC["DDAB"][0] > 0
    assert CATIONIC_QUAT_CC["LTAB"][0] < 0
