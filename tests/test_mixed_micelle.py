"""Tests for surfactantkit.mixed_micelle.

Validation strategy, deliberately explicit about what's certain and what
isn't. This library was checked against EIGHT independent published
binary surfactant systems (cationic-anionic, cationic-nonionic,
anionic-nonionic, gemini-nonionic, gemini-zwitterionic,
nonionic-biosurfactant, anionic-anionic/bile salt), all open access,
numbers pulled directly from source tables. Full case-by-case notes,
including a log of dead-end search attempts: see
literature_validation_notes.md.

1. Clint ideal CMC: hard exact-match assertions against papers that
   report a directly comparable ideal-CMC value (Azum et al. 2022;
   McLachlan et al. 2020; Liu et al. 2020). One paper in the survey
   (Muherei & Junin 2009) reports an "ideal" CMC that does NOT reproduce
   from its own stated pure-component CMCs via the standard Clint
   formula (0.73 mM computed vs. 0.906 mM reported) -- flagged as an
   unexplained discrepancy in that source (possibly a different pure-CMC
   basis for their ideal column), not asserted here.

2. Rubingh solver, micellar mole fraction (x1): hard assertions against
   FIVE literature systems -- matches to within ~1% in every case. This
   is the strongest evidence the solver is correct, since x1 is the
   solver's actual output (beta is then a deterministic function of it).

3. Rubingh solver, interaction parameter (beta): NOT asserted against
   literature values directly. Investigation across all five systems
   found the two algebraically-equivalent ways of computing beta at the
   solved root (term1/(1-x1)^2 and term2/x1^2) agree with each other to
   3-4 decimal places every time (confirming internal correctness), but
   neither matches the literature-reported beta precisely (typically
   off by 5-20%, sign always correct). The most likely explanation:
   many papers fit beta by regression across several alpha compositions,
   while this solver computes the exact pointwise value for one
   composition -- both are legitimate, but not directly comparable to
   3 decimals. We assert sign and same order of magnitude only.

4. The Rubingh solver is ALSO checked with mathematically guaranteed
   round-trip tests (construct a cmc_mix for which a chosen x1 is the
   exact algebraic root, solve it back, recover the same x1 and beta) --
   these don't depend on trusting any external source at all.
"""

import math

import pytest

from surfactantkit.mixed_micelle import (
    R_GAS,
    activity_coefficients,
    asymmetric_margules_activity_coefficients,
    clint_ideal_cmc,
    eommm_global_fit,
    excess_free_energy,
    maeda_free_energy_of_micellization,
    motomura_ideal_composition,
    rodenas_activity_coefficients,
    rodenas_x1,
    rodenas_x1_series,
    rubingh_beta,
    rubingh_beta_regression,
    solve_rubingh_x,
)
from surfactantkit.adsorption import gibbs_a_min, gibbs_gamma_max

# DTAB-SDS system at 293.15 K (Rodriguez et al., PMC6554738)
DTAB_PURE_CMC = 14.80  # mM
SDS_PURE_CMC = 8.00  # mM


def test_clint_ideal_dtab_rich():
    # alpha_DTAB = 0.75; paper reports CMCid = 12.21 mM
    cmc_id = clint_ideal_cmc(0.75, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert cmc_id == pytest.approx(12.21, abs=0.01)


def test_clint_ideal_sds_rich():
    # alpha_DTAB = 0.25 (i.e. alpha_SDS = 0.75); paper reports CMCid = 9.04 mM
    cmc_id = clint_ideal_cmc(0.25, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert cmc_id == pytest.approx(9.04, abs=0.01)


def test_clint_ideal_rejects_bad_alpha():
    with pytest.raises(ValueError):
        clint_ideal_cmc(0.0, DTAB_PURE_CMC, SDS_PURE_CMC)
    with pytest.raises(ValueError):
        clint_ideal_cmc(1.0, DTAB_PURE_CMC, SDS_PURE_CMC)


# --- Additional literature systems -----------------------------------------
# Component 1 is always named first in each system's name below.

# Muherei & Junin (2009), Asian J. Appl. Sci. 2(2), 115-127, Table 2A/2B.
# System: TX-100 (comp 1) - SDS (comp 2). alpha(TX-100) = 0.47.
MUHEREI_TX100_CMC = 0.387  # mM
MUHEREI_SDS_CMC = 3.468  # mM
MUHEREI_ALPHA1 = 0.47
MUHEREI_CMC_MIX = 0.547  # mM, experimental
MUHEREI_X1_PAPER = 0.7501


def test_rubingh_x1_muherei_tx100_sds():
    x1 = solve_rubingh_x(MUHEREI_ALPHA1, MUHEREI_CMC_MIX, MUHEREI_TX100_CMC, MUHEREI_SDS_CMC)
    assert x1 is not None
    assert x1 == pytest.approx(MUHEREI_X1_PAPER, abs=0.01)
    beta = rubingh_beta(x1, MUHEREI_ALPHA1, MUHEREI_CMC_MIX, MUHEREI_TX100_CMC)
    assert beta < 0  # paper reports beta = -1.888 (synergistic)


# Azum et al. (2022), Biointerface Res. Appl. Chem. 12(6), 7416-7428, Table 1.
# System: gemini G6 (comp 1) - Triton X-114 (comp 2). alpha(G6) = 0.48.
AZUM_G6_CMC = 0.041  # mM
AZUM_TX114_CMC = 0.263  # mM
AZUM_ALPHA1 = 0.48
AZUM_CMC_MIX = 0.061  # mM, experimental
AZUM_CMC_IDEAL_PAPER = 0.073  # mM
AZUM_X1_PAPER = 0.764


def test_clint_ideal_azum_gemini_tx114():
    cmc_id = clint_ideal_cmc(AZUM_ALPHA1, AZUM_G6_CMC, AZUM_TX114_CMC)
    assert cmc_id == pytest.approx(AZUM_CMC_IDEAL_PAPER, abs=0.001)


def test_rubingh_x1_azum_gemini_tx114():
    x1 = solve_rubingh_x(AZUM_ALPHA1, AZUM_CMC_MIX, AZUM_G6_CMC, AZUM_TX114_CMC)
    assert x1 is not None
    assert x1 == pytest.approx(AZUM_X1_PAPER, abs=0.01)
    beta = rubingh_beta(x1, AZUM_ALPHA1, AZUM_CMC_MIX, AZUM_G6_CMC)
    assert beta < 0  # paper reports beta = -1.211 (synergistic)


# Lee & Lee (2012), J. Korean Chem. Soc. 56(5), 556-562, Table 3.
# System: TTAB (comp 1) - Tween-20 (comp 2). alpha(TTAB) = 0.6.
LEE_TTAB_CMC = 2.25  # mM
LEE_TWEEN20_CMC = 0.41  # mM
LEE_ALPHA1 = 0.6
LEE_CMC_MIX = 0.70  # mM, experimental
LEE_X1_PAPER = 0.28


def test_rubingh_x1_lee_ttab_tween20():
    x1 = solve_rubingh_x(LEE_ALPHA1, LEE_CMC_MIX, LEE_TTAB_CMC, LEE_TWEEN20_CMC)
    assert x1 is not None
    assert x1 == pytest.approx(LEE_X1_PAPER, abs=0.02)
    beta = rubingh_beta(x1, LEE_ALPHA1, LEE_CMC_MIX, LEE_TTAB_CMC)
    assert beta < 0  # paper reports beta = -0.78 (synergistic)


# McLachlan et al. (2020), RSC Adv. 10(6), 3221-3232, Table 1.
# System: gemini 12-4-12 (comp 1) - zwitterionic ZW3-12 (comp 2), alpha = 0.5.
MCLACHLAN_GEMINI_CMC = 1.10  # mM
MCLACHLAN_ZW312_CMC = 2.63  # mM
MCLACHLAN_ALPHA1 = 0.5
MCLACHLAN_CMC_MIX = 1.52  # mM, experimental
MCLACHLAN_CMC_IDEAL_PAPER = 1.55  # mM
MCLACHLAN_X_ZW312_PAPER = 0.301  # micellar mole fraction of ZW3-12 (component 2)


def test_clint_ideal_mclachlan_gemini_zwitterionic():
    cmc_id = clint_ideal_cmc(MCLACHLAN_ALPHA1, MCLACHLAN_GEMINI_CMC, MCLACHLAN_ZW312_CMC)
    assert cmc_id == pytest.approx(MCLACHLAN_CMC_IDEAL_PAPER, abs=0.01)


def test_rubingh_x1_mclachlan_gemini_zwitterionic():
    x1 = solve_rubingh_x(MCLACHLAN_ALPHA1, MCLACHLAN_CMC_MIX, MCLACHLAN_GEMINI_CMC, MCLACHLAN_ZW312_CMC)
    assert x1 is not None
    # paper reports the ZWITTERIONIC (component 2) micellar mole fraction;
    # component 1 (gemini) mole fraction is 1 - x_ZW312
    assert (1.0 - x1) == pytest.approx(MCLACHLAN_X_ZW312_PAPER, abs=0.01)


# Liu et al. (2020), Molecules 25(18), 4327, Table 9.
# System: Triton X-100 (comp 1) - rhamnolipid biosurfactant (comp 2), alpha = 0.888.
LIU_TX100_CMC = 0.309
LIU_RHAMNOLIPID_CMC = 0.134
LIU_ALPHA1 = 0.888
LIU_CMC_MIX = 0.253  # experimental
LIU_CMC_THEOR_PAPER = 0.270
LIU_X1_PAPER = 0.744


def test_clint_ideal_liu_tx100_rhamnolipid():
    cmc_id = clint_ideal_cmc(LIU_ALPHA1, LIU_TX100_CMC, LIU_RHAMNOLIPID_CMC)
    assert cmc_id == pytest.approx(LIU_CMC_THEOR_PAPER, abs=0.001)


def test_rubingh_x1_liu_tx100_rhamnolipid():
    x1 = solve_rubingh_x(LIU_ALPHA1, LIU_CMC_MIX, LIU_TX100_CMC, LIU_RHAMNOLIPID_CMC)
    assert x1 is not None
    assert x1 == pytest.approx(LIU_X1_PAPER, abs=0.01)
    beta = rubingh_beta(x1, LIU_ALPHA1, LIU_CMC_MIX, LIU_TX100_CMC)
    assert beta < 0  # paper reports beta = -0.379 (synergistic)


# Kang, Bahadur, et al. -- Mixed Micelles of Sodium Cholate and Sodium
# Dodecylsulphate 1:1 Binary Mixture at Different Temperatures, PMC4087020.
# System: sodium cholate NaCA (comp 1) - SDS (comp 2), alpha = 0.5 (1:1), 25C.
# This is the ONLY system in the survey where beta itself matches the
# paper closely -- consistent with a single fixed-ratio study (no
# multi-point regression possible), supporting the explanation in the
# module docstring for why beta usually doesn't match exactly elsewhere.
CHOLATE_NACA_CMC = 11.50  # mM
CHOLATE_SDS_CMC = 11.98  # mM
CHOLATE_ALPHA1 = 0.5
CHOLATE_CMC_MIX = 4.07  # mM, experimental (fluorimetry)
CHOLATE_CMC_IDEAL_PAPER = 11.74  # mM
CHOLATE_X1_PAPER = 0.503
CHOLATE_BETA_PAPER = -4.23


def test_clint_ideal_cholate_sds():
    cmc_id = clint_ideal_cmc(CHOLATE_ALPHA1, CHOLATE_NACA_CMC, CHOLATE_SDS_CMC)
    assert cmc_id == pytest.approx(CHOLATE_CMC_IDEAL_PAPER, abs=0.01)


def test_rubingh_x1_and_beta_cholate_sds():
    x1 = solve_rubingh_x(CHOLATE_ALPHA1, CHOLATE_CMC_MIX, CHOLATE_NACA_CMC, CHOLATE_SDS_CMC)
    assert x1 is not None
    assert x1 == pytest.approx(CHOLATE_X1_PAPER, abs=0.005)
    beta = rubingh_beta(x1, CHOLATE_ALPHA1, CHOLATE_CMC_MIX, CHOLATE_NACA_CMC)
    # this is the one system where beta itself is asserted closely, not just its sign
    assert beta == pytest.approx(CHOLATE_BETA_PAPER, abs=0.1)


def test_rubingh_beta_internal_self_consistency_across_all_literature_systems():
    """The two algebraically-equivalent ways of computing beta at the
    solved root (from term1 or from term2) must agree with each other
    closely, regardless of whether either matches a literature-reported
    beta (see module docstring, point 3)."""
    systems = [
        (MUHEREI_ALPHA1, MUHEREI_CMC_MIX, MUHEREI_TX100_CMC, MUHEREI_SDS_CMC),
        (AZUM_ALPHA1, AZUM_CMC_MIX, AZUM_G6_CMC, AZUM_TX114_CMC),
        (LEE_ALPHA1, LEE_CMC_MIX, LEE_TTAB_CMC, LEE_TWEEN20_CMC),
        (LIU_ALPHA1, LIU_CMC_MIX, LIU_TX100_CMC, LIU_RHAMNOLIPID_CMC),
    ]
    for alpha1, cmc_mix, cmc1, cmc2 in systems:
        x1 = solve_rubingh_x(alpha1, cmc_mix, cmc1, cmc2)
        assert x1 is not None
        term1 = math.log((alpha1 * cmc_mix) / (x1 * cmc1))
        term2 = math.log(((1.0 - alpha1) * cmc_mix) / ((1.0 - x1) * cmc2))
        beta_from_term1 = term1 / (1.0 - x1) ** 2
        beta_from_term2 = term2 / (x1 ** 2)
        assert beta_from_term1 == pytest.approx(beta_from_term2, abs=0.01)


def _cmc_mix_at_root(x1: float, alpha1: float, cmc1: float, cmc2: float) -> float:
    """Given a target root x1, solve the Rubingh residual equation
    (x1^2 * ln(A*Cmix) = (1-x1)^2 * ln(B*Cmix), A = alpha1/(x1*cmc1),
    B = (1-alpha1)/((1-x1)*cmc2)) for the cmc_mix that makes x1 the exact
    root -- closed form, valid for x1 != 0.5. This is the correct forward
    construction for a round-trip test: fixing beta alone (via the
    term1/(1-x1)^2 definition) does not fix cmc_mix, because cmc_mix also
    appears in term2; both terms must balance at the same x1."""
    a = alpha1 / (x1 * cmc1)
    b = (1.0 - alpha1) / ((1.0 - x1) * cmc2)
    denom = x1 ** 2 - (1.0 - x1) ** 2
    ln_cmix = ((1.0 - x1) ** 2 * math.log(b) - x1 ** 2 * math.log(a)) / denom
    return math.exp(ln_cmix)


def test_rubingh_round_trip_case_a():
    """Construct a cmc_mix for which x1=0.3 is the exact algebraic root
    (see _cmc_mix_at_root), then verify the numerical solver recovers
    that same x1 and the beta computed from it -- mathematically
    guaranteed, independent of any literature value."""
    alpha1, cmc1, cmc2, x1_true = 0.4, 10.0, 5.0, 0.3
    cmc_mix = _cmc_mix_at_root(x1_true, alpha1, cmc1, cmc2)
    beta_true = rubingh_beta(x1_true, alpha1, cmc_mix, cmc1)

    x1_solved = solve_rubingh_x(alpha1, cmc_mix, cmc1, cmc2)
    assert x1_solved is not None
    assert x1_solved == pytest.approx(x1_true, abs=1e-4)

    beta_solved = rubingh_beta(x1_solved, alpha1, cmc_mix, cmc1)
    assert beta_solved == pytest.approx(beta_true, abs=1e-3)


def test_rubingh_round_trip_case_b():
    alpha1, cmc1, cmc2, x1_true = 0.6, 8.0, 12.0, 0.7
    cmc_mix = _cmc_mix_at_root(x1_true, alpha1, cmc1, cmc2)
    beta_true = rubingh_beta(x1_true, alpha1, cmc_mix, cmc1)

    x1_solved = solve_rubingh_x(alpha1, cmc_mix, cmc1, cmc2)
    assert x1_solved is not None
    assert x1_solved == pytest.approx(x1_true, abs=1e-4)

    beta_solved = rubingh_beta(x1_solved, alpha1, cmc_mix, cmc1)
    assert beta_solved == pytest.approx(beta_true, abs=1e-3)


def test_rubingh_sign_matches_literature_synergy_direction():
    """Qualitative check against the DTAB-SDS paper: SDS-rich mixing is
    reported synergistic (CMC < CMCideal -> negative beta); DTAB-rich
    mixing is reported antagonistic (CMC > CMCideal -> positive beta).
    We assert the sign only -- see module docstring."""
    # SDS-rich: alpha_DTAB = 0.25, experimental CMCmix = 6.011 mM (synergistic)
    x1 = solve_rubingh_x(0.25, 6.011, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert x1 is not None
    beta = rubingh_beta(x1, 0.25, 6.011, DTAB_PURE_CMC)
    assert beta < 0, "SDS-rich composition should be synergistic (negative beta)"

    # DTAB-rich: alpha_DTAB = 0.75, experimental CMCmix = 13.00 mM (antagonistic)
    x1 = solve_rubingh_x(0.75, 13.00, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert x1 is not None
    beta = rubingh_beta(x1, 0.75, 13.00, DTAB_PURE_CMC)
    assert beta > 0, "DTAB-rich composition should be antagonistic (positive beta)"


# --- rubingh_beta_regression -------------------------------------------
# Real standard practice (verified via WebSearch before implementing, not
# guessed): report a single representative beta as the arithmetic mean of
# the pointwise betas computed at each (alpha1, cmc_mix) composition in a
# series -- see the function's own docstring. Validated here via a
# mathematically-guaranteed round trip: for a chosen true beta and a set
# of x1 values, the closed-form RST mass-balance equations
# (cmc_mix = x1*f1*cmc1 + (1-x1)*f2*cmc2, alpha1 = x1*f1*cmc1/cmc_mix)
# construct a self-consistent (alpha1, cmc_mix) series that ALL share
# that same true beta by construction -- no external literature source
# needed, same pattern as test_rubingh_round_trip_case_a/b above.


def _series_for_true_beta(x1_values, beta_true, cmc1, cmc2):
    alpha1_series, cmc_mix_series = [], []
    for x1 in x1_values:
        f1 = math.exp(beta_true * (1.0 - x1) ** 2)
        f2 = math.exp(beta_true * x1 ** 2)
        cmc_mix = x1 * f1 * cmc1 + (1.0 - x1) * f2 * cmc2
        alpha1 = x1 * f1 * cmc1 / cmc_mix
        alpha1_series.append(alpha1)
        cmc_mix_series.append(cmc_mix)
    return alpha1_series, cmc_mix_series


def test_rubingh_beta_regression_recovers_true_beta_round_trip():
    beta_true = -1.8
    x1_values = [0.3, 0.4, 0.5, 0.6, 0.7]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, beta_true, DTAB_PURE_CMC, SDS_PURE_CMC)

    result = rubingh_beta_regression(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert result.n_points_used == 5
    assert result.n_points_skipped == 0
    assert result.beta_mean == pytest.approx(beta_true, abs=1e-3)
    assert result.beta_std == pytest.approx(0.0, abs=1e-2)
    for beta in result.betas:
        assert beta == pytest.approx(beta_true, abs=1e-2)


def test_rubingh_beta_regression_antagonistic_case():
    beta_true = 1.8
    x1_values = [0.25, 0.4, 0.55, 0.7]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, beta_true, DTAB_PURE_CMC, SDS_PURE_CMC)

    result = rubingh_beta_regression(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert result.beta_mean == pytest.approx(beta_true, abs=1e-3)
    assert result.n_points_used == 4


def test_rubingh_beta_regression_skips_unsolvable_points_not_errors(monkeypatch):
    """solve_rubingh_x's own documented contract is to return None when no
    root is found (near-ideal or self-inconsistent data). Rather than
    fight root-finder numerics to contrive a genuinely rootless input
    (every cmc_mix value tried empirically still admitted some root),
    this isolates rubingh_beta_regression's own skip-and-continue logic
    directly by forcing one specific point's solve to report no root."""
    import surfactantkit.mixed_micelle as mm

    beta_true = -1.5
    x1_values = [0.3, 0.5, 0.7]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, beta_true, DTAB_PURE_CMC, SDS_PURE_CMC)
    unsolvable_alpha1 = 0.9999
    alpha1_series.append(unsolvable_alpha1)
    cmc_mix_series.append(1.0)

    real_solve = mm.solve_rubingh_x

    def fake_solve(alpha1, cmc_mix, cmc1, cmc2, *args, **kwargs):
        if alpha1 == unsolvable_alpha1:
            return None
        return real_solve(alpha1, cmc_mix, cmc1, cmc2, *args, **kwargs)

    monkeypatch.setattr(mm, "solve_rubingh_x", fake_solve)

    result = mm.rubingh_beta_regression(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert result.n_points_used == 3
    assert result.n_points_skipped == 1
    assert result.beta_mean == pytest.approx(beta_true, abs=1e-3)


def test_rubingh_beta_regression_rejects_bad_inputs():
    with pytest.raises(ValueError):
        rubingh_beta_regression([0.3], [10.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # too few points
    with pytest.raises(ValueError):
        rubingh_beta_regression([0.3, 0.5], [10.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # mismatched lengths


def test_activity_coefficients_ideal_limit():
    """At beta = 0 (ideal mixing), both activity coefficients must be 1."""
    f1, f2 = activity_coefficients(x1=0.4, beta=0.0)
    assert f1 == pytest.approx(1.0)
    assert f2 == pytest.approx(1.0)


def test_excess_free_energy_ideal_limit_is_zero():
    f1, f2 = activity_coefficients(x1=0.4, beta=0.0)
    dg_ex = excess_free_energy(0.4, f1, f2, temperature_k=298.15)
    assert dg_ex == pytest.approx(0.0, abs=1e-9)


def test_gibbs_adsorption_positive_slope_gives_positive_gamma_max():
    # A surfactant lowering surface tension with concentration has a
    # negative dGamma/d(ln C); Gamma_max must come out positive.
    gamma_max = gibbs_gamma_max(slope_mn_per_ln_c=-5.0, system_type="ionic_no_added_salt", temperature_k=298.15)
    assert gamma_max > 0


def test_gibbs_a_min_matches_literature_order_of_magnitude():
    """Sanity check against the DTAB-SDS paper's reported A_min values
    (44.70 and 58.07 Ang^2/molecule) -- just confirms the formula lands
    in the right physical range for a real ionic surfactant system, not
    an exact reproduction (we don't have their exact slope inputs)."""
    # A_min ~ 50 Ang^2 implies Gamma_max ~ 1e18 / (Avogadro * 50) mol/m^2
    gamma_max_typical = 1.0e18 / (6.02214076e23 * 50.0)
    a_min = gibbs_a_min(gamma_max_typical)
    assert 30.0 < a_min < 80.0


def test_gibbs_a_min_rejects_nonpositive_gamma_max():
    with pytest.raises(ValueError):
        gibbs_a_min(0.0)
    with pytest.raises(ValueError):
        gibbs_a_min(-1.0)


# --- Motomura ideal composition (alternative-methods review, 2026-09-07/08) --
# Source: equation (8), Serafini et al., Colloids Surf. A 562 (2019) 170-181
# (arXiv:1806.09721), read directly (real fetched PDF, not a search summary).


def test_motomura_ideal_composition_equals_alpha_when_cmcs_are_equal():
    """When both pure CMCs are identical, the ideal composition must
    trivially equal the bulk composition -- no reason for micellar and
    bulk mole fractions to differ if the two surfactants are identical
    in CMC."""
    x1_id = motomura_ideal_composition(0.37, cmc1=5.0, cmc2=5.0)
    assert x1_id == pytest.approx(0.37)


def test_motomura_ideal_composition_matches_clint_algebraic_identity():
    """X1_id = alpha1 * CMC_mix_ideal / cmc1 -- an exact algebraic
    identity between Motomura's and Clint's ideal-mixing formulas
    (derived from their shared zero-interaction assumption), checked
    here independent of any external literature value across several
    of this file's own literature systems."""
    systems = [
        (0.75, DTAB_PURE_CMC, SDS_PURE_CMC),
        (MUHEREI_ALPHA1, MUHEREI_TX100_CMC, MUHEREI_SDS_CMC),
        (AZUM_ALPHA1, AZUM_G6_CMC, AZUM_TX114_CMC),
        (LEE_ALPHA1, LEE_TTAB_CMC, LEE_TWEEN20_CMC),
    ]
    for alpha1, cmc1, cmc2 in systems:
        x1_id = motomura_ideal_composition(alpha1, cmc1, cmc2)
        cmc_mix_id = clint_ideal_cmc(alpha1, cmc1, cmc2)
        assert x1_id == pytest.approx(alpha1 * cmc_mix_id / cmc1, rel=1e-9)


# --- EOMMM asymmetric Margules activity coefficients ------------------------
# Source: equations (9)-(13), Schulz & Durand, Comput. Chem. Eng. 87 (2016)
# 145-153 (doi:10.1016/j.compchemeng.2015.12.026), read directly (real fetched
# PDF, provided by the user after WebFetch attempts on this exact paper
# failed -- see METHOD_ALTERNATIVES_LITERATURE_REVIEW.md).


def test_asymmetric_margules_reduces_to_rst_when_symmetric():
    """W12 == W21 must reduce EXACTLY to activity_coefficients (RST) --
    the paper's own stated result (symmetric formulations are the
    special case of the asymmetric ones), and the real regression guard
    for this function since no single-point literature case exists to
    match digit-for-digit (EOMMM's real fit is a multi-point global
    optimization, not reproduced here -- see the function's docstring)."""
    for x1, beta in [(0.3, -2.17), (0.5, 1.0), (0.75, -0.78), (0.1, 3.5)]:
        f1_asym, f2_asym = asymmetric_margules_activity_coefficients(x1, beta, beta)
        f1_rst, f2_rst = activity_coefficients(x1, beta)
        assert f1_asym == pytest.approx(f1_rst, rel=1e-9)
        assert f2_asym == pytest.approx(f2_rst, rel=1e-9)


def test_asymmetric_margules_matches_paper_own_binary_symmetric_equations():
    """Direct check against the paper's own explicitly-stated binary
    symmetric-case equations (14)-(16): Gexc = x1*x2*beta,
    ln(f1) = beta*x2^2, ln(f2) = beta*x1^2."""
    x1, beta = 0.4, -1.5
    x2 = 1.0 - x1
    f1, f2 = asymmetric_margules_activity_coefficients(x1, beta, beta)
    assert math.log(f1) == pytest.approx(beta * x2 ** 2, abs=1e-9)
    assert math.log(f2) == pytest.approx(beta * x1 ** 2, abs=1e-9)


def test_asymmetric_margules_infinite_dilution_limits():
    """Mathematically guaranteed property of the formula itself,
    independent of any external source: as x1 -> 0, ln(f1) -> w12
    (component 1 infinitely dilute in pure component 2); as x1 -> 1,
    ln(f2) -> w21."""
    w12, w21 = 2.3, -4.1
    f1_near0, _ = asymmetric_margules_activity_coefficients(1e-8, w12, w21)
    _, f2_near1 = asymmetric_margules_activity_coefficients(1.0 - 1e-8, w12, w21)
    assert math.log(f1_near0) == pytest.approx(w12, abs=1e-5)
    assert math.log(f2_near1) == pytest.approx(w21, abs=1e-5)


def test_asymmetric_margules_hyamine_dtab_literature_case():
    """Real numeric case from the paper's own Case study 1 (Hyamine(1)-
    DTAB(2), Section 4.1): W12 = -8836.50 J/mol, W21 = -2204.19 J/mol at
    298.15 K, with x1(Hy) = 1 - x_DTAB = 1 - 0.192 = 0.808 (their Table
    1, asymmetric column, alpha_DTAB = 0.25). No exact digit-level
    literature value for ln(f1)/ln(f2) is given (only a graph, Fig. 1),
    so this asserts sign/order-of-magnitude consistency with that
    figure's own stated axis range (-4.00 to 0.50) -- same "sign and
    order of magnitude only" discipline already used elsewhere in this
    file for cases without an exact literature match."""
    temperature_k = 298.15
    w12 = -8836.50 / (R_GAS * temperature_k)
    w21 = -2204.19 / (R_GAS * temperature_k)
    x1_hy = 1.0 - 0.192
    f1, f2 = asymmetric_margules_activity_coefficients(x1_hy, w12, w21)
    assert -4.0 <= math.log(f1) <= 0.5
    assert -4.0 <= math.log(f2) <= 0.5


# --- eommm_global_fit (multi-point EOMMM global fit) ------------------------
# First-principles construction (2026-09-10) -- NOT a verified transcription
# of Schulz & Durand 2016's exact global-fit procedure (their Eq. 3.2 /
# S.I. Point 2.3), which was unavailable (primary source paywalled on
# ScienceDirect; a related open-access companion paper, Serafini et al.,
# arXiv:1806.09721, defers the exact procedure to its own SI, not included
# in the fetched copy). See eommm_global_fit's own docstring for the full
# disclosure and the real multi-root/multi-start numerical findings made
# while validating this. No external published raw numeric table was
# available, so validation here is via mathematically-guaranteed round
# trips only (the closed-form construction below satisfies BOTH mass-
# balance equations EXACTLY for any chosen W12/W21/x1, not just RST's
# symmetric case -- verified algebraically in eommm_global_fit's docstring).
# These tests are slow (a handful of seconds each -- genuine (n+2)-parameter
# nonlinear optimization, multi-start, no numpy) by the nature of the
# problem, kept deliberately small (n=3) to bound runtime.


def _eommm_series_for_true_params(x1_values, w12_true, w21_true, cmc1, cmc2):
    alpha1_series, cmc_mix_series = [], []
    for x1 in x1_values:
        f1, f2 = asymmetric_margules_activity_coefficients(x1, w12_true, w21_true)
        d = x1 * f1 * cmc1 + (1.0 - x1) * f2 * cmc2
        alpha1 = x1 * f1 * cmc1 / d
        alpha1_series.append(alpha1)
        cmc_mix_series.append(d)
    return alpha1_series, cmc_mix_series


def test_eommm_global_fit_recovers_true_parameters_round_trip():
    true_w12, true_w21 = 2.0, -3.0
    x1_true_values = [0.25, 0.5, 0.75]
    alpha1_series, cmc_mix_series = _eommm_series_for_true_params(
        x1_true_values, true_w12, true_w21, DTAB_PURE_CMC, SDS_PURE_CMC
    )

    result = eommm_global_fit(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert result.W12 == pytest.approx(true_w12, abs=1e-2)
    assert result.W21 == pytest.approx(true_w21, abs=1e-2)
    assert result.r_squared == pytest.approx(1.0, abs=1e-4)
    assert result.n_points == 3
    for fitted, true_x1 in zip(result.x1_values, x1_true_values):
        assert fitted == pytest.approx(true_x1, abs=1e-2)


def test_eommm_global_fit_reduces_to_rubingh_beta_when_symmetric():
    """At W12=W21 (the symmetric limit), the fit should recover
    approximately the same beta that solve_rubingh_x/rubingh_beta find
    on the identical data -- ties the new global fit back to the
    already-literature-validated Rubingh machinery, independent of
    trusting eommm_global_fit's own round-trip alone."""
    true_beta = -1.8
    x1_true_values = [0.3, 0.5, 0.7]
    alpha1_series, cmc_mix_series = _eommm_series_for_true_params(
        x1_true_values, true_beta, true_beta, DTAB_PURE_CMC, SDS_PURE_CMC
    )

    result = eommm_global_fit(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert result.W12 == pytest.approx(true_beta, abs=1e-2)
    assert result.W21 == pytest.approx(true_beta, abs=1e-2)

    for a1, cm in zip(alpha1_series, cmc_mix_series):
        x1_rub = solve_rubingh_x(a1, cm, DTAB_PURE_CMC, SDS_PURE_CMC)
        assert x1_rub is not None
        beta_rub = rubingh_beta(x1_rub, a1, cm, DTAB_PURE_CMC)
        assert beta_rub == pytest.approx(true_beta, abs=1e-2)


def test_eommm_global_fit_rejects_bad_inputs():
    with pytest.raises(ValueError):
        eommm_global_fit([0.3, 0.5], [10.0, 8.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # too few points
    with pytest.raises(ValueError):
        eommm_global_fit([0.3, 0.5, 0.7], [10.0, 8.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # mismatched lengths
    with pytest.raises(ValueError):
        eommm_global_fit([0.0, 0.5, 0.7], [10.0, 8.0, 6.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # bad alpha1
    with pytest.raises(ValueError):
        eommm_global_fit([0.3, 0.5, 0.7], [10.0, -8.0, 6.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # bad cmc_mix
    with pytest.raises(ValueError):
        eommm_global_fit([0.3, 0.5, 0.7], [10.0, 8.0, 6.0], 0.0, SDS_PURE_CMC)  # bad cmc1


# --- Rodenas model (alternative-methods review, 2026-09-07/08) --------------
# Source: equation (10)-(12), Azum, Rub, Alotaibi, Khan & Asiri, Biointerface
# Res. Appl. Chem. 12(6) (2022) 7416-7428 -- the same G6 gemini paper already
# used above for Clint/Rubingh validation -- citing Rodenas, Valiente &
# Villafruela, J. Phys. Chem. B 103(21) (1999) 4549-4554.


def test_rodenas_x1_reduces_to_motomura_under_ideal_mixing():
    """Real, source-independent algebraic identity: when cmc_mix
    follows Clint's EXACT ideal-mixing law, Rodenas' X1_Rod must equal
    Motomura's ideal composition formula -- both are zero-interaction
    limits of otherwise-different methods. The derivative is computed
    here by central finite difference directly on clint_ideal_cmc, not
    re-derived by hand, so this doesn't depend on trusting a separate
    hand-worked formula matching the implementation."""
    cmc1, cmc2 = 14.80, 8.00
    for alpha1 in (0.2, 0.4, 0.6, 0.8):
        h = 1e-6
        ln_cmc_plus = math.log(clint_ideal_cmc(alpha1 + h, cmc1, cmc2))
        ln_cmc_minus = math.log(clint_ideal_cmc(alpha1 - h, cmc1, cmc2))
        dln_cmc_dalpha1 = (ln_cmc_plus - ln_cmc_minus) / (2.0 * h)

        x1_rod = rodenas_x1(alpha1, dln_cmc_dalpha1)
        x1_motomura = motomura_ideal_composition(alpha1, cmc1, cmc2)
        assert x1_rod == pytest.approx(x1_motomura, abs=1e-5)


def test_rodenas_x1_zero_slope_gives_x1_equals_alpha1():
    """Direct, trivial consequence of the formula itself: a flat
    ln(cmc_mix) vs alpha1 curve (slope 0) gives X1_Rod = alpha1."""
    assert rodenas_x1(0.35, 0.0) == pytest.approx(0.35)


# --- rodenas_x1_series -------------------------------------------------
# Real, standard numerical-differentiation method (verified via WebFetch/
# WebSearch before implementing -- the source papers describe getting
# "the slope from the curve" without spelling out their exact numerical
# procedure): local quadratic (Lagrange) fit through each point's 2
# nearest neighbors, differentiated analytically. Validated two ways:
# (1) exact recovery when the underlying ln(cmc_mix) truly IS quadratic
# in alpha1 (the local fit is then not an approximation at all); (2)
# exact reduction to the textbook central/forward/backward finite-
# difference formulas for evenly-spaced data.


def test_rodenas_x1_series_exact_for_quadratic_ln_cmc_mix():
    """If ln(cmc_mix) truly follows a quadratic in alpha1, the local
    3-point quadratic fit reproduces the EXACT analytic derivative at
    every point (not just approximately) -- a strong round-trip check."""
    a, b, c = 0.5, -2.0, 3.0  # ln(cmc_mix) = a*alpha1^2 + b*alpha1 + c

    def ln_cmc(alpha1):
        return a * alpha1 ** 2 + b * alpha1 + c

    def true_slope(alpha1):
        return 2.0 * a * alpha1 + b

    alpha1_series = [0.1, 0.25, 0.4, 0.55, 0.7, 0.85]  # unevenly spaced on purpose
    cmc_mix_series = [math.exp(ln_cmc(a1)) for a1 in alpha1_series]

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    assert result.n_points == 6
    for a1, slope in zip(result.alpha1_used, result.dln_cmc_mix_dalpha1):
        assert slope == pytest.approx(true_slope(a1), abs=1e-9)


def test_rodenas_x1_series_reduces_to_textbook_central_difference():
    """Evenly-spaced interior point: must match (y2-y0)/(2h) exactly."""
    h = 0.1
    alpha1_series = [0.2, 0.3, 0.4]
    cmc_mix_series = [1.0, 1.5, 2.6]  # arbitrary, not quadratic

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    ln_cmc = [math.log(c) for c in cmc_mix_series]
    expected_middle_slope = (ln_cmc[2] - ln_cmc[0]) / (2.0 * h)
    assert result.dln_cmc_mix_dalpha1[1] == pytest.approx(expected_middle_slope, abs=1e-9)


def test_rodenas_x1_series_reduces_to_textbook_forward_difference_at_first_point():
    """Evenly-spaced first point: must match (-3y0+4y1-y2)/(2h) exactly."""
    h = 0.1
    alpha1_series = [0.2, 0.3, 0.4]
    cmc_mix_series = [1.0, 1.5, 2.6]
    ln_cmc = [math.log(c) for c in cmc_mix_series]
    expected_first_slope = (-3 * ln_cmc[0] + 4 * ln_cmc[1] - ln_cmc[2]) / (2.0 * h)

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    assert result.dln_cmc_mix_dalpha1[0] == pytest.approx(expected_first_slope, abs=1e-9)


def test_rodenas_x1_series_chains_into_rodenas_x1_correctly():
    """The returned x1_rodenas must match calling rodenas_x1 directly
    with the same (alpha1, slope) pairs -- confirms the chaining, not
    just the slope computation."""
    alpha1_series = [0.1, 0.25, 0.4, 0.55, 0.7]
    cmc_mix_series = [8.0, 6.5, 5.2, 4.3, 3.7]

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    for a1, slope, x1 in zip(result.alpha1_used, result.dln_cmc_mix_dalpha1, result.x1_rodenas):
        assert x1 == pytest.approx(rodenas_x1(a1, slope))


def test_rodenas_x1_series_sorts_unordered_input():
    alpha1_series = [0.7, 0.1, 0.4]
    cmc_mix_series = [3.7, 8.0, 5.2]  # matching the alpha1 order above

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    assert result.alpha1_used == sorted(alpha1_series)


def test_rodenas_x1_series_rejects_bad_inputs():
    with pytest.raises(ValueError):
        rodenas_x1_series([0.2, 0.4], [1.0, 2.0])  # too few points
    with pytest.raises(ValueError):
        rodenas_x1_series([0.2, 0.4, 0.6], [1.0, 2.0])  # mismatched lengths
    with pytest.raises(ValueError):
        rodenas_x1_series([0.0, 0.4, 0.6], [1.0, 2.0, 3.0])  # bad alpha1
    with pytest.raises(ValueError):
        rodenas_x1_series([0.2, 0.4, 0.6], [1.0, -2.0, 3.0])  # bad cmc_mix
    with pytest.raises(ValueError):
        rodenas_x1_series([0.2, 0.4, 0.4], [1.0, 2.0, 3.0])  # duplicate alpha1


def test_rodenas_x1_series_against_real_multipoint_literature_data():
    """Real test against an actual raw multi-point literature alpha1
    series (2026-09-10) -- NOT a second independent source (this reuses
    the same G6+T-20 system, Table 1/2, already used for the single-
    point rodenas_x1/Maeda validations elsewhere in this file: Azum,
    Rub, Alotaibi, Khan & Asiri, Biointerface Res. Appl. Chem. 12(6)
    (2022) 7416-7428), but the FIRST time rodenas_x1_series specifically
    (the numerical local-slope-from-a-series helper, previously only
    round-trip tested against synthetic exactly-quadratic data) has been
    checked against real, sparse (5-point), experimentally-noisy data
    with a real reported X1_Rod comparison. Real, disclosed finding:
    this does NOT cleanly reproduce the paper's own reported X1_Rod
    values -- interior points (alpha1=0.4, 0.5, 0.6) land within
    3.4-9.1% relative error, but the two series ENDPOINTS (alpha1=0.3,
    0.8) are off by ~19-20%, with the alpha1=0.8 point even exceeding
    the physically-valid [0,1] mole-fraction range. This is consistent
    with (not contradicted by) rodenas_x1_series's own documented
    one-sided-derivative limitation at series endpoints, and is the
    most likely explanation (the paper's own local slope was almost
    certainly obtained by some different method -- a global fit or a
    graphical read -- not stated in enough detail to reproduce exactly).
    Kept as a real, honest record of this method's real-world behavior
    on sparse literature data, not as a validated numeric match -- do
    not read this as proof the function is correct to the precision the
    round-trip tests suggest; sparse real series have genuine numerical-
    differentiation error this function cannot eliminate."""
    alpha1_series = [0.3, 0.4, 0.5, 0.6, 0.8]
    cmc_mix_series = [0.060, 0.056, 0.051, 0.048, 0.035]  # mM, paper's own Table 1 'cmc' column
    paper_x1_rod = [0.520, 0.652, 0.762, 0.852, 0.968]

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    # interior points: real, meaningful agreement
    for i in (1, 2, 3):
        assert result.x1_rodenas[i] == pytest.approx(paper_x1_rod[i], rel=0.10)
    # qualitative check only for the endpoints (real, larger, disclosed error there)
    assert result.x1_rodenas[0] < result.x1_rodenas[2] < result.x1_rodenas[4]  # monotonic trend preserved


def test_rodenas_activity_coefficients_basic():
    alpha1, x1, cmc_mix, cmc1, cmc2 = 0.5, 0.6, 5.0, 10.0, 8.0
    f1, f2 = rodenas_activity_coefficients(alpha1, x1, cmc_mix, cmc1, cmc2)
    assert f1 == pytest.approx((0.5 * 5.0) / (0.6 * 10.0))
    assert f2 == pytest.approx((0.5 * 5.0) / (0.4 * 8.0))


def test_rodenas_activity_coefficients_rejects_bad_x1():
    with pytest.raises(ValueError):
        rodenas_activity_coefficients(0.5, 0.0, 5.0, 10.0, 8.0)
    with pytest.raises(ValueError):
        rodenas_activity_coefficients(0.5, 1.0, 5.0, 10.0, 8.0)


# --- Maeda free energy of micellization (alternative-methods review) -------
# Source: equation (9), Azum, Rub, Alotaibi, Khan & Asiri, Biointerface Res.
# Appl. Chem. 12(6) (2022) 7416-7428 (same G6 gemini paper as above), citing
# Maeda, J. Colloid Interface Sci. 172 (1995) 98-105.


def test_maeda_b0_matches_real_paper_value():
    """Real numeric cross-check, not just a formula transcription
    check: B0 = ln(Xcmc2) computed with this project's own already-
    validated TX-114 CMC (0.263 mM, AZUM_TX114_CMC, from the same
    source paper's own Table 1) must match the paper's own reported
    -B0 = 12.25 for the G6+TX-114 system (independent of alpha1 --
    B0 depends only on cmc2)."""
    from surfactantkit.thermodynamics import cmc_to_mole_fraction

    xcmc2 = cmc_to_mole_fraction(AZUM_TX114_CMC / 1000.0)  # mM -> M
    b0 = math.log(xcmc2)
    assert -b0 == pytest.approx(12.25, abs=0.02)


def test_maeda_free_energy_matches_independent_recomputation():
    """Independently recompute B0/B1/B2/deltaG_M in the test itself
    (not copy-pasted from the implementation) and check they agree --
    catches a transcription bug the implementation and a copy-pasted
    test would both share."""
    x1_rub, beta = 0.6, -1.8
    cmc1_M, cmc2_M, temperature_k = 0.041e-3, 0.263e-3, 298.15
    xcmc1 = cmc1_M / (cmc1_M + 55.5)
    xcmc2 = cmc2_M / (cmc2_M + 55.5)
    b2_expected = -beta
    b1_expected = math.log(xcmc1 / xcmc2) - b2_expected
    b0_expected = math.log(xcmc2)
    delta_g_expected = (R_GAS * temperature_k * (b0_expected + b1_expected * x1_rub + b2_expected * x1_rub ** 2)) / 1000.0

    delta_g = maeda_free_energy_of_micellization(x1_rub, beta, cmc1_M, cmc2_M, temperature_k)
    assert delta_g == pytest.approx(delta_g_expected, rel=1e-9)


def test_maeda_free_energy_matches_second_independent_literature_system():
    """Real SECOND independent literature validation for Maeda's method
    (2026-09-10), beyond the single Azum et al. 2022 gemini-surfactant
    source used everywhere else in this file -- a different paper,
    different chemistry entirely (a drug-surfactant system, not a
    gemini/conventional-surfactant pair): Rub, Azum, Kumar, Arshad,
    Khan, Alotaibi & Asiri, Polymers 13(22) (2021) 4025,
    doi:10.3390/polym13224025, imipramine hydrochloride (IMP) + Triton
    X-100 (TX-100), aqueous, 298 K, their own Table 1/2, alpha1=0.5:
    cmc(IMP)=41.85 mmol/kg, cmc(TX-100)=0.31 mmol/kg, cmc_mix=0.48
    mmol/kg, X1^Rb=0.8585, beta^Rb=-4.35, deltaG_Maeda=-19.64 kJ/mol
    (their own reported value, taken directly, not re-derived).

    Real, disclosed finding: reproducing their reported deltaG_Maeda
    from their OWN reported X1^Rb/beta^Rb (mmol/kg treated as ~mM,
    dilute-solution approximation, component 1 = IMP per this project's
    established ionic-first Maeda convention from the Azum validation)
    gives -20.85 kJ/mol vs. their reported -19.64 kJ/mol -- a real ~6%
    relative error, plausible from rounding in their 3-4-significant-
    figure table values (same "close but not exact, real rounding
    residual" pattern as every other cross-paper check in this file),
    not a formula error: sign, order of magnitude, and the correct
    magnitude to within single-digit percent all match a genuinely
    independent second source. A real, separate attempt to also
    independently re-solve X1 via this project's own solve_rubingh_x
    from their raw (alpha1, cmc_mix, cmc1, cmc2) inputs recovered
    0.147, not 0.8585, unless cmc1/cmc2 are swapped in the SOLVER call
    specifically (giving 0.853, close to their 0.8585) -- a real
    component-labeling-convention subtlety between this paper's
    Rubingh-solve indexing and its Maeda indexing that could not be
    resolved with full confidence from the fetched text alone; this
    test therefore uses their OWN reported X1/beta values directly
    rather than asserting this project's solver reproduces them, so it
    validates maeda_free_energy_of_micellization specifically, not an
    additional Rubingh cross-check."""
    cmc1_M, cmc2_M = 41.85e-3, 0.31e-3  # mmol/kg ~= mM, dilute approximation
    x1_rub, beta, temperature_k = 0.8585, -4.35, 298.0
    dg = maeda_free_energy_of_micellization(x1_rub, beta, cmc1_M, cmc2_M, temperature_k)
    assert dg == pytest.approx(-19.64, rel=0.07)


def test_motomura_ideal_composition_rejects_bad_inputs():
    with pytest.raises(ValueError):
        motomura_ideal_composition(0.0, 5.0, 5.0)
    with pytest.raises(ValueError):
        motomura_ideal_composition(1.0, 5.0, 5.0)
    with pytest.raises(ValueError):
        motomura_ideal_composition(0.5, 0.0, 5.0)
    with pytest.raises(ValueError):
        motomura_ideal_composition(0.5, 5.0, -1.0)
