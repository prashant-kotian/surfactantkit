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
    activity_coefficients,
    clint_ideal_cmc,
    excess_free_energy,
    rubingh_beta,
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
