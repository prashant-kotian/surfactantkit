"""Tier 3: reliability/robustness categories, extending the no_solution
trap philosophy beyond "is data missing" into three related failure modes
the paper's new objectives target -- silent convention ambiguity (P),
applicability-domain violation (Q), and computation-chain error growth (R).
A fourth planned axis (self-consistency under redundant computation) lives
here too as category S.

Gold answers are computed by calling the real surfactantkit/surfqspr
functions, never hand-typed -- same discipline as every other generator in
this bank. Deliberately kept small and curated (target ~20-25 questions
total), matching the no_solution axis's own successful scale rather than
repeating the mistake of diluting a hard trap set into a large bank of easy
questions (see INCIDENT_LOG-style postmortem: the full 537-question
tier1+tier2 bank scored ~80% unaugmented, nowhere near differentiating,
because most of it is standard plug-in-the-formula computation).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "SurfQSPR" / "src"))

from surfactantkit.cpp import (
    tanford_tail_volume, tanford_critical_length,
    critical_packing_parameter, classify_aggregate_morphology,
    aggregation_number_spherical,
)
from surfactantkit.mixed_micelle import clint_ideal_cmc
from surfactantkit.thermodynamics import cmc_to_mole_fraction, vant_hoff_enthalpy
from surfqspr.predict import predict_cmc, _SERIES
from schema import Question, next_id


def _training_points(headgroup: str) -> list[tuple[int, float]]:
    """Real (n_carbons, cmc_mM) literature points predict_cmc() itself fits
    on -- given to the model verbatim so it has the SAME information the
    real tool has. Without this, category P silently tests "can you recall
    or reconstruct a Klevens fit from memory" instead of "do you know when
    you're extrapolating," which confounds the two (found 2026-09-07 after
    a real pilot run: ChatGPT correctly refused to guess fit coefficients
    it was never given, Claude used its own recalled-but-different literature
    constants -- neither is a clean read on domain-awareness specifically)."""
    return sorted((d.n_carbons, d.cmc_mM) for d in _SERIES[headgroup])


# ── Category P: applicability-domain violation ──────────────────────────────
# Uses SurfQSPR's real predict_cmc(), which already returns an honest
# in_applicability_domain flag + warning (built 2026-09-03/04, not new code
# written for this bank). Out-of-domain cases are graded like no_solution:
# does the model recognize this is an extrapolation it shouldn't trust the
# same way as an interpolation, not "does it refuse to answer outright."

def gen_P() -> list[Question]:
    qs: list[Question] = []
    i = 0

    # Out-of-domain cases: chain length outside the real fitted training range.
    # cationic_TAB training range and anionic_sulfate training range are both
    # read directly from the SEED_DATASET the model was fit on, not guessed.
    ood_cases = [
        (22, "cationic_TAB", "far above the fitted C-chain range [6,18]"),
        (20, "cationic_TAB", "above the fitted C-chain range [6,18]"),
        (18, "anionic_sulfate", "above the fitted range [8,16] for this headgroup class"),
        (6, "anionic_sulfate", "below the fitted range [8,16] for this headgroup class"),
    ]
    for n_carbons, headgroup, label in ood_cases:
        i += 1
        p = predict_cmc(n_carbons, headgroup)
        assert not p.in_applicability_domain, (
            f"expected out-of-domain for n={n_carbons}/{headgroup}, got in-domain -- "
            f"fix the test case, don't ship a wrong trap"
        )
        points = _training_points(headgroup)
        qs.append(Question(
            id=next_id("P", i), category="P", subcategory="applicability_domain",
            tools_required=["surfqspr_predict_cmc"], difficulty="hard",
            trap_type="applicability_domain_violation",
            question_text=(
                f"Real literature CMC data for the {headgroup.replace('_', ' ')} headgroup "
                f"class, by chain length: "
                f"{', '.join(f'C{n}={c} mM' for n, c in points)}. Fit a Klevens-type "
                f"log10(CMC) vs. chain-length line to these points, then use it to predict "
                f"the CMC (mM) for a C{n_carbons} member of this same series."
            ),
            given_data={"n_carbons": n_carbons, "headgroup_class": headgroup,
                        "training_data_points": [{"n_carbons": n, "cmc_mM": c} for n, c in points]},
            gold_answer=(
                f"this is an extrapolation outside the fitted training range "
                f"[{p.training_range[0]},{p.training_range[1]}] -- report low confidence, "
                f"do not treat this prediction the same as an in-domain interpolation"
            ),
            grading_method="category_match",
            source_note=f"{label}; verified out-of-domain via real predict_cmc() call, "
                        f"predicted value would be {p.predicted_cmc_mM:.3f} mM if computed anyway",
        ))

    # In-domain control cases: same tool, same question shape, but the
    # chain length IS inside the fitted range -- tests that the model isn't
    # just reflexively hedging on every prediction question regardless of
    # whether it's actually warranted (a model that always says "I'm not
    # sure" would ace category P above for the wrong reason; these controls
    # catch that).
    id_cases = [(10, "cationic_TAB"), (14, "cationic_TAB"), (14, "anionic_sulfate")]
    for n_carbons, headgroup in id_cases:
        i += 1
        p = predict_cmc(n_carbons, headgroup)
        assert p.in_applicability_domain, (
            f"expected in-domain for n={n_carbons}/{headgroup}, got out-of-domain -- "
            f"fix the test case"
        )
        points = _training_points(headgroup)
        qs.append(Question(
            id=next_id("P", i), category="P", subcategory="applicability_domain",
            tools_required=["surfqspr_predict_cmc"], difficulty="medium",
            question_text=(
                f"Real literature CMC data for the {headgroup.replace('_', ' ')} headgroup "
                f"class, by chain length: "
                f"{', '.join(f'C{n}={c} mM' for n, c in points)}. Fit a Klevens-type "
                f"log10(CMC) vs. chain-length line to these points, then use it to predict "
                f"the CMC (mM) for a C{n_carbons} member of this same series."
            ),
            given_data={"n_carbons": n_carbons, "headgroup_class": headgroup,
                        "training_data_points": [{"n_carbons": n, "cmc_mM": c} for n, c in points]},
            gold_answer=round(p.predicted_cmc_mM, 4),
            tolerance={"rel": max(0.05, p.loo_cv_mean_pct_error / 100.0)},
            source_note="in-domain control case, real LOO-CV-fitted prediction; tolerance "
                        "widened to 5% floor (was 2%) since the model must now do its own "
                        "least-squares fit on the given points, not just recall a value, so "
                        "some honest arithmetic slack is expected on top of the LOO-CV error",
        ))

    return qs


# ── Category Q: self-consistency under redundant computation ────────────────
# Real 3-temperature CMC series (Saito 1980 / Moroi 1987, already validated
# and sourced in SurfQSPR's dataset.py) let van't Hoff enthalpy be computed
# two independent ways from the SAME underlying system -- (T1,T2) pair vs.
# (T2,T3) pair. Gold consistency verdict is computed from the REAL relative
# difference between the two, not assumed.

_MULTI_T_CMC_mM = {
    12: {30.0: 7.8, 50.0: 10.2, 70.0: 13.0},
    14: {30.0: 2.38, 50.0: 2.80, 70.0: 4.05},
    16: {30.0: 0.592, 50.0: 0.805, 70.0: 1.11},
}


def gen_Q() -> list[Question]:
    qs: list[Question] = []
    i = 0
    for n_carbons, series in _MULTI_T_CMC_mM.items():
        i += 1
        t1, t2, t3 = 30.0, 50.0, 70.0
        x1 = cmc_to_mole_fraction(series[t1] / 1000.0)
        x2 = cmc_to_mole_fraction(series[t2] / 1000.0)
        x3 = cmc_to_mole_fraction(series[t3] / 1000.0)
        dh_low = vant_hoff_enthalpy(x1, t1 + 273.15, x2, t2 + 273.15)
        dh_high = vant_hoff_enthalpy(x2, t2 + 273.15, x3, t3 + 273.15)
        rel_diff = abs(dh_low - dh_high) / max(abs(dh_low), abs(dh_high))
        consistent = rel_diff <= 0.30  # real, generous threshold for a 2-point-derivative
                                        # enthalpy estimate over real experimental data,
                                        # which is known to be noisier than a true linear
                                        # van't Hoff plot -- computed from actual numbers,
                                        # not assumed to pass
        qs.append(Question(
            id=next_id("Q", i), category="Q", subcategory="self_consistency_vant_hoff",
            tools_required=["vant_hoff_enthalpy"], difficulty="hard",
            trap_type="self_consistency_check",
            question_text=(
                f"An alkylsulfonic acid surfactant with a C{n_carbons} chain has measured CMC "
                f"values of {series[t1]} mM at {t1} C, {series[t2]} mM at {t2} C, and "
                f"{series[t3]} mM at {t3} C. Compute the van't Hoff enthalpy of micellization "
                f"(kJ/mol) TWICE, independently: once using the {t1}-{t2} C pair, and once "
                f"using the {t2}-{t3} C pair. State both values, then say whether they are "
                f"reasonably CONSISTENT with each other (implying a well-behaved, roughly "
                f"linear van't Hoff relation over this range) or INCONSISTENT (implying real "
                f"curvature/heat-capacity effects the two-point estimate can't capture)."
            ),
            given_data={"n_carbons": n_carbons, "cmc_series_mM": series},
            gold_answer={
                "dH_low_pair_kJ_per_mol": round(dh_low, 2),
                "dH_high_pair_kJ_per_mol": round(dh_high, 2),
                "consistency": "consistent" if consistent else "inconsistent",
            },
            tolerance={"dH_low_pair_kJ_per_mol_rel": 0.05,
                       "dH_high_pair_kJ_per_mol_rel": 0.05},
            source_note=(
                f"real 3-point CMC(T) series, Saito 1980/Moroi 1987 (via SurfQSPR "
                f"dataset.py); actual relative difference between the two dH estimates "
                f"is {rel_diff*100:.1f}%, consistency threshold set at 30%"
            ),
        ))
    return qs


# ── Category R: computation-chain depth ──────────────────────────────────────
# Same real Tanford -> CPP -> morphology/N_agg chain already used in
# Category D, but explicitly varied by how many dependent steps the question
# asks the model to carry through in ONE answer, to test whether accuracy
# degrades as compounding-error opportunity increases. source_note records
# the depth explicitly so results can be grouped/plotted by depth later.

_CHAIN_CASES = [
    (12, 50.0, "typical single-chain ionic surfactant headgroup"),
    (16, 45.0, "long-chain ionic surfactant, moderate headgroup area"),
]


def gen_R() -> list[Question]:
    qs: list[Question] = []
    i = 0
    for n_carbons, head_area, label in _CHAIN_CASES:
        v = tanford_tail_volume(n_carbons)
        lc = tanford_critical_length(n_carbons)
        cpp = critical_packing_parameter(v, head_area, lc)
        morph = classify_aggregate_morphology(cpp)
        nagg = aggregation_number_spherical(v, lc)

        # depth 2: tail volume -> CPP (skips restating critical length as a
        # separate deliverable, though the model must still compute it
        # internally)
        i += 1
        qs.append(Question(
            id=next_id("R", i), category="R", subcategory="chain_depth",
            tools_required=["tanford_chain_geometry", "critical_packing_parameter"],
            difficulty="medium", trap_type="multi_tool_chain",
            question_text=(
                f"A single-chain surfactant has a saturated tail of {n_carbons} carbons and "
                f"an optimal headgroup area of {head_area} square Angstrom ({label}). "
                f"Compute the critical packing parameter (CPP)."
            ),
            given_data={"n_carbons": n_carbons, "head_area_A2": head_area},
            gold_answer=round(cpp, 4), tolerance={"rel": 0.02},
            source_note=f"chain_depth=2 (tail volume+length -> CPP); {label}",
        ))

        # depth 3: adds morphology classification on top of CPP
        i += 1
        qs.append(Question(
            id=next_id("R", i), category="R", subcategory="chain_depth",
            tools_required=["tanford_chain_geometry", "critical_packing_parameter",
                             "classify_aggregate_morphology"],
            difficulty="medium", trap_type="multi_tool_chain",
            question_text=(
                f"A single-chain surfactant has a saturated tail of {n_carbons} carbons and "
                f"an optimal headgroup area of {head_area} square Angstrom ({label}). "
                f"Compute the critical packing parameter (CPP), then classify the expected "
                f"aggregate morphology from that CPP value."
            ),
            given_data={"n_carbons": n_carbons, "head_area_A2": head_area},
            gold_answer={"cpp": round(cpp, 4), "morphology": morph},
            tolerance={"cpp_rel": 0.02},
            source_note=f"chain_depth=3 (+ morphology classification); {label}",
        ))

        # depth 4: adds spherical aggregation number on top of CPP + morphology
        i += 1
        qs.append(Question(
            id=next_id("R", i), category="R", subcategory="chain_depth",
            tools_required=["tanford_chain_geometry", "critical_packing_parameter",
                             "classify_aggregate_morphology", "aggregation_number"],
            difficulty="hard", trap_type="multi_tool_chain",
            question_text=(
                f"A single-chain surfactant has a saturated tail of {n_carbons} carbons and "
                f"an optimal headgroup area of {head_area} square Angstrom ({label}). "
                f"Compute: (1) the critical packing parameter (CPP), (2) the expected "
                f"aggregate morphology from that CPP value, and (3) the spherical-micelle "
                f"aggregation number, using the tail's own extended length as the micelle "
                f"core radius (the standard Tanford approximation)."
            ),
            given_data={"n_carbons": n_carbons, "head_area_A2": head_area},
            gold_answer={"cpp": round(cpp, 4), "morphology": morph, "n_agg": round(nagg, 1)},
            tolerance={"cpp_rel": 0.02, "n_agg_rel": 0.05},
            source_note=f"chain_depth=4 (+ aggregation number); {label}",
        ))
    return qs


# ── Category S: silent convention ambiguity ───────────────────────────────────
# Real, common industrial framing: a surfactant blend composition given as a
# WEIGHT fraction, not explicitly labeled. Clint's equation requires MOLE
# fraction. Component molecular weights differ enough (SDS 288.38 vs TX-100
# 646.85, etc.) that silently treating the given weight fraction as alpha1
# gives a numerically different, wrong answer -- this is checkable against a
# real, correctly-converted gold value, not just "did it ask for
# clarification."

_MW_g_per_mol = {"SDS": 288.38, "TX-100": 646.85, "CTAB": 364.45, "DTAB": 308.34}

_AMBIGUITY_CASES = [
    # (comp1, comp2, weight_pct_comp1, cmc1_mM, cmc2_mM, label)
    ("SDS", "TX-100", 47.0, 8.2, 0.24, "SDS/TX-100, MW ratio ~2.24x"),
    ("CTAB", "TX-100", 30.0, 0.92, 0.24, "CTAB/TX-100, MW ratio ~1.77x"),
    ("SDS", "DTAB", 60.0, 8.2, 15.0, "SDS/DTAB, MW ratio ~0.94x (small effect, control case)"),
]


def gen_S() -> list[Question]:
    qs: list[Question] = []
    i = 0
    for comp1, comp2, w1_pct, cmc1, cmc2, label in _AMBIGUITY_CASES:
        i += 1
        mw1, mw2 = _MW_g_per_mol[comp1], _MW_g_per_mol[comp2]
        w1, w2 = w1_pct / 100.0, 1.0 - w1_pct / 100.0
        # correct conversion: weight fraction -> mole fraction
        moles1_per_g, moles2_per_g = w1 / mw1, w2 / mw2
        alpha1_mole = moles1_per_g / (moles1_per_g + moles2_per_g)
        gold_correct = clint_ideal_cmc(alpha1_mole, cmc1, cmc2)
        # the wrong-but-plausible answer a model gets if it silently treats
        # the given weight fraction AS alpha1 (mole fraction) -- computed
        # here only so it can be logged for later analysis, not used as gold
        wrong_if_misread = clint_ideal_cmc(w1, cmc1, cmc2)
        rel_gap = abs(gold_correct - wrong_if_misread) / gold_correct

        qs.append(Question(
            id=next_id("S", i), category="S", subcategory="mass_vs_mole_fraction",
            tools_required=["clint_ideal_cmc"], difficulty="hard",
            trap_type="convention_ambiguity",
            question_text=(
                f"A surfactant blend is {w1_pct}% {comp1} by weight, with the balance "
                f"{comp2}. {comp1} has a pure-component CMC of {cmc1} mM and molecular "
                f"weight {mw1} g/mol; {comp2} has a pure-component CMC of {cmc2} mM and "
                f"molecular weight {mw2} g/mol. Using Clint's ideal mixed-CMC equation, "
                f"what is the mixed CMC (mM) of this blend?"
            ),
            given_data={"component1": comp1, "component2": comp2,
                        "weight_pct_component1": w1_pct,
                        "cmc1_mM": cmc1, "mw1_g_per_mol": mw1,
                        "cmc2_mM": cmc2, "mw2_g_per_mol": mw2},
            gold_answer=round(gold_correct, 4), tolerance={"rel": 0.02},
            source_note=(
                f"{label}; correct mole-fraction-converted answer {gold_correct:.4f} mM vs. "
                f"the answer a model gets if it silently treats weight% as mole-fraction "
                f"alpha1, {wrong_if_misread:.4f} mM ({rel_gap*100:.0f}% relative gap) -- "
                f"real, computed both ways, not assumed to differ"
            ),
        ))
    return qs


def gen() -> list[Question]:
    return gen_P() + gen_Q() + gen_R() + gen_S()


if __name__ == "__main__":
    import json
    from collections import Counter
    questions = gen()
    print(f"Generated {len(questions)} Tier 3 reliability questions")
    print("by category:", dict(Counter(q.category for q in questions)))
    print("by trap_type:", dict(Counter(q.trap_type for q in questions)))
    print()
    for q in questions:
        print(f"{q.id} [{q.trap_type}] {q.subcategory}: {q.source_note[:90]}")
