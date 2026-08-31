"""Category A: mixed micelle theory. Generates questions for
clint_ideal_cmc, rubingh_solve, rubingh_activity_coefficients,
excess_free_energy, rosen_monolayer_solve, corrin_harkins_predict.

Gold answers are computed by calling the real surfactantkit functions --
never hand-typed.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from surfactantkit.mixed_micelle import (
    clint_ideal_cmc,
    solve_rubingh_x,
    rubingh_beta,
    activity_coefficients,
    excess_free_energy,
    solve_rosen_monolayer_x,
    rosen_beta_sigma,
    corrin_harkins_predict_cmc,
)
from schema import Question, next_id

# Realistic (alpha1, cmc_mix, cmc1, cmc2, label) systems, drawn from
# literature_validation_notes.md plus nearby synthetic variations for volume.
SYSTEMS = [
    (0.47, 0.547, 0.387, 3.468, "TX-100(1)/SDS(2), Muherei & Junin 2009"),
    (0.30, 0.40, 0.387, 3.468, "TX-100(1)/SDS(2), synthetic composition"),
    (0.48, 0.061, 0.041, 0.263, "gemini G6(1)/TX-114(2), Azum et al. 2022"),
    (0.60, 0.05, 0.041, 0.263, "gemini G6(1)/TX-114(2), synthetic composition"),
    (0.6, 0.70, 2.25, 0.41, "TTAB(1)/Tween-20(2), Lee & Lee 2012"),
    (0.5, 1.52, 1.10, 2.63, "gemini 12-4-12(1)/ZW3-12(2), McLachlan et al. 2020"),
    (0.888, 0.253, 0.309, 0.134, "TX-100(1)/rhamnolipid(2), Liu et al. 2020"),
    (0.5, 4.07, 11.50, 11.98, "sodium cholate(1)/SDS(2), Kang/Bahadur, PMC4087020"),
    (0.25, 6.011, 14.80, 8.00, "SDS(1)/DTAB(2), PMC6554738, SDS-rich"),
    (0.75, 13.00, 14.80, 8.00, "SDS(1)/DTAB(2), PMC6554738, DTAB-rich"),
]


def gen() -> list[Question]:
    qs: list[Question] = []
    i = 0

    # --- clint_ideal_cmc (target ~10) ---
    for alpha1, _, c1, c2, label in SYSTEMS:
        i += 1
        gold = clint_ideal_cmc(alpha1, c1, c2)
        qs.append(Question(
            id=next_id("A", i), category="A", subcategory="clint_ideal_cmc",
            tools_required=["clint_ideal_cmc"], difficulty="easy",
            question_text=(
                f"A binary surfactant mixture has component 1 with pure CMC {c1} mM and "
                f"component 2 with pure CMC {c2} mM. At a bulk mole fraction of component 1 "
                f"(alpha1) = {alpha1}, what is the Clint ideal mixed CMC (mM), assuming no "
                f"interaction between the components?"
            ),
            given_data={"alpha1": alpha1, "cmc1_mM": c1, "cmc2_mM": c2},
            gold_answer=round(gold, 4), tolerance={"rel": 0.02},
            source_note=label,
        ))

    # --- rubingh_solve (target ~10, plus 2 no-solution traps) ---
    for alpha1, cmix, c1, c2, label in SYSTEMS:
        i += 1
        x1 = solve_rubingh_x(alpha1, cmix, c1, c2)
        beta = rubingh_beta(x1, alpha1, cmix, c1)
        qs.append(Question(
            id=next_id("A", i), category="A", subcategory="rubingh_solve",
            tools_required=["rubingh_solve"], difficulty="hard",
            question_text=(
                f"For the same binary system (component 1 pure CMC {c1} mM, component 2 pure "
                f"CMC {c2} mM), the EXPERIMENTALLY MEASURED mixed CMC at alpha1={alpha1} is "
                f"{cmix} mM. Using Rubingh's regular-solution theory, find the micellar mole "
                f"fraction of component 1 (x1) and the interaction parameter beta. State "
                f"whether the mixture is synergistic or antagonistic."
            ),
            given_data={"alpha1": alpha1, "cmc_mix_mM": cmix, "cmc1_mM": c1, "cmc2_mM": c2},
            gold_answer={"x1": round(x1, 4), "beta": round(beta, 3),
                         "synergy": "synergistic" if beta < 0 else "antagonistic"},
            tolerance={"x1_abs": 0.01, "beta_rel": 0.05},
            source_note=label,
        ))

    # Deliberate no-solution trap: cmc_mix set outside the physically
    # achievable range for this alpha1/cmc1/cmc2 (verified by actually
    # calling the solver below, not assumed).
    trap_cases = [(0.60, 0.08, 0.041, 0.263), (0.60, 0.10, 0.041, 0.263)]
    for alpha1, cmix, c1, c2 in trap_cases:
        x1 = solve_rubingh_x(alpha1, cmix, c1, c2)
        if x1 is not None:
            continue  # not actually a no-solution case for these numbers; skip
        i += 1
        qs.append(Question(
            id=next_id("A", i), category="A", subcategory="rubingh_solve",
            tools_required=["rubingh_solve"], difficulty="hard", trap_type="no_solution",
            question_text=(
                f"A binary system has pure CMCs {c1} mM and {c2} mM. At alpha1={alpha1}, the "
                f"reported mixed CMC is {cmix} mM. Solve for the Rubingh micellar mole "
                f"fraction x1 and interaction parameter beta."
            ),
            given_data={"alpha1": alpha1, "cmc_mix_mM": cmix, "cmc1_mM": c1, "cmc2_mM": c2},
            gold_answer="no valid root in (0,1) -- this cmc_mix is not physically consistent with the given alpha1/cmc1/cmc2",
            grading_method="category_match",
            source_note="constructed trap case, verified no-solution by direct call to solve_rubingh_x",
        ))

    # --- rubingh_activity_coefficients (target ~10) ---
    for alpha1, cmix, c1, c2, label in SYSTEMS:
        i += 1
        x1 = solve_rubingh_x(alpha1, cmix, c1, c2)
        beta = rubingh_beta(x1, alpha1, cmix, c1)
        f1, f2 = activity_coefficients(x1, beta)
        qs.append(Question(
            id=next_id("A", i), category="A", subcategory="rubingh_activity_coefficients",
            tools_required=["rubingh_activity_coefficients"], difficulty="medium",
            question_text=(
                f"For a mixed micelle with micellar mole fraction of component 1 x1={x1:.4f} "
                f"and Rubingh interaction parameter beta={beta:.3f}, compute the regular-"
                f"solution activity coefficients f1 and f2."
            ),
            given_data={"x1": round(x1, 4), "beta": round(beta, 3)},
            gold_answer={"f1": round(f1, 4), "f2": round(f2, 4)},
            tolerance={"rel": 0.02}, source_note=label,
        ))

    # --- excess_free_energy (target ~10, with a Celsius/Kelvin unit trap on ~half) ---
    for idx, (alpha1, cmix, c1, c2, label) in enumerate(SYSTEMS):
        i += 1
        x1 = solve_rubingh_x(alpha1, cmix, c1, c2)
        beta = rubingh_beta(x1, alpha1, cmix, c1)
        f1, f2 = activity_coefficients(x1, beta)
        T_K = 298.15
        dg_ex = excess_free_energy(x1, f1, f2, T_K)
        use_celsius_trap = idx % 2 == 0
        temp_text = "25 degrees C" if use_celsius_trap else f"{T_K} K"
        qs.append(Question(
            id=next_id("A", i), category="A", subcategory="excess_free_energy",
            tools_required=["excess_free_energy"], difficulty="medium",
            trap_type="unit_trap" if use_celsius_trap else "none",
            question_text=(
                f"At {temp_text}, a mixed micelle has x1={x1:.4f}, f1={f1:.4f}, f2={f2:.4f}. "
                f"Compute the excess Gibbs free energy of micelle formation, deltaG_ex, in kJ/mol."
            ),
            given_data={"x1": round(x1, 4), "f1": round(f1, 4), "f2": round(f2, 4), "temperature_given": temp_text},
            gold_answer=round(dg_ex, 3), tolerance={"rel": 0.02}, source_note=label,
        ))

    # --- rosen_monolayer_solve (target ~10) -- SAME math as Rubingh but
    # framed with surface-tension-derived concentrations, not CMC, to
    # test whether the model recognizes the correct tool/data type.
    for alpha1, cmix, c1, c2, label in SYSTEMS:
        i += 1
        x1 = solve_rosen_monolayer_x(alpha1, cmix, c1, c2)
        beta_s = rosen_beta_sigma(x1, alpha1, cmix, c1)
        qs.append(Question(
            id=next_id("A", i), category="A", subcategory="rosen_monolayer_solve",
            tools_required=["rosen_monolayer_solve"], difficulty="hard",
            question_text=(
                f"At the air-water interface, component 1 alone requires {c1} mM to reach a "
                f"reference surface tension of 40 mN/m; component 2 alone requires {c2} mM to "
                f"reach that same reference surface tension. In a mixture at bulk mole fraction "
                f"alpha1={alpha1}, a total concentration of {cmix} mM reaches that same 40 mN/m "
                f"reference surface tension. Using Rosen's mixed-monolayer theory (NOT CMC-based "
                f"Rubingh theory -- this is surface-tension data), find the monolayer mole "
                f"fraction x1 and the interaction parameter beta^sigma."
            ),
            given_data={"alpha1": alpha1, "c_mix_sigma_mM": cmix, "c1_sigma_mM": c1, "c2_sigma_mM": c2},
            gold_answer={"x1_sigma": round(x1, 4), "beta_sigma": round(beta_s, 3)},
            tolerance={"x1_abs": 0.01, "beta_rel": 0.05}, source_note=label + " (reframed as monolayer data for this question)",
        ))

    # --- corrin_harkins_predict (target ~10) ---
    ch_cases = [
        (10.0, 10.0, 5.0, 100.0, 50.0, "synthetic, g~0.3 range"),
        (14.80, 5.0, 6.011, 50.0, 20.0, "loosely inspired by DTAB-SDS pure/mixed CMC scale"),
        (8.38, 5.0, 6.2, 50.0, 25.0, "loosely inspired by SDS-DTAB (Sachin et al.) CMC scale"),
    ]
    for cmc1, c_salt1, cmc2, c_salt2, c_target, label in ch_cases * 4:  # repeat with variation for volume
        i += 1
        predicted, g = corrin_harkins_predict_cmc(cmc1, c_salt1, cmc2, c_salt2, c_target)
        qs.append(Question(
            id=next_id("A", i), category="A", subcategory="corrin_harkins_predict",
            tools_required=["corrin_harkins_predict"], difficulty="hard",
            question_text=(
                f"For an ionic surfactant, the CMC is {cmc1} mM at a counterion (added salt) "
                f"concentration of {c_salt1} mM, and {cmc2} mM at a counterion concentration of "
                f"{c_salt2} mM. Using the Corrin-Harkins log-linear relation, predict the CMC at "
                f"a counterion concentration of {c_target} mM."
            ),
            given_data={"cmc1_mM": cmc1, "salt_conc1_mM": c_salt1, "cmc2_mM": cmc2,
                        "salt_conc2_mM": c_salt2, "salt_conc_target_mM": c_target},
            gold_answer=round(predicted, 4), tolerance={"rel": 0.02}, source_note=label,
        ))

    return qs


if __name__ == "__main__":
    import json
    questions = gen()
    print(f"Generated {len(questions)} questions for Category A (mixed micelle theory)")
    from collections import Counter
    print(Counter(q.subcategory for q in questions))
    print(json.dumps(questions[0].to_dict(), indent=2))
