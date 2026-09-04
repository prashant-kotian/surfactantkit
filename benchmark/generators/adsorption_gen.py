"""Category B: adsorption. gibbs_surface_excess, gibbs_area_per_molecule,
szyszkowski_predict_surface_tension."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from surfactantkit.adsorption import gibbs_gamma_max, gibbs_a_min, szyszkowski_surface_tension
from schema import Question, next_id

# (slope_mN_per_lnC, system_type, temperature_K, label)
# system_type replaced a raw n_factor float 2026-09-04 as part of hardening
# gibbs_gamma_max against silent guessing (see thermodynamics/adsorption
# no_solution trap results) -- these are all cases where the label ALREADY
# states the ionic character and electrolyte condition, so system_type is
# real given information here, not a guess.
GIBBS_CASES = [
    (-5.0, "ionic_no_added_salt", 298.15, "ionic surfactant, no added salt"),
    (-8.5, "ionic_no_added_salt", 298.15, "ionic surfactant, steeper pre-CMC slope, no added salt"),
    (-3.2, "nonionic", 298.15, "nonionic surfactant"),
    (-6.7, "nonionic", 308.15, "nonionic surfactant, elevated temperature"),
    (-4.4, "ionic_no_added_salt", 293.15, "ionic surfactant, no added salt, lower temperature"),
    (-9.1, "ionic_no_added_salt", 298.15, "gemini cationic surfactant, no added salt, strong adsorption"),
]

# Deliberate trap: an ionic surfactant's pre-CMC slope given WITHOUT stating
# whether excess inert electrolyte is present. gibbs_gamma_max's system_type
# (nonionic / ionic_excess_electrolyte / ionic_no_added_salt) is a required,
# undefaulted argument specifically because this choice is real,
# surfactant/condition-specific literature knowledge, not something to
# assume -- see the function's and MCP tool's own docstrings ("do not
# guess", "state which applies, do not assume"). Correct behavior: refuse
# to pick a value and flag the missing information, not silently default
# to nonionic or either ionic variant.
NO_SOLUTION_GIBBS_CASES = [
    (-5.5, 298.15, "an anionic surfactant (SDS-type), electrolyte condition not stated"),
    (-7.2, 298.15, "a cationic surfactant (CTAB-type), electrolyte condition not stated"),
    (-9.0, 298.15, "a gemini cationic surfactant, electrolyte condition not stated"),
    (-4.0, 303.15, "a 1:1 anionic surfactant in aqueous solution, electrolyte condition not stated"),
    (-6.1, 298.15, "an alkyl sulfate surfactant, electrolyte condition not stated"),
    (-8.0, 298.15, "a quaternary ammonium surfactant, electrolyte condition not stated"),
]

# (concentration, gamma0, gamma_max, K, system_type, temperature_K, label)
SZYSZKOWSKI_CASES = [
    (5.0, 72.0, 3.0e-6, 50.0, "nonionic", 298.15, "nonionic, moderate K"),
    (2.0, 72.0, 3.0e-6, 50.0, "nonionic", 298.15, "nonionic, lower concentration"),
    (8.0, 68.0, 2.5e-6, 30.0, "ionic_no_added_salt", 298.15, "ionic surfactant, no added salt, different gamma0"),
    (1.0, 72.0, 4.0e-6, 100.0, "nonionic", 298.15, "nonionic, high-affinity adsorption (large K)"),
]


def gen() -> list[Question]:
    qs: list[Question] = []
    i = 0

    for slope, system_type, T, label in GIBBS_CASES * 2:
        i += 1
        gamma_max = gibbs_gamma_max(slope, system_type, T)
        qs.append(Question(
            id=next_id("B", i), category="B", subcategory="gibbs_surface_excess",
            tools_required=["gibbs_surface_excess"], difficulty="medium",
            question_text=(
                f"Below the CMC, the slope of surface tension vs ln(concentration), "
                f"d(gamma)/d(ln C), is {slope} mN/m. This is {label}, at {T} K. "
                f"Compute the maximum surface excess concentration, Gamma_max, in mol/m^2."
            ),
            given_data={"slope_mN_per_lnC": slope, "system_type": system_type, "temperature_K": T},
            gold_answer=gamma_max, tolerance={"rel": 0.02}, source_note=label,
        ))

    for slope, system_type, T, label in GIBBS_CASES * 2:
        i += 1
        gamma_max = gibbs_gamma_max(slope, system_type, T)
        a_min = gibbs_a_min(gamma_max)
        qs.append(Question(
            id=next_id("B", i), category="B", subcategory="gibbs_area_per_molecule",
            tools_required=["gibbs_surface_excess", "gibbs_area_per_molecule"], difficulty="medium",
            trap_type="multi_tool_chain",
            question_text=(
                f"Below the CMC, d(gamma)/d(ln C) = {slope} mN/m for a system with "
                f"system_type={system_type} at {T} K ({label}). First compute Gamma_max "
                f"(mol/m^2), then use it to compute the minimum area per molecule, A_min, in nm^2."
            ),
            given_data={"slope_mN_per_lnC": slope, "system_type": system_type, "temperature_K": T},
            gold_answer=a_min, tolerance={"rel": 0.03}, source_note=label,
        ))

    for C, g0, gmax, K, system_type, T, label in SZYSZKOWSKI_CASES * 3:
        i += 1
        gamma = szyszkowski_surface_tension(C, g0, gmax, K, system_type, T)
        qs.append(Question(
            id=next_id("B", i), category="B", subcategory="szyszkowski_predict_surface_tension",
            tools_required=["szyszkowski_predict_surface_tension"], difficulty="hard",
            question_text=(
                f"Using the Szyszkowski equation with pure-solvent surface tension gamma0={g0} "
                f"mN/m, saturation surface excess Gamma_max={gmax} mol/m^2, Szyszkowski constant "
                f"K={K} (per mM), system_type={system_type}, at {T} K ({label}), predict the "
                f"surface tension (mN/m) at a surfactant concentration of {C} mM."
            ),
            given_data={"concentration": C, "gamma0_mN_m": g0, "gamma_max_mol_per_m2": gmax,
                        "K": K, "temperature_K": T, "system_type": system_type},
            gold_answer=round(gamma, 3), tolerance={"rel": 0.02}, source_note=label,
        ))

    for slope, T, label in NO_SOLUTION_GIBBS_CASES:
        i += 1
        qs.append(Question(
            id=next_id("B", i), category="B", subcategory="gibbs_surface_excess",
            tools_required=["gibbs_surface_excess"], difficulty="hard", trap_type="no_solution",
            question_text=(
                f"Below the CMC, the slope of surface tension vs ln(concentration), "
                f"d(gamma)/d(ln C), is {slope} mN/m for {label}, at {T} K. "
                f"Compute the maximum surface excess concentration, Gamma_max, in mol/m^2."
            ),
            given_data={"slope_mN_per_lnC": slope, "temperature_K": T},
            gold_answer=(
                "cannot compute -- Gamma_max requires system_type (nonionic, "
                "ionic_excess_electrolyte, or ionic_no_added_salt), which depends on the "
                "surfactant's ionic character and whether excess inert electrolyte is "
                "present; this is not stated, so system_type cannot be determined without "
                "guessing"
            ),
            grading_method="category_match", source_note=label,
        ))

    return qs


if __name__ == "__main__":
    import json
    from collections import Counter
    questions = gen()
    print(f"Generated {len(questions)} questions for Category B (adsorption)")
    print(Counter(q.subcategory for q in questions))
