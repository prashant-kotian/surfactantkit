"""Category B: adsorption. gibbs_surface_excess, gibbs_area_per_molecule,
szyszkowski_predict_surface_tension."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from surfactantkit.adsorption import gibbs_gamma_max, gibbs_a_min, szyszkowski_surface_tension
from schema import Question, next_id

# (slope_mN_per_lnC, n_factor, temperature_K, label)
GIBBS_CASES = [
    (-5.0, 2.0, 298.15, "ionic surfactant, no added salt (n=2)"),
    (-8.5, 2.0, 298.15, "ionic surfactant, steeper pre-CMC slope"),
    (-3.2, 1.0, 298.15, "nonionic surfactant or ionic with excess electrolyte (n=1)"),
    (-6.7, 1.0, 308.15, "nonionic surfactant, elevated temperature"),
    (-4.4, 2.0, 293.15, "ionic surfactant, lower temperature"),
    (-9.1, 2.0, 298.15, "gemini cationic surfactant, strong adsorption"),
]

SZYSZKOWSKI_CASES = [
    (5.0, 72.0, 3.0e-6, 50.0, 298.15, 1.0, "nonionic, moderate K"),
    (2.0, 72.0, 3.0e-6, 50.0, 298.15, 1.0, "nonionic, lower concentration"),
    (8.0, 68.0, 2.5e-6, 30.0, 298.15, 2.0, "ionic surfactant (n=2), different gamma0"),
    (1.0, 72.0, 4.0e-6, 100.0, 298.15, 1.0, "high-affinity adsorption (large K)"),
]


def gen() -> list[Question]:
    qs: list[Question] = []
    i = 0

    for slope, n, T, label in GIBBS_CASES * 2:
        i += 1
        gamma_max = gibbs_gamma_max(slope, n, T)
        qs.append(Question(
            id=next_id("B", i), category="B", subcategory="gibbs_surface_excess",
            tools_required=["gibbs_surface_excess"], difficulty="medium",
            question_text=(
                f"Below the CMC, the slope of surface tension vs ln(concentration), "
                f"d(gamma)/d(ln C), is {slope} mN/m. This is {label}, at {T} K. "
                f"Compute the maximum surface excess concentration, Gamma_max, in mol/m^2."
            ),
            given_data={"slope_mN_per_lnC": slope, "n_factor": n, "temperature_K": T},
            gold_answer=gamma_max, tolerance={"rel": 0.02}, source_note=label,
        ))

    for slope, n, T, label in GIBBS_CASES * 2:
        i += 1
        gamma_max = gibbs_gamma_max(slope, n, T)
        a_min = gibbs_a_min(gamma_max)
        qs.append(Question(
            id=next_id("B", i), category="B", subcategory="gibbs_area_per_molecule",
            tools_required=["gibbs_surface_excess", "gibbs_area_per_molecule"], difficulty="medium",
            trap_type="multi_tool_chain",
            question_text=(
                f"Below the CMC, d(gamma)/d(ln C) = {slope} mN/m for a system with n_factor={n} "
                f"at {T} K ({label}). First compute Gamma_max (mol/m^2), then use it to compute "
                f"the minimum area per molecule, A_min, in nm^2."
            ),
            given_data={"slope_mN_per_lnC": slope, "n_factor": n, "temperature_K": T},
            gold_answer=a_min, tolerance={"rel": 0.03}, source_note=label,
        ))

    for C, g0, gmax, K, T, n, label in SZYSZKOWSKI_CASES * 3:
        i += 1
        gamma = szyszkowski_surface_tension(C, g0, gmax, K, T, n)
        qs.append(Question(
            id=next_id("B", i), category="B", subcategory="szyszkowski_predict_surface_tension",
            tools_required=["szyszkowski_predict_surface_tension"], difficulty="hard",
            question_text=(
                f"Using the Szyszkowski equation with pure-solvent surface tension gamma0={g0} "
                f"mN/m, saturation surface excess Gamma_max={gmax} mol/m^2, Szyszkowski constant "
                f"K={K} (per mM), n_factor={n}, at {T} K ({label}), predict the surface tension "
                f"(mN/m) at a surfactant concentration of {C} mM."
            ),
            given_data={"concentration": C, "gamma0_mN_m": g0, "gamma_max_mol_per_m2": gmax,
                        "K": K, "temperature_K": T, "n_factor": n},
            gold_answer=round(gamma, 3), tolerance={"rel": 0.02}, source_note=label,
        ))

    return qs


if __name__ == "__main__":
    import json
    from collections import Counter
    questions = gen()
    print(f"Generated {len(questions)} questions for Category B (adsorption)")
    print(Counter(q.subcategory for q in questions))
