"""Category I: solubilization. molar_solubilization_ratio."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from surfactantkit.solubilization import molar_solubilization_ratio
from schema import Question, next_id

# (total_solubilized_M, water_solubility_M, surfactant_conc_M, cmc_M, label)
CASES = [
    (0.51e-3, 0.01e-3, 10e-3, 2e-3, "PAH solubilization above CMC"),
    (0.20e-3, 0.005e-3, 5e-3, 1e-3, "moderate solubilization capacity"),
    (1.0e-3, 0.02e-3, 20e-3, 3e-3, "high solubilization capacity gemini system"),
    (0.08e-3, 0.01e-3, 8e-3, 2e-3, "weak solubilization enhancement"),
    (0.35e-3, 0.015e-3, 15e-3, 4e-3, "bile-salt mixed micelle solubilization"),
]

# deliberate trap: surfactant concentration below CMC -- no micelles, MSR undefined
NO_SOLUTION_CASES = [
    (0.05e-3, 0.01e-3, 1e-3, 2e-3, "surfactant concentration below CMC"),
]


def gen() -> list[Question]:
    qs: list[Question] = []
    i = 0

    for total, water_sol, surf_c, cmc, label in CASES * 3:
        i += 1
        msr = molar_solubilization_ratio(total, water_sol, surf_c, cmc)
        qs.append(Question(
            id=next_id("I", i), category="I", subcategory="molar_solubilization_ratio",
            tools_required=["molar_solubilization_ratio"], difficulty="medium",
            question_text=(
                f"A hydrophobic solubilizate has total solubilized concentration "
                f"{total*1000:.4f} mM in a surfactant solution at {surf_c*1000:.2f} mM "
                f"surfactant concentration ({label}). The intrinsic water solubility of the "
                f"solubilizate is {water_sol*1000:.4f} mM, and the surfactant's CMC is "
                f"{cmc*1000:.2f} mM. Compute the Molar Solubilization Ratio (MSR)."
            ),
            given_data={"total_solubilized_M": total, "intrinsic_water_solubility_M": water_sol,
                        "surfactant_concentration_M": surf_c, "cmc_M": cmc},
            gold_answer=round(msr, 4), tolerance={"rel": 0.01}, source_note=label,
        ))

    for total, water_sol, surf_c, cmc, label in NO_SOLUTION_CASES * 3:
        i += 1
        qs.append(Question(
            id=next_id("I", i), category="I", subcategory="molar_solubilization_ratio",
            difficulty="hard", trap_type="no_solution", tools_required=["molar_solubilization_ratio"],
            question_text=(
                f"A surfactant solution at {surf_c*1000:.2f} mM (below its CMC of "
                f"{cmc*1000:.2f} mM, {label}) shows {total*1000:.4f} mM total solubilized "
                f"hydrophobic compound (intrinsic water solubility {water_sol*1000:.4f} mM). "
                f"Compute the Molar Solubilization Ratio."
            ),
            given_data={"total_solubilized_M": total, "intrinsic_water_solubility_M": water_sol,
                        "surfactant_concentration_M": surf_c, "cmc_M": cmc},
            gold_answer="undefined -- surfactant concentration is below the CMC, so no micelles exist to solubilize into; MSR is not meaningful here",
            grading_method="category_match", source_note=label,
        ))

    return qs


if __name__ == "__main__":
    from collections import Counter
    questions = gen()
    print(f"Generated {len(questions)} questions for Category I (solubilization)")
    print(Counter(q.trap_type for q in questions))
