"""Category D: molecular geometry. tanford_chain_geometry,
critical_packing_parameter, aggregation_number."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from surfactantkit.cpp import (
    tanford_tail_volume, tanford_critical_length,
    critical_packing_parameter, classify_aggregate_morphology,
    aggregation_number_spherical,
)
from schema import Question, next_id

CHAIN_LENGTHS = [10, 12, 14, 16, 18]

# (n_carbons, head_area_A2, label) -- head areas span spherical-micelle
# to vesicle-forming ionic/gemini regimes
CPP_CASES = [
    (12, 50.0, "typical single-chain ionic surfactant headgroup"),
    (12, 35.0, "tightly packed ionic headgroup (high salt/screening)"),
    (14, 80.0, "loosely packed, bulky headgroup"),
    (16, 45.0, "long-chain ionic surfactant, moderate headgroup area"),
    (18, 30.0, "long-chain surfactant with a small headgroup -- bilayer/vesicle-prone"),
    (12, 100.0, "very bulky headgroup (e.g. large nonionic/gemini with wide spacer)"),
]


def gen() -> list[Question]:
    qs: list[Question] = []
    i = 0

    for nc in CHAIN_LENGTHS * 3:
        i += 1
        v = tanford_tail_volume(nc)
        lc = tanford_critical_length(nc)
        qs.append(Question(
            id=next_id("D", i), category="D", subcategory="tanford_chain_geometry",
            tools_required=["tanford_chain_geometry"], difficulty="easy",
            question_text=(
                f"Using Tanford's formulas, compute the hydrophobic tail volume (cubic "
                f"Angstrom) and the maximum extended (critical) chain length (Angstrom) for a "
                f"saturated, unbranched C{nc} alkyl chain."
            ),
            given_data={"n_carbons": nc},
            gold_answer={"tail_volume_A3": round(v, 2), "critical_length_A": round(lc, 3)},
            tolerance={"rel": 0.01}, source_note=f"C{nc} chain",
        ))

    for nc, head_area, label in CPP_CASES * 2:
        i += 1
        v = tanford_tail_volume(nc)
        lc = tanford_critical_length(nc)
        cpp = critical_packing_parameter(v, head_area, lc)
        morphology = classify_aggregate_morphology(cpp)
        qs.append(Question(
            id=next_id("D", i), category="D", subcategory="critical_packing_parameter",
            tools_required=["tanford_chain_geometry", "critical_packing_parameter"],
            difficulty="hard", trap_type="multi_tool_chain",
            question_text=(
                f"A surfactant has a saturated C{nc} alkyl tail and an optimal headgroup area of "
                f"{head_area} square Angstrom ({label}). First compute the tail volume and "
                f"critical chain length via Tanford's formulas, then the critical packing "
                f"parameter (CPP), then state the predicted aggregate morphology."
            ),
            given_data={"n_carbons": nc, "head_area_A2": head_area},
            gold_answer={"cpp": round(cpp, 4), "morphology": morphology},
            tolerance={"cpp_rel": 0.02}, source_note=label,
        ))

    for nc, _, label in CPP_CASES[:4] * 2:
        i += 1
        v = tanford_tail_volume(nc)
        lc = tanford_critical_length(nc)
        nagg = aggregation_number_spherical(v, lc)
        qs.append(Question(
            id=next_id("D", i), category="D", subcategory="aggregation_number",
            tools_required=["tanford_chain_geometry", "aggregation_number"],
            difficulty="hard", trap_type="multi_tool_chain",
            question_text=(
                f"For a spherical micelle formed by a C{nc} surfactant, approximate the core "
                f"radius as the critical chain length (from Tanford's formula) and use the "
                f"tail volume to estimate the geometric aggregation number."
            ),
            given_data={"n_carbons": nc},
            gold_answer=round(nagg, 2), tolerance={"rel": 0.02}, source_note=label,
        ))

    return qs


if __name__ == "__main__":
    from collections import Counter
    questions = gen()
    print(f"Generated {len(questions)} questions for Category D (geometry)")
    print(Counter(q.subcategory for q in questions))
