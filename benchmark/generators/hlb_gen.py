"""Category C: HLB. hlb_from_mw, hlb_from_groups."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from surfactantkit.hlb import hlb_griffin, hlb_davies
from schema import Question, next_id

GRIFFIN_CASES = [
    (400.0, 600.0, "a nonionic ethoxylate with a moderate EO chain"),
    (150.0, 450.0, "a nonionic surfactant with a short EO chain"),
    (550.0, 700.0, "a highly ethoxylated nonionic surfactant"),
    (80.0, 350.0, "a low-HLB sorbitan ester"),
    (300.0, 320.0, "an almost entirely hydrophilic small molecule"),
    (20.0, 400.0, "a mostly lipophilic surfactant"),
]

# Davies group-count sets built purely from the verified group table
DAVIES_CASES = [
    ({"SO4Na": 1, "CH2": 11, "CH3": 1}, "SDS (dodecyl sulfate, sodium salt)"),
    ({"COONa": 1, "CH2": 10, "CH3": 1}, "sodium laurate-like soap"),
    ({"O_ether": 4, "OH_free": 1, "CH2": 11, "CH3": 1}, "a short-EO nonionic ethoxylate"),
    ({"ester_free": 1, "OH_free": 2, "CH2": 15, "CH3": 1}, "a glycerol monoester-like nonionic"),
    ({"COOH": 1, "CH2": 15, "CH3": 1}, "a free fatty acid"),
]

# Deliberate trap: a group with no verified Davies number (quaternary
# ammonium) -- gold answer is "cannot compute, group not verified"
NO_SOLUTION_DAVIES = [
    ({"quaternary_ammonium": 1, "CH2": 11, "CH3": 1}, "a quaternary ammonium cationic surfactant"),
    ({"amide": 1, "CH2": 11, "CH3": 1}, "an amide-linked surfactant"),
]


def gen() -> list[Question]:
    qs: list[Question] = []
    i = 0

    for mw_h, mw_t, label in GRIFFIN_CASES * 2:
        i += 1
        hlb = hlb_griffin(mw_h, mw_t)
        qs.append(Question(
            id=next_id("C", i), category="C", subcategory="hlb_from_mw",
            tools_required=["hlb_from_mw"], difficulty="easy",
            question_text=(
                f"A surfactant molecule has a total molecular weight of {mw_t} g/mol, of which "
                f"{mw_h} g/mol belongs to the hydrophilic portion ({label}). Using Griffin's "
                f"method, compute the HLB value."
            ),
            given_data={"mw_hydrophilic": mw_h, "mw_total": mw_t},
            gold_answer=round(hlb, 3), tolerance={"rel": 0.02}, source_note=label,
        ))

    for groups, label in DAVIES_CASES * 2:
        i += 1
        hlb = hlb_davies(groups)
        qs.append(Question(
            id=next_id("C", i), category="C", subcategory="hlb_from_groups",
            tools_required=["hlb_from_groups"], difficulty="medium",
            question_text=(
                f"Using Davies' group-contribution method, compute the HLB of a surfactant with "
                f"the following structural groups ({label}): {groups}. "
                f"(HLB = 7 + sum(hydrophilic numbers) - sum(lipophilic numbers).)"
            ),
            given_data={"group_counts": groups},
            gold_answer=round(hlb, 3), tolerance={"rel": 0.02}, source_note=label,
        ))

    for groups, label in NO_SOLUTION_DAVIES:
        i += 1
        qs.append(Question(
            id=next_id("C", i), category="C", subcategory="hlb_from_groups",
            tools_required=["hlb_from_groups"], difficulty="hard", trap_type="no_solution",
            question_text=(
                f"Using Davies' group-contribution method, compute the HLB of {label} with "
                f"structural groups: {groups}."
            ),
            given_data={"group_counts": groups},
            gold_answer="cannot compute -- one or more groups (e.g. quaternary ammonium, amide) have no verified Davies group number in the literature; do not guess a value",
            grading_method="category_match",
            source_note="deliberately includes a group with no cited Davies number",
        ))

    return qs


if __name__ == "__main__":
    from collections import Counter
    questions = gen()
    print(f"Generated {len(questions)} questions for Category C (HLB)")
    print(Counter(q.subcategory for q in questions))
