"""Orchestrator for Tier 2 (non-tool control categories): runs all 6 hand-
curated category modules, checks for ID collisions, and writes the combined
bank to question_bank_tier2.json. Mirrors build_question_bank.py's discipline
for Tier 1 (verify, don't just assume, before writing the final file)."""

from __future__ import annotations
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import category_J_classification as cat_j
import category_K_property_recall as cat_k
import category_L_unit_conversion as cat_l
import category_M_conceptual_theory as cat_m
import category_N_gemini_structural as cat_n
import category_O_real_world_application as cat_o

MODULES = [
    ("J", cat_j, 40), ("K", cat_k, 40), ("L", cat_l, 30),
    ("M", cat_m, 40), ("N", cat_n, 30), ("O", cat_o, 26),
]


def main():
    all_questions = []
    print("Per-category counts:")
    for label, mod, target in MODULES:
        qs = mod.QUESTIONS
        status = "OK" if len(qs) == target else f"MISMATCH (target {target})"
        print(f"  Category {label}: {len(qs)} questions -- {status}")
        all_questions.extend(qs)

    total = len(all_questions)
    print(f"\nTotal Tier 2 questions: {total}")
    assert total == 206, f"expected 206 total, got {total}"

    # ID collision check
    ids = [q.id for q in all_questions]
    dupes = [i for i, c in Counter(ids).items() if c > 1]
    assert not dupes, f"duplicate question IDs found: {dupes}"
    print("No duplicate IDs.")

    # gold-answer grading-safety check: no string gold answer over 8 words
    # (a longer one breaks category_match's substring-containment grading --
    # this was a real bug caught and fixed twice while building this bank)
    long_answers = [q.id for q in all_questions
                    if isinstance(q.gold_answer, str) and len(q.gold_answer.split()) > 8]
    assert not long_answers, f"gold answers too long for reliable grading: {long_answers}"
    print("All string gold answers are grading-safe length (<=8 words).")

    # every Tier 2 question must have tools_required == [] (that's the whole point)
    tooled = [q.id for q in all_questions if q.tools_required]
    assert not tooled, f"Tier 2 questions must have no tools_required, found some on: {tooled}"
    print("Confirmed: zero tools_required across all Tier 2 questions (as designed).")

    out_path = Path(__file__).resolve().parent.parent / "question_bank_tier2.json"
    out_path.write_text(json.dumps([q.to_dict() for q in all_questions], indent=2))
    print(f"\nWrote {total} questions to {out_path}")

    print("\nGrading method distribution:", dict(Counter(q.grading_method for q in all_questions)))
    print("Difficulty distribution:", dict(Counter(q.difficulty for q in all_questions)))


if __name__ == "__main__":
    main()
