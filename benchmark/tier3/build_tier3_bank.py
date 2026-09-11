"""Orchestrator for Tier 3 (reliability/robustness categories P/Q/R/S):
runs reliability_gen.py, checks for ID collisions, and writes the bank to
question_bank_tier3.json. Mirrors build_tier2_bank.py's verification
discipline, scaled to this bank's much smaller, deliberately curated size.
"""

from __future__ import annotations
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "generators"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from reliability_gen import gen


def main():
    all_questions = gen()
    total = len(all_questions)
    print(f"Total Tier 3 questions: {total}")

    by_cat = Counter(q.category for q in all_questions)
    print("Per-category counts:", dict(by_cat))
    for cat in "PQRS":
        assert by_cat.get(cat, 0) > 0, f"category {cat} produced zero questions"

    ids = [q.id for q in all_questions]
    dupes = [i for i, c in Counter(ids).items() if c > 1]
    assert not dupes, f"duplicate question IDs found: {dupes}"
    print("No duplicate IDs.")

    # Category-match string gold answers routed through a DEDICATED trap_type
    # branch in grading.py (not the generic strict-substring category_match
    # path) are exempt from the 8-word constraint documented in
    # build_tier2_bank.py -- verify every string-gold question here actually
    # has a dedicated trap_type branch, so a future new trap_type added
    # without updating grading.py fails loudly here instead of silently
    # grading wrong later.
    DEDICATED_TRAP_TYPES = {"no_solution", "applicability_domain_violation"}
    for q in all_questions:
        if isinstance(q.gold_answer, str) and q.grading_method == "category_match":
            assert q.trap_type in DEDICATED_TRAP_TYPES, (
                f"{q.id}: string gold answer with grading_method=category_match but "
                f"trap_type={q.trap_type!r} has no dedicated grading.py branch -- it "
                f"would fall through to the generic strict-substring path, which "
                f"requires <=8 words in gold_answer (this one has "
                f"{len(q.gold_answer.split())})"
            )
    print("All category_match string gold answers route through a dedicated grading.py branch.")

    out_path = Path(__file__).resolve().parent.parent / "question_bank_tier3.json"
    out_path.write_text(json.dumps([q.to_dict() for q in all_questions], indent=2))
    print(f"\nWrote {total} questions to {out_path}")

    print("\nTrap type distribution:", dict(Counter(q.trap_type for q in all_questions)))
    print("Difficulty distribution:", dict(Counter(q.difficulty for q in all_questions)))


if __name__ == "__main__":
    main()
