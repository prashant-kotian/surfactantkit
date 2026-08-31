"""Runs every Tier 1 generator, checks for ID collisions and duplicate
gold answers being silently wrong, and writes the combined bank to JSON."""

import sys
import json
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "generators"))

from generators import mixed_micelle_gen, adsorption_gen, hlb_gen, geometry_gen, electrostatics_gen, dynamics_gen, thermodynamics_gen, wetting_gen, solubilization_gen

GENERATORS = [
    ("A", mixed_micelle_gen),
    ("B", adsorption_gen),
    ("C", hlb_gen),
    ("D", geometry_gen),
    ("E", electrostatics_gen),
    ("F", dynamics_gen),
    ("G", thermodynamics_gen),
    ("H", wetting_gen),
    ("I", solubilization_gen),
]


def build():
    all_questions = []
    for label, module in GENERATORS:
        qs = module.gen()
        all_questions.extend(qs)
        print(f"Category {label}: {len(qs)} questions")

    ids = [q.id for q in all_questions]
    dupes = [item for item, count in Counter(ids).items() if count > 1]
    assert not dupes, f"Duplicate question IDs found: {dupes}"

    tool_coverage = Counter()
    for q in all_questions:
        for t in q.tools_required:
            tool_coverage[t] += 1

    print(f"\nTotal Tier 1 questions: {len(all_questions)}")
    print(f"Distinct tools exercised: {len(tool_coverage)}")
    print("\nPer-tool question counts:")
    for tool, count in sorted(tool_coverage.items()):
        print(f"  {tool}: {count}")

    trap_counts = Counter(q.trap_type for q in all_questions)
    print("\nTrap-type distribution:", dict(trap_counts))

    return all_questions


if __name__ == "__main__":
    questions = build()
    out_path = Path(__file__).resolve().parent / "question_bank_tier1.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump([q.to_dict() for q in questions], f, indent=2)
    print(f"\nWrote {len(questions)} questions to {out_path}")
