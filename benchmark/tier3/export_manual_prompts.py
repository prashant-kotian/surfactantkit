"""Exports Tier 3's question bank as one combined, manual-paste-ready text
file for a Claude/ChatGPT web-UI unaugmented run -- reuses prompts.py's real
build_system()/build_user() so the wording matches every other unaugmented
run in this project exactly (same system instructions, same FINAL_JSON
convention), rather than a fresh ad-hoc format.
"""

from __future__ import annotations
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "harness"))
from prompts import build_system, build_user

BANK_PATH = Path(__file__).resolve().parent.parent / "question_bank_tier3.json"
OUT_PATH = Path(__file__).resolve().parent / "tier3_manual_prompt.txt"


def main():
    questions = json.loads(BANK_PATH.read_text())
    system = build_system("unaugmented")

    separator = "=" * 70
    header = (
        f"{system}\n\n"
        + separator + "\n"
        + "There are 19 questions below, each starting with a line like "
          "'### P-001'. Answer them ALL in this same conversation, in order. "
          "For EACH question, start your answer with its own '### <ID>' line, "
          "then your reasoning, then end that question's answer with its own "
          "'FINAL_JSON: {...}' line before moving to the next question. Do not "
          "skip any question.\n"
        + separator + "\n"
    )

    blocks = [header]
    for q in questions:
        blocks.append(f"\n### {q['id']}\n{build_user(q)}\n")

    OUT_PATH.write_text("".join(blocks), encoding="utf-8")
    print(f"Wrote {len(questions)} questions to {OUT_PATH}")


if __name__ == "__main__":
    main()
