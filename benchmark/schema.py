"""Shared question schema for all SurfBench generators.

Every Tier 1 question is produced by calling the real SurfactantKit
function to compute the gold answer -- generators never hand-type an
expected numeric result.
"""

from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Any, Literal

TrapType = Literal["none", "unit_trap", "no_solution", "multi_tool_chain"]
GradingMethod = Literal["numeric_tolerance", "exact_match", "category_match"]


@dataclass
class Question:
    id: str
    category: str  # "A".."I" (Tier 1) or "J".."O" (Tier 2)
    subcategory: str  # tool name, or a short label for Tier 2
    tools_required: list[str]
    difficulty: Literal["easy", "medium", "hard"]
    question_text: str
    given_data: dict[str, Any]
    gold_answer: Any
    tolerance: dict[str, float] = field(default_factory=dict)  # e.g. {"rel": 0.02}
    grading_method: GradingMethod = "numeric_tolerance"
    trap_type: TrapType = "none"
    source_note: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


def next_id(category: str, index: int) -> str:
    return f"{category}-{index:03d}"
