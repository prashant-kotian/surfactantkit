"""Renders category_B_questions.json into unaugmented + augmented manual-
paste prompt files, same PASTE-marker convention as Round 3/Category A."""
from __future__ import annotations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
QUESTIONS = json.loads((HERE / "category_B_questions.json").read_text(encoding="utf-8"))

FIELD_DESCRIPTIONS = {
    "clint_ideal_cmc_mM": '"clint_ideal_cmc_mM": <number or null>',
    "rubingh_x1": '"rubingh_x1": <number or null, the REGULAR-SOLUTION MICELLAR MOLE FRACTION OF COMPONENT 1>',
    "rubingh_x1_component1": '"rubingh_x1": <number or null, the REGULAR-SOLUTION MICELLAR MOLE FRACTION OF COMPONENT 1 SPECIFICALLY -- not component 2>',
    "rubingh_beta": '"rubingh_beta": <number or null, the Rubingh interaction parameter>',
    "synergy_classification": '"synergy_classification": "synergistic" or "antagonistic" or "ideal" or null',
    "literature_validated_fields": '"literature_validated_fields": [<list of field names above that you believe ARE independently literature-cross-checkable for this specific system, vs. tool-only>]',
}


def schema_block(ask):
    lines = [FIELD_DESCRIPTIONS[a] for a in ask]
    lines.append('"undeterminable": ["short phrase + reason", ...]')
    return "{\n  " + ",\n  ".join(lines) + "\n}"


UNAUG_TEMPLATE = """{qid}
[reference only, do not paste -- system: {system}, source: {source}]

>>> PASTE BELOW >>>
You are analyzing a real binary surfactant mixture using the data below. Work entirely from first principles -- do not search the internet, and do not recall or assume a literature value for any of the quantities asked below for THIS specific system; derive everything from the numbers given.

Component 1: {c1_name} (pure CMC = {c1_cmc} mM)
Component 2: {c2_name} (pure CMC = {c2_cmc} mM)
Bulk mole fraction of component 1 (alpha1): {alpha1}
Experimentally measured mixed CMC: {cmc_mix} mM
{prompt_note}
Determine: {ask_list}

End your response with a block in exactly this form (use null for anything not genuinely determined):

FINAL_JSON:
{schema}
<<< PASTE ABOVE <<<
"""

AUG_TEMPLATE = """{qid}
[reference only, do not paste -- system: {system}, source: {source}]

>>> PASTE BELOW >>>
You have access to the surfactantkit-mcp server's clint_ideal_cmc and rubingh_solve tools -- call whichever are relevant to analyze the real binary mixture below. Do not compute anything by hand; use the tool(s), and report exactly what they return.

Component 1: {c1_name} (pure CMC = {c1_cmc} mM)
Component 2: {c2_name} (pure CMC = {c2_cmc} mM)
Bulk mole fraction of component 1 (alpha1): {alpha1}
Experimentally measured mixed CMC: {cmc_mix} mM
{prompt_note}
Determine: {ask_list}

Call rubingh_solve(alpha1, cmc_mix_mM, cmc1_mM, cmc2_mM) with component 1's pure CMC as cmc1_mM -- its returned micellar_mole_fraction_x1 is ALWAYS for component 1 as you passed it, regardless of how any literature source you might recall frames it. Report the tool's own returned values faithfully -- do not override, "correct", or second-guess what it returns.

End your response with a block in exactly this form:

FINAL_JSON:
{schema}
<<< PASTE ABOVE <<<
"""

ASK_LABELS = {
    "clint_ideal_cmc_mM": "the Clint ideal mixed CMC",
    "rubingh_x1": "the Rubingh micellar mole fraction of component 1 (x1) and the interaction parameter beta",
    "rubingh_x1_component1": "the Rubingh micellar mole fraction of component 1 SPECIFICALLY (not component 2)",
    "rubingh_beta": None,  # covered by rubingh_x1's label
    "synergy_classification": "whether the mixture is synergistic, antagonistic, or ideal",
    "literature_validated_fields": "which of your own reported fields are independently literature-cross-checkable for this exact system",
}


def ask_list_text(ask):
    labels = [ASK_LABELS[a] for a in ask if ASK_LABELS.get(a)]
    return "; ".join(labels)


def render(template, q):
    ask = q["ask"]
    prompt_note = q.get("prompt_note")
    prompt_note_block = f"\n{prompt_note}\n" if prompt_note else ""
    return template.format(
        qid=q["id"], system=q["system"], source=q["source"],
        c1_name=q["component1"]["name"], c1_cmc=q["component1"]["pure_cmc_mM"],
        c2_name=q["component2"]["name"], c2_cmc=q["component2"]["pure_cmc_mM"],
        alpha1=q["alpha1"], cmc_mix=q["cmc_mix_mM"],
        prompt_note=prompt_note_block, ask_list=ask_list_text(ask),
        schema=schema_block(ask),
    )


def main():
    unaug = "\n\n".join(render(UNAUG_TEMPLATE, q) for q in QUESTIONS)
    (HERE / "CategoryB_Unaugmented.txt").write_text(unaug, encoding="utf-8")
    aug = "\n\n".join(render(AUG_TEMPLATE, q) for q in QUESTIONS)
    (HERE / "CategoryB_Augmented.txt").write_text(aug, encoding="utf-8")
    print(f"Rendered {len(QUESTIONS)} questions to CategoryB_Unaugmented.txt / CategoryB_Augmented.txt")


if __name__ == "__main__":
    main()
