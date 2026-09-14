"""Builds the MCP-augmented (tool-available) manual-paste companion to
Round3_Questions_Unaugmented.txt -- same 8 real questions, same underlying
data (reads the already-built round3_questions.json directly rather than
regenerating, so the two conditions are guaranteed to use byte-identical
SMILES/curves/temperatures), different instruction: use the real
SurfactantKit MCP server's ONE new orchestrator tool
(derive_all_properties_from_smiles_and_curve), not the atomic per-property
tools -- matching the 2026-09-13 ROADMAP decision that the SurfMCP-
augmented condition for this round tests the orchestrator specifically,
not a bag of atomic functions a human would have to chain by hand.

Deliberately does NOT tell the model when to pass electrolyte_condition or
counterion_dissociation_alpha to the tool call -- the tool itself only
fixes the ARITHMETIC once the right parameters are supplied; recognizing
prose cues like R3-07's "100 mM NaCl" and translating that into the right
tool argument is still a real judgment call the model has to make itself.
Having the tool available does not remove that trap, by design.

Usage:
  python build_round3_augmented.py
"""
from __future__ import annotations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

QUESTION_TEMPLATE = """R3-{n:02d}
[reference only, do not paste -- data provenance: {provenance}]

>>> PASTE BELOW >>>
You have access to the surfactantkit-mcp server's derive_all_properties_from_smiles_and_curve tool -- call it to analyze the real surfactant below. Do not compute anything by hand; use the tool, and report exactly what it returns. Read its own docstring for what optional parameters it accepts and when passing them is appropriate (e.g. electrolyte_condition, counterion_dissociation_alpha) -- only pass those if this question's own text actually states the relevant real-world condition; do not guess a value for either just to get a more complete-looking answer, and do not search the internet for this compound's condition either.

SMILES: {smiles}
Temperature: {temperature_K} K
Raw tensiometry data (concentration in mM, surface tension in mN/m):
concentration_mM = {concs}
surface_tension_mN_per_m = {gammas}
{narrative}
Call the tool with this SMILES, curve, and temperature. Report the tool's own returned values faithfully -- do not override, "correct", or second-guess a value the tool returns, and do not fill in a value the tool itself left null. Copy the tool's own 'gaps' list into "undeterminable" below (paraphrasing is fine, but don't drop or add gaps the tool didn't itself report).

End your response with a block in exactly this form:

FINAL_JSON:
{{
  "charge_type": "anionic|cationic|nonionic|zwitterionic|ambiguous_pH_dependent",
  "cmc_mM": <number or null>,
  "gamma_max_mol_per_m2": <number or null>,
  "a_min_nm2": <number or null>,
  "isotherm_model": "langmuir" or "frumkin" or null,
  "frumkin_a": <number or null, only meaningful if isotherm_model is "frumkin">,
  "delta_g_mic_kJ_per_mol": <number or null>,
  "undeterminable": ["short phrase + reason", ...]
}}
<<< PASTE ABOVE <<<
"""


def main():
    questions = json.loads((HERE / "round3_questions.json").read_text(encoding="utf-8"))

    blocks = []
    for i, q in enumerate(questions, 1):
        narrative = (q["narrative_extra"] + "\n") if q["narrative_extra"] else ""
        blocks.append(QUESTION_TEMPLATE.format(
            n=i, provenance=q["data_provenance"], smiles=q["smiles"], temperature_K=q["temperature_K"],
            concs=q["concentrations_mM"], gammas=[round(g, 4) for g in q["surface_tensions_mN_per_m"]],
            narrative=narrative,
        ))
    text = "\n\n".join(blocks)
    (HERE / "Round3_Questions_Augmented.txt").write_text(text, encoding="utf-8")
    print(f"Built {len(questions)} augmented-condition question prompts -> Round3_Questions_Augmented.txt")


if __name__ == "__main__":
    main()
