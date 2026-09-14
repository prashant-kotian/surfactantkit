"""Renders category_D_questions.json into unaugmented + augmented prompt
files. Paste block contains only real input data + task instruction --
no design commentary (same discipline applied since the Category B leak)."""
from __future__ import annotations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
QUESTIONS = json.loads((HERE / "category_D_questions.json").read_text(encoding="utf-8"))

UNAUG_TEMPLATE = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
You are analyzing a real surfactant's micellization free energy using the data below. Work entirely from first principles -- do not search the internet, and do not recall or assume a deltaG_mic value for this compound from your own training; derive it from the numbers and formula given.

Compound: {compound_name}
Real critical micelle concentration (CMC): {cmc_mM} mM
Temperature: {temperature_C} C
Real, independently-measured degree of counterion binding (beta), from conductivity: {beta}
Formula to use for this compound's architecture (its own source paper's stated convention, not a generic default -- use this exact formula, not a different counterion_factor convention you may know for a different surfactant architecture): {formula_given}

Compute the counterion_factor implied by this formula, then compute deltaG_mic (kJ/mol) using the real CMC and temperature given. X_cmc = CMC / (CMC + 55.5 mol/L), i.e. mole fraction relative to water's molarity.

End your response with a block in exactly this form:

FINAL_JSON:
{{
  "counterion_factor": <number>,
  "delta_g_mic_kJ_per_mol": <number>,
  "undeterminable": ["short phrase + reason", ...]
}}
<<< PASTE ABOVE <<<
"""

AUG_TEMPLATE = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
You have access to the surfactantkit-mcp server's gibbs_free_energy_micellization tool -- call it to analyze the real data below. Do not compute anything by hand; use the tool, and report exactly what it returns.

Compound: {compound_name}
Real critical micelle concentration (CMC): {cmc_mM} mM
Temperature: {temperature_C} C
Real, independently-measured degree of counterion binding (beta), from conductivity: {beta}
Formula to use for this compound's architecture (its own source paper's stated convention, not a generic default -- use this exact formula, not a different counterion_factor convention you may know for a different surfactant architecture): {formula_given}

Compute the counterion_factor implied by this formula yourself, then pass it to the tool along with the real CMC (converted to mole fraction or however the tool expects it -- check its own docstring) and temperature. Report the tool's own returned deltaG_mic faithfully -- do not override, "correct", or second-guess what it returns.

End your response with a block in exactly this form:

FINAL_JSON:
{{
  "counterion_factor": <number>,
  "delta_g_mic_kJ_per_mol": <number>,
  "undeterminable": ["short phrase + reason", ...]
}}
<<< PASTE ABOVE <<<
"""


def render(template, q):
    return template.format(
        qid=q["id"], source=q["source"], compound_name=q["compound_name"],
        cmc_mM=q["cmc_mM"], temperature_C=q["temperature_C"], beta=q["beta"],
        formula_given=q["formula_given"],
    )


def main():
    unaug = "\n\n".join(render(UNAUG_TEMPLATE, q) for q in QUESTIONS)
    (HERE / "CategoryD_Unaugmented.txt").write_text(unaug, encoding="utf-8")
    aug = "\n\n".join(render(AUG_TEMPLATE, q) for q in QUESTIONS)
    (HERE / "CategoryD_Augmented.txt").write_text(aug, encoding="utf-8")
    print(f"Rendered {len(QUESTIONS)} questions.")


if __name__ == "__main__":
    main()
