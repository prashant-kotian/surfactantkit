"""Renders category_H_questions.json -- same template pattern as
Category G (identical task type), separate category because the source
compound deliberately overlaps with Round 3's R3-04. No gold leakage and
no mention of the cross-reference in the paste block."""
from __future__ import annotations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
QUESTIONS = json.loads((HERE / "category_H_questions.json").read_text(encoding="utf-8"))

UNAUG_TEMPLATE = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
Compute the real work of adhesion and spreading coefficient for the liquid-solid (glass) system below, using the Young-Dupre equation. Work entirely from first principles -- do not search the internet or recall a value from training; derive it from the numbers given.

Compound: {compound_name}
Solvent system: {solvent_system}
Real liquid surface tension (gamma_LV): {gamma_LV_mN_per_m} mN/m
Real measured contact angle on glass (theta): {contact_angle_deg} degrees

Use: work of adhesion W_A (mJ/m^2) = gamma_LV * (1 + cos(theta)); spreading coefficient S (mJ/m^2) = gamma_LV * (cos(theta) - 1).

End your response with a block in exactly this form:

FINAL_JSON:
{{
  "work_of_adhesion_mJ_per_m2": <number>,
  "spreading_coefficient_mJ_per_m2": <number>
}}
<<< PASTE ABOVE <<<
"""

AUG_TEMPLATE = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
You have access to the surfactantkit-mcp server's wetting_work_of_adhesion and wetting_spreading_coefficient tools -- call them to analyze the real data below. Do not compute anything by hand.

Compound: {compound_name}
Solvent system: {solvent_system}
Real liquid surface tension (gamma_LV): {gamma_LV_mN_per_m} mN/m
Real measured contact angle on glass (theta): {contact_angle_deg} degrees

Report the tools' own returned values faithfully.

End your response with a block in exactly this form:

FINAL_JSON:
{{
  "work_of_adhesion_mJ_per_m2": <number>,
  "spreading_coefficient_mJ_per_m2": <number>
}}
<<< PASTE ABOVE <<<
"""


def render(template, q):
    return template.format(qid=q["id"], source=q["source"], compound_name=q["compound_name"],
                            solvent_system=q["solvent_system"], gamma_LV_mN_per_m=q["gamma_LV_mN_per_m"],
                            contact_angle_deg=q["contact_angle_deg"])


def main():
    unaug = "\n\n".join(render(UNAUG_TEMPLATE, q) for q in QUESTIONS)
    (HERE / "CategoryH_Unaugmented.txt").write_text(unaug, encoding="utf-8")
    aug = "\n\n".join(render(AUG_TEMPLATE, q) for q in QUESTIONS)
    (HERE / "CategoryH_Augmented.txt").write_text(aug, encoding="utf-8")
    print(f"Rendered {len(QUESTIONS)} questions.")


if __name__ == "__main__":
    main()
