"""Renders category_C_questions.json into unaugmented + augmented prompt
files. Lesson applied directly from the Category B leak bug: the prompt
block contains ONLY real input data and the task instruction -- no design
commentary, no hint about which case is the "trap" or "the messy one"."""
from __future__ import annotations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
QUESTIONS = json.loads((HERE / "category_C_questions.json").read_text(encoding="utf-8"))

UNAUG_TEMPLATE = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
You are analyzing a real surfactant's micellization thermodynamics using the data below. Work entirely from first principles -- do not search the internet, and do not recall or assume a deltaH_mic value for this compound from your own training; derive it from the real van't Hoff (Gibbs-Helmholtz) relationship applied to the numbers given.

Compound: {compound_name}
SMILES: {smiles}
Real CMC measured at each of the following temperatures (all from the same source):
{temp_block}

Compute the enthalpy of micellization (deltaH_mic, kJ/mol) using the WIDEST available temperature span (the lowest and highest temperature given). If 3 or more temperatures are given, ALSO compute deltaH_mic for each adjacent pair of temperatures separately (sub-intervals), and assess whether these sub-interval values are reasonably self-consistent with each other (i.e. whether the "deltaH_mic is constant over this temperature range" assumption a two-point estimate always makes actually holds here) -- do not claim a self-consistency assessment if only 2 temperatures are given, since no sub-interval comparison is possible with 2 points.

End your response with a block in exactly this form (use null for anything not genuinely determined):

FINAL_JSON:
{{
  "delta_h_endpoint_kJ_per_mol": <number or null>,
  "endpoint_temps_C": [<lower>, <upper>],
  "sub_intervals": [{{"t_lo_C": <number>, "t_hi_C": <number>, "delta_h_kJ_per_mol": <number>}}, ...] or null,
  "self_consistent": true, false, or null,
  "undeterminable": ["short phrase + reason", ...]
}}
<<< PASTE ABOVE <<<
"""

AUG_TEMPLATE = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
You have access to the surfactantkit-mcp server's vant_hoff_enthalpy tool -- call it to analyze the real data below. Do not compute anything by hand; use the tool, and report exactly what it returns.

Compound: {compound_name}
SMILES: {smiles}
Real CMC measured at each of the following temperatures (all from the same source):
{temp_block}

Call vant_hoff_enthalpy(cmc1_M, temperature1_K, cmc2_M, temperature2_K) -- CMC in molarity (mol/L, i.e. mM/1000), temperature in Kelvin (C + 273.15) -- once using the lowest and highest given temperature (the endpoint value), and ALSO once per adjacent pair of temperatures if 3 or more are given (sub-intervals). Do not attempt a sub-interval calculation if only 2 temperatures are given. Compare the sub-interval values (if computed) to assess whether they are reasonably self-consistent with each other. Report the tool's own returned values faithfully -- do not override, "correct", or second-guess what it returns.

End your response with a block in exactly this form:

FINAL_JSON:
{{
  "delta_h_endpoint_kJ_per_mol": <number or null>,
  "endpoint_temps_C": [<lower>, <upper>],
  "sub_intervals": [{{"t_lo_C": <number>, "t_hi_C": <number>, "delta_h_kJ_per_mol": <number>}}, ...] or null,
  "self_consistent": true, false, or null,
  "undeterminable": ["short phrase + reason", ...]
}}
<<< PASTE ABOVE <<<
"""


def render(template, q):
    temp_block = "\n".join(f"  T = {t} C, CMC = {c} mM" for t, c in sorted(q["given_temps_C_cmc_mM"].items(), key=lambda kv: float(kv[0])))
    return template.format(qid=q["id"], source=q["source"], compound_name=q["compound_name"],
                            smiles=q["smiles"], temp_block=temp_block)


def main():
    unaug = "\n\n".join(render(UNAUG_TEMPLATE, q) for q in QUESTIONS)
    (HERE / "CategoryC_Unaugmented.txt").write_text(unaug, encoding="utf-8")
    aug = "\n\n".join(render(AUG_TEMPLATE, q) for q in QUESTIONS)
    (HERE / "CategoryC_Augmented.txt").write_text(aug, encoding="utf-8")
    print(f"Rendered {len(QUESTIONS)} questions.")


if __name__ == "__main__":
    main()
