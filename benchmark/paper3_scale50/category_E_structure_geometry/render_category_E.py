"""Renders category_E_questions.json into unaugmented + augmented prompt
files. Two question kinds (cpp_chain, hlb_griffin) get different templates.
Paste block contains only real input data + the real, universal formula
constants (Tanford's regression coefficients, the EO unit mass) -- these
are fixed textbook/literature constants, not a per-paper convention choice,
so supplying them is fair scoping, not a gold leak. No design commentary,
no paper-comparison numbers anywhere in the paste block."""
from __future__ import annotations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
QUESTIONS = json.loads((HERE / "category_E_questions.json").read_text(encoding="utf-8"))

CPP_UNAUG = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
You are analyzing a real surfactant's packing geometry using the data below. Work entirely from first principles -- do not search the internet, and do not recall or assume a critical packing parameter (CPP) or morphology for this compound from your own training; derive it from the numbers and formulas given.

Compound: {compound_name}
Alkyl chain length: {n_carbons} carbons
Real, independently-measured optimal headgroup area (from Gibbs adsorption isotherm surface-tension data): A_min = {a_min_A2} A^2
Temperature: {temperature_K} K

Use Tanford's standard empirical formulas:
  tail volume V (A^3) = 27.4 + 26.9 * n_carbons
  critical chain length l_c (A) = 1.5 + 1.265 * n_carbons
  critical packing parameter P = V / (A_min * l_c)

Classify the predicted aggregate morphology from P using the standard thresholds: P <= 1/3 -> spherical micelle; 1/3 < P <= 1/2 -> cylindrical/rodlike micelle; 1/2 < P <= 1 -> vesicle/bilayer; P > 1 -> inverted structures.

End your response with a block in exactly this form:

FINAL_JSON:
{{
  "tail_volume_A3": <number>,
  "critical_length_A": <number>,
  "cpp": <number>,
  "predicted_morphology": "<short string>",
  "undeterminable": ["short phrase + reason", ...]
}}
<<< PASTE ABOVE <<<
"""

CPP_AUG = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
You have access to the surfactantkit-mcp server's tanford_chain_geometry and critical_packing_parameter tools -- call them to analyze the real data below. Do not compute anything by hand; use the tools, and report exactly what they return.

Compound: {compound_name}
Alkyl chain length: {n_carbons} carbons
Real, independently-measured optimal headgroup area (from Gibbs adsorption isotherm surface-tension data): A_min = {a_min_A2} A^2
Temperature: {temperature_K} K

Call tanford_chain_geometry(n_carbons) first to get the real tail volume and critical length, then pass those (with A_min above) to critical_packing_parameter to get the CPP and predicted morphology. Report the tools' own returned values faithfully -- do not override, "correct", or second-guess what they return.

End your response with a block in exactly this form:

FINAL_JSON:
{{
  "tail_volume_A3": <number>,
  "critical_length_A": <number>,
  "cpp": <number>,
  "predicted_morphology": "<short string>",
  "undeterminable": ["short phrase + reason", ...]
}}
<<< PASTE ABOVE <<<
"""

HLB_UNAUG = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
You are analyzing a real nonionic surfactant's HLB (hydrophilic-lipophilic balance) using the data below. Work entirely from first principles -- do not search the internet, and do not recall or assume an HLB value for this compound from your own training; derive it from the numbers and formula given.

Compound: {compound_name}
Hydrophobe (non-ethoxylated portion) molecular weight: {hydrophobe_mw_g_per_mol} g/mol
Number of ethylene oxide (EO) units: {n_eo_units}
Real molar mass of one EO addition unit: 44.053 g/mol

Use Griffin's method: HLB = 20 * (M_hydrophilic / M_total), where M_hydrophilic = n_EO * 44.053, and M_total = hydrophobe MW + M_hydrophilic.

End your response with a block in exactly this form:

FINAL_JSON:
{{
  "mw_hydrophilic_g_per_mol": <number>,
  "mw_total_g_per_mol": <number>,
  "hlb": <number>,
  "undeterminable": ["short phrase + reason", ...]
}}
<<< PASTE ABOVE <<<
"""

HLB_AUG = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
You have access to the surfactantkit-mcp server's hlb_from_mw tool -- call it to analyze the real data below. Do not compute anything by hand; use the tool, and report exactly what it returns.

Compound: {compound_name}
Hydrophobe (non-ethoxylated portion) molecular weight: {hydrophobe_mw_g_per_mol} g/mol
Number of ethylene oxide (EO) units: {n_eo_units}
Real molar mass of one EO addition unit: 44.053 g/mol

Compute M_hydrophilic = n_EO * 44.053 and M_total = hydrophobe MW + M_hydrophilic yourself, then call hlb_from_mw with those two values. Report the tool's own returned HLB faithfully -- do not override, "correct", or second-guess what it returns.

End your response with a block in exactly this form:

FINAL_JSON:
{{
  "mw_hydrophilic_g_per_mol": <number>,
  "mw_total_g_per_mol": <number>,
  "hlb": <number>,
  "undeterminable": ["short phrase + reason", ...]
}}
<<< PASTE ABOVE <<<
"""


def render(q, aug):
    if q["kind"] == "cpp_chain":
        template = CPP_AUG if aug else CPP_UNAUG
        return template.format(qid=q["id"], source=q["source"], compound_name=q["compound_name"],
                                n_carbons=q["n_carbons"], a_min_A2=q["a_min_A2"], temperature_K=q["temperature_K"])
    else:
        template = HLB_AUG if aug else HLB_UNAUG
        return template.format(qid=q["id"], source=q["source"], compound_name=q["compound_name"],
                                hydrophobe_mw_g_per_mol=q["hydrophobe_mw_g_per_mol"], n_eo_units=q["n_eo_units"])


def main():
    unaug = "\n\n".join(render(q, False) for q in QUESTIONS)
    (HERE / "CategoryE_Unaugmented.txt").write_text(unaug, encoding="utf-8")
    aug = "\n\n".join(render(q, True) for q in QUESTIONS)
    (HERE / "CategoryE_Augmented.txt").write_text(aug, encoding="utf-8")
    print(f"Rendered {len(QUESTIONS)} questions.")


if __name__ == "__main__":
    main()
