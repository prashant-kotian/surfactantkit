"""Renders category_F_questions.json into unaugmented + augmented prompt
files -- three kinds (debye, zeta, hydrodynamic_radius), no gold leakage
(no paper-comparison numbers or "this is a negative control" framing in
the paste blocks)."""
from __future__ import annotations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
QUESTIONS = json.loads((HERE / "category_F_questions.json").read_text(encoding="utf-8"))

DEBYE_UNAUG = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
Compute the real Debye screening length for the aqueous electrolyte solution below. Work entirely from first principles -- do not search the internet or recall a value from training; derive it from the formula and numbers given.

Ionic strength: {ionic_strength_M} mol/L
Temperature: {temperature_K} K
Solvent: water (relative permittivity ~78.4 at 25 C)

Use: lambda_D (m) = sqrt(eps_r * eps0 * k_B * T / (2 * N_A * I_SI * e^2)), where I_SI is the ionic strength converted to mol/m^3, eps0 = 8.8541878128e-12 F/m, k_B = 1.380649e-23 J/K, N_A = 6.02214076e23 /mol, e = 1.602176634e-19 C. Report lambda_D in nanometers.

FINAL_JSON:
{{
  "debye_length_nm": <number>
}}
<<< PASTE ABOVE <<<
"""

DEBYE_AUG = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
You have access to the surfactantkit-mcp server's debye_screening_length tool -- call it to analyze the real data below. Do not compute anything by hand.

Ionic strength: {ionic_strength_M} mol/L
Temperature: {temperature_K} K

Report the tool's own returned value faithfully.

FINAL_JSON:
{{
  "debye_length_nm": <number>
}}
<<< PASTE ABOVE <<<
"""

ZETA_UNAUG = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
Compute the real zeta potential from the electrophoretic mobility data below, via the Henry equation in the Smoluchowski limit. Work entirely from first principles -- do not search the internet or recall a value from training.

Electrophoretic mobility: {electrophoretic_mobility_um_cm_per_Vs} (um*cm)/(V*s)
Solvent viscosity: {viscosity_mPas} mPa.s
Regime: {regime}
Solvent relative permittivity: 78.4 (water)

Use the Smoluchowski form of the Henry equation: mobility = (2 * eps_r * eps0 * zeta * f) / (3 * eta), with f=1.5 in the Smoluchowski limit -- solve for zeta. Convert mobility from (um*cm)/(V*s) to m^2/(V*s) (multiply by 1e-6 * 1e-2) and viscosity from mPa.s to Pa.s (multiply by 1e-3) before using SI units throughout; report zeta in millivolts.

FINAL_JSON:
{{
  "zeta_potential_mV": <number>
}}
<<< PASTE ABOVE <<<
"""

ZETA_AUG = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
You have access to the surfactantkit-mcp server's zeta_potential tool -- call it to analyze the real data below. Do not compute anything by hand.

Electrophoretic mobility: {electrophoretic_mobility_um_cm_per_Vs} (um*cm)/(V*s)
Solvent viscosity: {viscosity_mPas} mPa.s
Regime: {regime}

Report the tool's own returned value faithfully.

FINAL_JSON:
{{
  "zeta_potential_mV": <number>
}}
<<< PASTE ABOVE <<<
"""

RH_UNAUG = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
Compute the real hydrodynamic radius from the diffusion coefficient below, via the Stokes-Einstein equation. Work entirely from first principles -- do not search the internet or recall a value from training, and do not adjust your answer to match any other hydrodynamic-size measurement you may know of for a similar system -- report exactly what this formula gives for the D actually provided below.

Diffusion coefficient: {diffusion_coefficient_cm2_per_s} cm^2/s
Solvent viscosity: {viscosity_mPas} mPa.s
Temperature: {temperature_K} K

Use: R_h (m) = k_B * T / (6 * pi * eta * D), with D and eta converted to SI units first (D: cm^2/s -> m^2/s, multiply by 1e-4; eta: mPa.s -> Pa.s, multiply by 1e-3). Report R_h in nanometers.

FINAL_JSON:
{{
  "hydrodynamic_radius_nm": <number>
}}
<<< PASTE ABOVE <<<
"""

RH_AUG = """{qid}
[reference only, do not paste -- source: {source}]

>>> PASTE BELOW >>>
You have access to the surfactantkit-mcp server's hydrodynamic_radius tool -- call it to analyze the real data below. Do not compute anything by hand.

Diffusion coefficient: {diffusion_coefficient_cm2_per_s} cm^2/s
Solvent viscosity: {viscosity_mPas} mPa.s
Temperature: {temperature_K} K

Report the tool's own returned value faithfully -- do not adjust it to match any other hydrodynamic-size measurement you may know of for a similar system.

FINAL_JSON:
{{
  "hydrodynamic_radius_nm": <number>
}}
<<< PASTE ABOVE <<<
"""

TEMPLATES = {
    "debye": (DEBYE_UNAUG, DEBYE_AUG),
    "zeta": (ZETA_UNAUG, ZETA_AUG),
    "hydrodynamic_radius": (RH_UNAUG, RH_AUG),
}


def render(q, aug):
    unaug_t, aug_t = TEMPLATES[q["kind"]]
    template = aug_t if aug else unaug_t
    fields = {k: v for k, v in q.items() if isinstance(v, (str, int, float))}
    fields["qid"] = q["id"]
    return template.format(**fields)


def main():
    unaug = "\n\n".join(render(q, False) for q in QUESTIONS)
    (HERE / "CategoryF_Unaugmented.txt").write_text(unaug, encoding="utf-8")
    aug = "\n\n".join(render(q, True) for q in QUESTIONS)
    (HERE / "CategoryF_Augmented.txt").write_text(aug, encoding="utf-8")
    print(f"Rendered {len(QUESTIONS)} questions.")


if __name__ == "__main__":
    main()
