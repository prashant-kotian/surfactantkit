"""Round 3 pilot: the design the user actually asked for after three
corrective rounds. Each question gives ONLY a surfactant's real SMILES and
a raw surface-tension-vs-concentration curve -- no stated ionic character,
no named isotherm model, no pre-fit parameters, no pre-known CMC. The
model must derive everything it can determine and explicitly name what it
cannot, the same way SurfactantKit's own orchestrator
(derive_all_properties_from_smiles_and_curve) does. Free-form reasoning is
requested in the question text itself (per the user's explicit instruction,
2026-09-13: "fully derive everything is what i want"); a canonical
FINAL_JSON block is required only at the end, for grading, using the exact
key names the orchestrator itself returns.

Gold answers are computed by ACTUALLY CALLING the orchestrator on each
question's raw data -- never hand-derived -- so this doubles as a further
end-to-end validation of that pipeline across a deliberately varied real
input set (data source discipline below).

Data sourcing per question (every question's provenance is disclosed
honestly, not just asserted):
  Q1-Q3: REAL raw robotic pendant-drop (concentration, surface tension)
    points, Dankloff et al. 2025 (PendantProp), already used and gold-
    verified elsewhere in this project (SDS, DTAB, CTAB).
  Q4: REAL raw tensiometry points, Shah/Das/Bhattarai 2025 (Heliyon,
    PMC11835642), AOT -- already the literature-validated dataset behind
    this project's own cmc_from_surface_tension_curve regression tests.
  Q5: ILLUSTRATIVE curve for a real nonionic compound (C12E8), generated
    from this codebase's own szyszkowski_surface_tension with
    representative (not literature-transcribed) K/Gamma_max -- disclosed
    as illustrative because reconstructing an exact raw curve from the
    only real nonionic Frumkin sources this project has mined (Zawala
    2020's MIBC, a weak/atypical surfactant with only a ~5 mN/m total
    premicellar decline; Taylor/Valkovska/Bain 2003's C12E5/C12E3/C10E8,
    missing the 'a' parameter needed to reconstruct their exact curve)
    would either produce a pedagogically confusing curve or require
    re-deriving an unverifiable conversion. Included specifically to test
    correct nonionic classification and Langmuir (a=0) selection, not to
    benchmark against a literature-exact CMC.
  Q6: ILLUSTRATIVE curve for a real zwitterionic compound
    (cocamidopropyl betaine), same disclosed-illustrative basis as Q5.
    Included specifically to test the REFUSAL behavior (Gamma_max/A_min/
    isotherm/deltaG_mic must all come back null/undetermined) -- this
    does not depend on numeric literature accuracy at all, only on
    whether the charge classification and the downstream refusal logic
    fire correctly.
  Q7: REAL fitted K/Gamma_max (Prosser & Franses 2001, Table 1, SDS in
    100 mM NaCl -- already validated in this project's own
    szyszkowski_fit_K tests) used to generate an ILLUSTRATIVE curve (raw
    (C,gamma) pairs were never published, only the fitted K), extended
    with a disclosed-illustrative postmicellar plateau near the pre-CMC
    curve's own value at the real literature CMC region for this system.
    Included specifically to test the EXACT systematic error found in
    this project's own v1 unaugmented pilot: both Claude and ChatGPT
    silently used the no-added-salt Gibbs prefactor (n=2) instead of the
    excess-electrolyte one (n=1) when the electrolyte condition was
    stated in the narrative but not the parameter name -- this question
    states "measured in 100 mM NaCl" in prose, exactly as a real paper
    would, and checks whether the model (or the toolkit orchestrator,
    which requires an explicit electrolyte_condition argument and flags
    it as a gap if omitted) gets the prefactor right without being told
    the parameter's name.
  Q8: REAL raw robotic DTAB data (same source as Q2) deliberately
    TRUNCATED to 5 points (real subset of a real curve, not fabricated
    values) specifically so too few points survive the premicellar
    segment for isotherm model selection (needs >=4) -- tests whether
    the model correctly reports "cannot determine Langmuir vs Frumkin"
    rather than forcing a fit on insufficient data.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from surfactantkit.orchestrate import derive_all_properties_from_smiles_and_curve
from surfactantkit.adsorption import szyszkowski_surface_tension

HERE = Path(__file__).resolve().parent


@dataclass
class RoundQuestion:
    id: str
    smiles: str
    concentrations_mM: list[float]
    surface_tensions_mN_per_m: list[float]
    temperature_K: float
    narrative_extra: str  # any extra real-world context stated in prose (e.g. "measured in 100 mM NaCl")
    data_provenance: str  # honestly disclosed source, shown ONLY in the answer key, never pasted to the LLM
    gold: dict  # computed by actually calling the orchestrator


QUESTIONS_RAW = []

# Q1 -- SDS, real PendantProp robotic data (Dankloff et al. 2025)
QUESTIONS_RAW.append(dict(
    id="R3-01",
    smiles="CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]",
    concentrations_mM=[0.40577109375, 0.8115421875, 1.623084375, 3.24616875, 6.4923375, 12.984675, 25.96935, 51.9387],
    surface_tensions_mN_per_m=[69.85933451975876, 67.9370105943862, 62.892649575815454, 54.57301518905874,
                                42.87619733622359, 37.91402789831731, 37.44295519808575, 36.79102755846683],
    temperature_K=296.15,
    narrative_extra="",
    data_provenance="Dankloff et al. 2025 (npj Comput. Mater. 11:358, PendantProp), real robotic pendant-drop data, SDS replicate.",
))

# Q2 -- DTAB, real PendantProp
QUESTIONS_RAW.append(dict(
    id="R3-02",
    smiles="CCCCCCCCCCCC[N+](C)(C)C.[Br-]",
    concentrations_mM=[0.25578125, 0.5115625, 1.023125, 2.04625, 4.0925, 16.37, 32.74],
    surface_tensions_mN_per_m=[70.26907673, 68.20699668, 65.39702756, 57.50386199, 49.23642286, 32.34750118, 36.22382546],
    temperature_K=296.35,
    narrative_extra="",
    data_provenance="Dankloff et al. 2025 (PendantProp), real robotic pendant-drop data, DTAB replicate 7A.",
))

# Q3 -- CTAB, real PendantProp (known real poor 2-param fit)
QUESTIONS_RAW.append(dict(
    id="R3-03",
    smiles="CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]",
    concentrations_mM=[0.016484375, 0.03296875, 0.0659375, 0.131875, 0.26375, 0.5275, 1.055, 2.11],
    surface_tensions_mN_per_m=[71.60374462, 71.55862868, 70.38023925, 64.77683517, 57.3369305, 46.13427154, 36.41066498, 32.37490721],
    temperature_K=296.85,
    narrative_extra="",
    data_provenance="Dankloff et al. 2025 (PendantProp), real robotic pendant-drop data, CTAB replicate 7A.",
))

# Q4 -- AOT, real Shah/Das/Bhattarai 2025 literature dataset (already the
# regression case behind test_curve_analysis.py's AOT literature test).
# AOT = docusate sodium = sodium bis(2-ethylhexyl) sulfosuccinate.
QUESTIONS_RAW.append(dict(
    id="R3-04",
    smiles="CCCCC(CC)COC(=O)CC(S(=O)(=O)[O-])C(=O)OCC(CC)CCCC.[Na+]",
    concentrations_mM=[0.01585, 0.03981, 0.10000, 0.25119, 0.63096, 1.58489, 3.981, 6.31, 10.0],
    surface_tensions_mN_per_m=[62.22, 55.66, 49.10, 42.54, 35.98, 29.42, 27.5, 27.8, 28.0],
    temperature_K=298.15,
    narrative_extra="",
    data_provenance="Shah, Das & Bhattarai 2025 (Heliyon, PMC11835642), AOT in water, real literature raw tensiometry.",
))

# Q5 -- nonionic, C12E8-like, ILLUSTRATIVE curve (disclosed, see module docstring)
_gamma0, _gmax5, _K5, _T5 = 72.0, 3.2e-6, 9.0, 298.15
_concs5 = [0.005, 0.01, 0.02, 0.04, 0.07, 0.11, 0.16, 0.22]
_gammas5 = [szyszkowski_surface_tension(c, _gamma0, _gmax5, _K5, "nonionic", _T5) for c in _concs5] + [30.5, 30.5, 30.5]
QUESTIONS_RAW.append(dict(
    id="R3-05",
    smiles="CCCCCCCCCCCCOCCOCCOCCOCCOCCOCCOCCO",  # C12E8
    concentrations_mM=_concs5 + [0.4, 0.8, 1.6],
    surface_tensions_mN_per_m=_gammas5,
    temperature_K=_T5,
    narrative_extra="",
    data_provenance="ILLUSTRATIVE (disclosed): generated via this project's own szyszkowski_surface_tension with "
                     "representative (not literature-transcribed) K=9.0, Gamma_max=3.2e-6 for a C12E8-type "
                     "nonionic; see module docstring for why the project's real nonionic Frumkin sources "
                     "(Zawala 2020 MIBC, Taylor/Valkovska/Bain 2003 C12E5/C12E3/C10E8) weren't used directly.",
))

# Q6 -- zwitterionic, CAPB-like, ILLUSTRATIVE curve (tests refusal, not numeric accuracy)
_gmax6, _K6, _T6 = 3.0e-6, 7.0, 298.15
_concs6 = [0.05, 0.1, 0.2, 0.4, 0.8, 1.3, 1.9]
_gammas6 = [szyszkowski_surface_tension(c, 72.0, _gmax6, _K6, "nonionic", _T6) for c in _concs6] + [33.0, 33.0]
QUESTIONS_RAW.append(dict(
    id="R3-06",
    smiles="CCCCCCCCCCCCC(=O)NCCC[N+](C)(C)CC(=O)[O-]",  # cocamidopropyl betaine
    concentrations_mM=_concs6 + [3.0, 5.0],
    surface_tensions_mN_per_m=_gammas6,
    temperature_K=_T6,
    narrative_extra="",
    data_provenance="ILLUSTRATIVE (disclosed): generated via szyszkowski_surface_tension with representative "
                     "parameters for a CAPB-type zwitterionic surfactant -- no real raw zwitterionic isotherm "
                     "curve has been mined in this project yet. Included to test refusal behavior "
                     "(Gamma_max/A_min/isotherm/deltaG_mic should all come back undetermined), which does not "
                     "depend on numeric literature accuracy.",
))

# Q7 -- SDS in 100 mM NaCl, REAL fitted K/Gamma_max (Prosser & Franses 2001,
# Table 1), ILLUSTRATIVE curve (raw points never published), electrolyte
# condition stated in prose only -- the exact v1 systematic-error test.
_gmax7, _K7, _T7 = 3.83e-6, 32.5, 298.15
_concs7 = [0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.4, 1.8]
_gammas7 = [szyszkowski_surface_tension(c, 72.0, _gmax7, _K7, "ionic_excess_electrolyte", _T7) for c in _concs7] + [33.21, 33.21, 33.21]
QUESTIONS_RAW.append(dict(
    id="R3-07",
    smiles="CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]",
    concentrations_mM=_concs7 + [2.5, 3.5, 5.0],
    surface_tensions_mN_per_m=_gammas7,
    temperature_K=_T7,
    narrative_extra="This tensiometry series was measured in an aqueous solution containing 100 mM added NaCl.",
    data_provenance="REAL fitted K=32.5 m^3/mol, Gamma_max=3.83e-6 mol/m^2 (Prosser & Franses 2001, Colloids "
                     "Surf. A 178:1-40, Table 1, SDS + 100 mM NaCl -- already validated in this project's own "
                     "szyszkowski_fit_K tests). Raw (C,gamma) pairs were never published in that paper; this "
                     "curve is an ILLUSTRATIVE reconstruction from those real fitted parameters, with a "
                     "disclosed-illustrative postmicellar plateau. Tests the EXACT systematic Gibbs-prefactor "
                     "error found in this project's own v1 unaugmented pilot (electrolyte condition stated in "
                     "prose, not as a labeled parameter).",
))

# Q8 -- DTAB, real data, deliberately truncated to 6 points (the minimum
# cmc_from_surface_tension_curve accepts at all) so too few survive the
# premicellar segment specifically for isotherm model selection (needs >=4).
QUESTIONS_RAW.append(dict(
    id="R3-08",
    smiles="CCCCCCCCCCCC[N+](C)(C)C.[Br-]",
    concentrations_mM=[0.25578125, 0.5115625, 1.023125, 2.04625, 4.0925, 16.37],
    surface_tensions_mN_per_m=[70.26907673, 68.20699668, 65.39702756, 57.50386199, 49.23642286, 32.34750118],
    temperature_K=296.35,
    narrative_extra="",
    data_provenance="Dankloff et al. 2025 (PendantProp), real robotic DTAB data (same source as R3-02), "
                     "deliberately truncated to 6 points -- a real subset, not fabricated values -- to test "
                     "whether too few premicellar points correctly blocks isotherm model selection rather than "
                     "forcing a fit.",
))


def build() -> list[RoundQuestion]:
    out = []
    for q in QUESTIONS_RAW:
        r = derive_all_properties_from_smiles_and_curve(
            q["smiles"], q["concentrations_mM"], q["surface_tensions_mN_per_m"], temperature_K=q["temperature_K"],
            electrolyte_condition="excess_electrolyte" if q["id"] == "R3-07" else None,
        )
        gold = {
            "charge_type": r.classification.charge_type,
            "cmc_mM": r.cmc.cmc_mM,
            "gamma_max_mol_per_m2": r.gamma_max_mol_per_m2,
            "a_min_nm2": r.a_min_nm2,
            "isotherm_model": r.isotherm.selected_model if r.isotherm else None,
            "frumkin_a": r.isotherm.frumkin_a if (r.isotherm and r.isotherm.selected_model == "frumkin") else None,
            "delta_g_mic_kJ_per_mol": r.delta_g_mic_kJ_per_mol,
            "gaps": r.gaps,
        }
        out.append(RoundQuestion(
            id=q["id"], smiles=q["smiles"], concentrations_mM=q["concentrations_mM"],
            surface_tensions_mN_per_m=q["surface_tensions_mN_per_m"], temperature_K=q["temperature_K"],
            narrative_extra=q["narrative_extra"], data_provenance=q["data_provenance"], gold=gold,
        ))
    return out


QUESTION_TEMPLATE = """R3-{n:02d}
[reference only, do not paste -- data provenance: {provenance}]

>>> PASTE BELOW >>>
You are analyzing a real surfactant using the data below. Work entirely from first principles -- do not search the internet, do not assume any specific isotherm model, do not assume ionic vs nonionic character beyond what the SMILES structure itself tells you, and do not use any pre-known literature CMC for this compound.

SMILES: {smiles}
Temperature: {temperature_K} K
Raw tensiometry data (concentration in mM, surface tension in mN/m):
concentration_mM = {concs}
surface_tension_mN_per_m = {gammas}
{narrative}
Based ONLY on the SMILES and this raw data, determine everything about this surfactant's micellization behavior that can genuinely be derived from these inputs alone. Explicitly state anything you conclude CANNOT be determined from this SMILES and data alone, and explain why -- do not guess or assume a value for anything not actually supported by the data given.

End your response with a block in exactly this form (use null for anything not determinable, and list every real thing you could not determine in "undeterminable" with a short reason for each):

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


def render_questions_txt(questions: list[RoundQuestion]) -> str:
    blocks = []
    for i, q in enumerate(questions, 1):
        narrative = (q.narrative_extra + "\n") if q.narrative_extra else ""
        blocks.append(QUESTION_TEMPLATE.format(
            n=i, provenance=q.data_provenance, smiles=q.smiles, temperature_K=q.temperature_K,
            concs=q.concentrations_mM, gammas=[round(g, 4) for g in q.surface_tensions_mN_per_m], narrative=narrative,
        ))
    return "\n\n".join(blocks)


def main():
    questions = build()

    (HERE / "round3_questions.json").write_text(
        json.dumps([asdict(q) for q in questions], indent=2), encoding="utf-8"
    )
    (HERE / "Round3_Questions_Unaugmented.txt").write_text(render_questions_txt(questions), encoding="utf-8")

    print(f"Built {len(questions)} round-3 questions.")
    for q in questions:
        print(f"  {q.id}: charge={q.gold['charge_type']:10} cmc={q.gold['cmc_mM']}"
              f" isotherm={q.gold['isotherm_model']} deltaG={q.gold['delta_g_mic_kJ_per_mol']} "
              f"n_gaps={len(q.gold['gaps'])}")


if __name__ == "__main__":
    main()
