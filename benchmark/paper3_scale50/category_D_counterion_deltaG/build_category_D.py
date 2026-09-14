"""Category D -- counterion-binding-SUPPLIED deltaG_mic (the mirror-image
trap to Round 3's refusal traps).

Every Round 3 ionic-surfactant question had gold deltaG_mic=null because
no independent counterion-binding measurement (alpha/beta) was supplied --
correctly, since a surface-tension curve and SMILES alone can't provide
it. This category is the other half of that same judgment call: here a
REAL, independently-measured beta (from real conductivity data, published
alongside the CMC) IS given, so deltaG_mic SHOULD be computed, not
refused. A model that reflexively refuses every ionic-surfactant deltaG_mic
question (matching Round 3's correct pattern there) would be WRONG here --
this category exists specifically to check the refusal isn't just a fixed
reflex, but actually tracks whether the needed data is present.

Real double-sourced gold: SurfactantKit's own gibbs_free_energy_micellization
computed fresh from each paper's own real (CMC, T, beta) triple, using
EXACTLY that paper's own stated counterion_factor formula (never guessed,
always given explicitly in the prompt, matching the module's own "not
auto-detected, do not guess a convention" design) -- cross-checked against
each paper's own reported deltaG_mic (agreement 0.01-1.09%, both numbers
kept, not smoothed over).

Real entangled variety across the 5:
  D-01: Fu et al. 2019, [C12mpy][Br] pyridinium ionic liquid surfactant --
        the STANDARD ionic counterion_factor = (2-beta). Already the
        strongest single validation case in this project's own literature
        notes (0.03% agreement).
  D-02..D-05: Xie et al. 2018, dendritic cationic TETRAMERIC surfactants
        (4C12/14/16tetraQ) -- a DIFFERENT counterion_factor formula,
        (0.25+beta), specific to this 4-headgroup architecture (their own
        Eq. 5). Same underlying computation, different multiplicative
        prefactor -- tests whether the model correctly uses the SPECIFIC
        formula given for THIS question rather than defaulting to the
        more commonly-seen (2-beta) ionic convention out of habit. Two
        different chain lengths (C12, C16) and two temperatures for C16
        (25C, 45C) give real within-series and within-compound variety,
        not the same single data point four times.

Usage:
  python build_category_D.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(r"H:\CodeProjects\SurfactantKit\src")))
from surfactantkit import thermodynamics as thermo


def gold_dg(cmc_mM, temp_C, counterion_factor):
    x = thermo.cmc_to_mole_fraction(cmc_mM / 1000.0)
    return thermo.gibbs_free_energy_micellization(x, temp_C + 273.15, counterion_factor)


QUESTIONS = []

# D-01: Fu et al. 2019, standard ionic (2-beta)
cmc, T, beta = 9.74, 25.0, 0.268
cf = 2 - beta
dg = gold_dg(cmc, T, cf)
QUESTIONS.append({
    "id": "D-01",
    "compound_name": "1-dodecyl-3-methylpyridinium bromide, [C12mpy][Br] (pyridinium ionic-liquid surfactant)",
    "source": "Fu et al., RSC Advances 2019, 9, 28799-28807 (open access, PMC9071189), Tables 3-4",
    "architecture_class": "monomeric ionic (single headgroup)",
    "cmc_mM": cmc, "temperature_C": T, "beta": beta,
    "formula_given": "deltaG_mic = R*T*(2 - beta)*ln(X_cmc)  [standard ionic 1:1 surfactant convention]",
    "counterion_factor": cf,
    "gold": {"counterion_factor": cf, "delta_g_mic_kJ_per_mol": dg},
    "literature_cross_check": {"paper_deltaG_mic_kJ_per_mol": -37.12, "tool": dg, "pct_diff": 100*abs(dg-(-37.12))/37.12},
})

# D-02..D-05: Xie et al. 2018 tetrameric, (0.25+beta)
TETRA_ROWS = [
    ("D-02", "4C12tetraQ", 25.0, 0.412, 0.66, -26.64),
    ("D-03", "4C14tetraQ", 25.0, 0.0484, 0.69, -32.16),
    ("D-04", "4C16tetraQ", 25.0, 0.0063, 0.71, -38.05),
    ("D-05", "4C16tetraQ", 45.0, 0.0070, 0.67, -38.65),
]
for qid, name, T, cmc, beta, paper_dg in TETRA_ROWS:
    cf = 0.25 + beta
    dg = gold_dg(cmc, T, cf)
    QUESTIONS.append({
        "id": qid,
        "compound_name": f"{name} (dendritic cationic tetrameric quaternary-ammonium surfactant, 4 headgroups)",
        "source": "Xie, Li, Li, Sun, Wang & Qu, RSC Advances 2018, 8(63), 36015-36024 (open access), Table 2",
        "architecture_class": "dendritic tetrameric (4 headgroups)",
        "cmc_mM": cmc, "temperature_C": T, "beta": beta,
        "formula_given": "deltaG_mic = R*T*(0.25 + beta)*ln(X_cmc)  [this paper's own stated equation for its dendritic tetrameric architecture, Eq. 5]",
        "counterion_factor": cf,
        "gold": {"counterion_factor": cf, "delta_g_mic_kJ_per_mol": dg},
        "literature_cross_check": {"paper_deltaG_mic_kJ_per_mol": paper_dg, "tool": dg, "pct_diff": 100*abs(dg-paper_dg)/abs(paper_dg)},
    })

if __name__ == "__main__":
    (HERE / "category_D_questions.json").write_text(json.dumps(QUESTIONS, indent=2), encoding="utf-8")
    print(f"Built {len(QUESTIONS)} Category D questions.")
    for q in QUESTIONS:
        print(f"  {q['id']}: {q['compound_name'][:40]}... cf={q['gold']['counterion_factor']:.3f} "
              f"dG={q['gold']['delta_g_mic_kJ_per_mol']:.3f} vs paper={q['literature_cross_check']['paper_deltaG_mic_kJ_per_mol']} "
              f"({q['literature_cross_check']['pct_diff']:.2f}% diff)")
