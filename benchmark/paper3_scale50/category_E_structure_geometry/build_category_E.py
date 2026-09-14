"""Category E -- structure/geometry, NO tensiometry curve at all.

Different real tool family again: Tanford packing-geometry (chain volume
-> critical length -> critical packing parameter -> predicted aggregate
morphology) and Griffin's HLB. No curve-fitting anywhere in this category
-- entanglement here is a real multi-step GEOMETRIC chain, not a curve
derivation.

E-01..E-04: real chained CPP worked example, this project's own hardest-
won literature validation case (Sachin, Karpe, Singh & Bhattarai 2019,
R. Soc. Open Sci. 6, 181979, PMC6458362, Table 2) -- SDS and DTAB are both
dodecyl (C12) chains, so Tanford volume/length are shared; only the real,
independently-MEASURED headgroup area (Amin, from the paper's own Gibbs-
adsorption surface-tension-slope measurement -- a real experimental input,
not assumed) differs. Chaining tanford_tail_volume -> tanford_critical_
length -> critical_packing_parameter -> classify_aggregate_morphology
reproduces the paper's own P to <1.2% and its own stated "cylindrical/
rod-shaped micelles" finding at every point. 4 real points across both
surfactants and multiple real temperatures (293.15-303.15K) -- not the
same single data point four times.

E-05..E-07: real Griffin HLB worked examples (Kosswig, "Surfactants,"
Ullmann's Encyclopedia of Industrial Chemistry Vol. 35, Table 13), 2
chemically different real nonionic hydrophobe classes (aliphatic
octadecanol vs. aromatic nonylphenol ethoxylates) at 3 real EO counts.
Gold carries an honest, real, systematic ~4-7% offset (Griffin's formula
assumes a monodisperse pure compound; real industrial ethoxylates are
polydisperse with some unreacted starting alcohol, per the same source
chapter) -- disclosed as a real, explained, EXPECTED gap, not something a
correct answer needs to somehow erase.

Usage:
  python build_category_E.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(r"H:\CodeProjects\SurfactantKit\src")))
from surfactantkit.cpp import tanford_tail_volume, tanford_critical_length, critical_packing_parameter, classify_aggregate_morphology
from surfactantkit.hlb import hlb_griffin

QUESTIONS = []

V12 = tanford_tail_volume(12)
LC12 = tanford_critical_length(12)
CPP_SOURCE = "Sachin, Karpe, Singh & Bhattarai, R. Soc. Open Sci. 6, 181979 (2019), open access, PMC6458362, Table 2"
CPP_ROWS = [
    ("E-01", "sodium dodecyl sulfate (SDS), SDS-rich mixed system", 44.70, 293.15, 0.47),
    ("E-02", "sodium dodecyl sulfate (SDS), SDS-rich mixed system", 49.37, 298.15, 0.43),
    ("E-03", "dodecyltrimethylammonium bromide (DTAB), DTAB-rich mixed system", 59.62, 298.15, 0.35),
    ("E-04", "dodecyltrimethylammonium bromide (DTAB), DTAB-rich mixed system", 61.59, 303.15, 0.34),
]
for qid, name, amin, T_K, p_paper in CPP_ROWS:
    p = critical_packing_parameter(V12, amin, LC12)
    morph = classify_aggregate_morphology(p)
    QUESTIONS.append({
        "id": qid, "kind": "cpp_chain",
        "compound_name": name, "source": CPP_SOURCE,
        "n_carbons": 12, "a_min_A2": amin, "temperature_K": T_K,
        "gold": {"tail_volume_A3": V12, "critical_length_A": LC12, "cpp": p, "predicted_morphology": morph},
        "literature_cross_check": {"paper_P": p_paper, "tool_P": p, "pct_diff": 100*abs(p-p_paper)/p_paper,
                                    "paper_morphology": "cylindrical or rod-shaped micelles"},
    })

HLB_SOURCE = "Kosswig, K., \"Surfactants,\" Ullmann's Encyclopedia of Industrial Chemistry Vol. 35 (Wiley-VCH, 2000), Table 13"
EO_MASS = 44.053
HLB_ROWS = [
    ("E-05", "octadecanol", 270.49, 12, 14.0),
    ("E-06", "isononylphenol", 220.35, 6, 11.5),
    ("E-07", "isononylphenol", 220.35, 5, 10.7),
]
for qid, hydrophobe, mw_hydrophobe, n_eo, table_hlb in HLB_ROWS:
    mh = n_eo * EO_MASS
    m_total = mw_hydrophobe + mh
    hlb = hlb_griffin(mh, m_total)
    QUESTIONS.append({
        "id": qid, "kind": "hlb_griffin",
        "compound_name": f"{hydrophobe} + {n_eo} ethylene oxide (EO) units (polyethoxylated nonionic surfactant)",
        "source": HLB_SOURCE,
        "hydrophobe_mw_g_per_mol": mw_hydrophobe, "n_eo_units": n_eo,
        "gold": {"mw_hydrophilic_g_per_mol": mh, "mw_total_g_per_mol": m_total, "hlb": hlb},
        "literature_cross_check": {"table_hlb": table_hlb, "tool_hlb": hlb, "pct_diff": 100*abs(hlb-table_hlb)/table_hlb},
    })

if __name__ == "__main__":
    (HERE / "category_E_questions.json").write_text(json.dumps(QUESTIONS, indent=2), encoding="utf-8")
    print(f"Built {len(QUESTIONS)} Category E questions.")
    for q in QUESTIONS:
        print(f"  {q['id']} ({q['kind']}): gold={q['gold']}")
