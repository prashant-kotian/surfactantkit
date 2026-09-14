"""Category C -- van't Hoff enthalpy, multi-temperature CMC series.

Directly re-tests the real historical trap this project already found and
fixed (2026-09-04): vant_hoff_enthalpy() had a genuine sign bug (dividing
by (1/T1-1/T2) instead of the correct (1/T2-1/T1)), giving the exact
opposite sign of the true enthalpy every time, undetected until an
unaugmented model's own from-scratch Gibbs-Helmholtz reasoning disagreed
with the tool's answer. Same real data source as SurfQSPR's own already-
built (uncommitted) thermodynamics_analysis.py module: Saito, Moroi,
Matuura 1980, J. Colloid Interface Sci. 76(1), 256-258, Table I --
n-alkylsulfonic acid CMC at 3 real temperatures (30/50/70C) for C12/C14/
C16, and 2 temperatures (50/70C) for C18. Gold computed here directly from
SurfactantKit's own real (now-correctly-signed) vant_hoff_enthalpy()
function -- cross-checked against SurfQSPR's independent reimplementation,
which gives identical numbers to 3 decimal places (a genuine, useful
double-implementation check, not just a paper cross-check).

Deliberately NOT using vant_hoff_multi_point_fit (the newer 3+-point
Gibbs-Helmholtz/deltaCp tool) here: its deltaG(Tj) term needs a real
counterion_factor, and no independently-measured counterion binding degree
for alkylsulfonic acids specifically was found in the literature searched
for this project -- using it would require guessing that parameter, the
exact "do not guess a system parameter" trap this project's own discipline
exists to avoid. Plain two-point vant_hoff_enthalpy needs no such factor
(it works directly on ln(X_cmc) vs. 1/T), so that's what's used and asked
for here, applied both to the endpoint temperature pair AND (where 3
temperatures exist) each adjacent sub-interval, as a real self-consistency
check on the "deltaH constant over this T range" assumption every two-
point estimate makes.

Real entangled variety:
  C-01 (C12): 3 temps, sub-intervals agree closely -- self-consistent case.
  C-02 (C14): 3 temps, sub-intervals disagree by ~90% of the endpoint value
       -- a REAL, disclosed self-consistency failure (not a data error);
       tests whether the model actually computes and compares sub-
       intervals rather than just reporting one plausible-looking endpoint
       number and calling it done.
  C-03 (C16): 3 temps, sub-intervals agree reasonably (moderate case,
       between C12's tight agreement and C14's clear failure).
  C-04 (C18): only 2 temperatures exist -- no sub-interval check is
       POSSIBLE at all. Tests whether the model correctly says so rather
       than fabricating a self-consistency claim it has no data to make.

Usage:
  python build_category_C.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(r"H:\CodeProjects\SurfactantKit\src")))
from surfactantkit import thermodynamics as thermo

SERIES = {
    12: {30.0: 7.8, 50.0: 10.2, 70.0: 13.0},
    14: {30.0: 2.38, 50.0: 2.80, 70.0: 4.05},
    16: {30.0: 0.592, 50.0: 0.805, 70.0: 1.11},
    18: {50.0: 0.226, 70.0: 0.317},
}
SMILES_BY_N = {
    12: "CCCCCCCCCCCCS(=O)(=O)O",
    14: "CCCCCCCCCCCCCCS(=O)(=O)O",
    16: "CCCCCCCCCCCCCCCCS(=O)(=O)O",
    18: "CCCCCCCCCCCCCCCCCCS(=O)(=O)O",
}
SOURCE = "Saito, Moroi, Matuura 1980, J. Colloid Interface Sci. 76(1), 256-258, Table I (real conductivity CMC, n-alkylsulfonic acid series)"
SELF_CONSISTENCY_THRESHOLD_PCT = 50.0  # |sub1-sub2| relative to |endpoint dH|, stated explicitly


def dh_kJ(cmc1_mM, t1_C, cmc2_mM, t2_C):
    x1 = thermo.cmc_to_mole_fraction(cmc1_mM / 1000.0)
    x2 = thermo.cmc_to_mole_fraction(cmc2_mM / 1000.0)
    return thermo.vant_hoff_enthalpy(x1, t1_C + 273.15, x2, t2_C + 273.15)


def build_question(qid, n_carbons):
    temps_dict = SERIES[n_carbons]
    ts = sorted(temps_dict)
    t_lo, t_hi = ts[0], ts[-1]
    dh_endpoint = dh_kJ(temps_dict[t_lo], t_lo, temps_dict[t_hi], t_hi)

    sub_intervals = None
    self_consistent = None
    if len(ts) >= 3:
        sub_intervals = []
        for i in range(len(ts) - 1):
            a, b = ts[i], ts[i + 1]
            sub_intervals.append({"t_lo_C": a, "t_hi_C": b, "delta_h_kJ_per_mol": dh_kJ(temps_dict[a], a, temps_dict[b], b)})
        spread = abs(sub_intervals[0]["delta_h_kJ_per_mol"] - sub_intervals[1]["delta_h_kJ_per_mol"])
        self_consistent = (100 * spread / abs(dh_endpoint)) < SELF_CONSISTENCY_THRESHOLD_PCT

    return {
        "id": qid, "n_carbons": n_carbons,
        "compound_name": f"1-{'dodecane' if n_carbons==12 else 'tetradecane' if n_carbons==14 else 'hexadecane' if n_carbons==16 else 'octadecane'}sulfonic acid (C{n_carbons})",
        "smiles": SMILES_BY_N[n_carbons],
        "source": SOURCE,
        "given_temps_C_cmc_mM": temps_dict,
        "gold": {
            "delta_h_endpoint_kJ_per_mol": dh_endpoint,
            "endpoint_temps_C": [t_lo, t_hi],
            "sub_intervals": sub_intervals,
            "self_consistent": self_consistent,
        },
    }


QUESTIONS = [build_question(f"C-{i:02d}", n) for i, n in enumerate([12, 14, 16, 18], 1)]

if __name__ == "__main__":
    (HERE / "category_C_questions.json").write_text(json.dumps(QUESTIONS, indent=2), encoding="utf-8")
    print(f"Built {len(QUESTIONS)} Category C questions.")
    for q in QUESTIONS:
        g = q["gold"]
        print(f"  {q['id']} (C{q['n_carbons']}): endpoint dH={g['delta_h_endpoint_kJ_per_mol']:.2f} kJ/mol, "
              f"self_consistent={g['self_consistent']}")
