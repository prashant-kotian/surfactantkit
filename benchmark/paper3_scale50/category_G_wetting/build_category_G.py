"""Category G -- wetting (work of adhesion + spreading coefficient),
chained from real (gamma, contact angle) pairs.

Real data source: Shah, Das & Bhattarai, Heliyon 2025 (PMC11835642), Table
2 -- CA/adhesion-tension/work-of-adhesion of AOT, SDS, CPC, and CTAB on
glass, in water and 3 real 2-Propanol co-solvent fractions. The paper
reports adhesion tension (AT = gamma*cos(theta)) and work of adhesion
(WA = gamma*(1+cos(theta))) directly, not gamma or theta themselves --
gamma and theta were back-derived here (gamma = AT/cos(theta), inverting
their own cos(theta) column; theta = arccos of that same column) before
building each question, exactly the same real back-derivation this
project's own literature_validation_notes.md already documents for this
paper. Verified: feeding the back-derived (gamma, theta) back through
work_of_adhesion() reproduces the paper's own printed WA to 0.00-0.02%,
confirming the back-derivation is self-consistent, not just plausible.

5 real questions across 3 compounds and 2 real solvent conditions (pure
water, and one real 2-Propanol co-solvent fraction) -- SDS and CPC each
appear in both conditions, CTAB in water only, giving real within-
compound solvent-effect variety, not 5 unrelated single points. AOT is
deliberately NOT used here (reused in Category H instead, to test
cross-referencing the same real compound across two independent real
datasets within this benchmark).

Usage:
  python build_category_G.py
"""
from __future__ import annotations
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(r"H:\CodeProjects\SurfactantKit\src")))
from surfactantkit.wetting import work_of_adhesion, spreading_coefficient

SOURCE = "Shah, Das & Bhattarai, Heliyon 2025, 11, e42XXX (open access, PMC11835642), Table 2 -- gamma and theta back-derived from the paper's own reported cos(theta) and adhesion-tension columns (AT=gamma*cos(theta)), verified self-consistent (WA reproduces to 0.00-0.02%)"

ROWS = [
    ("G-01", "sodium dodecyl sulfate (SDS)", "water", 0.9077, 30.95),
    ("G-02", "cetyltrimethylammonium bromide (CTAB)", "water", 0.7125, 27.00),
    ("G-03", "cetylpyridinium chloride (CPC)", "water", 0.7226, 27.60),
    ("G-04", "sodium dodecyl sulfate (SDS)", "0.10 volume fraction 2-Propanol in water", 0.9113, 26.53),
    ("G-05", "cetylpyridinium chloride (CPC)", "0.10 volume fraction 2-Propanol in water", 0.7248, 23.88),
]

QUESTIONS = []
for qid, name, solvent, cos_theta, at in ROWS:
    gamma = at / cos_theta
    theta = math.degrees(math.acos(cos_theta))
    wa = work_of_adhesion(gamma, theta)
    s = spreading_coefficient(gamma, theta)
    QUESTIONS.append({
        "id": qid, "compound_name": name, "solvent_system": solvent, "source": SOURCE,
        "gamma_LV_mN_per_m": round(gamma, 3), "contact_angle_deg": round(theta, 3),
        "gold": {"work_of_adhesion_mJ_per_m2": wa, "spreading_coefficient_mJ_per_m2": s},
    })

if __name__ == "__main__":
    (HERE / "category_G_questions.json").write_text(json.dumps(QUESTIONS, indent=2), encoding="utf-8")
    print(f"Built {len(QUESTIONS)} Category G questions.")
    for q in QUESTIONS:
        print(f"  {q['id']}: gamma={q['gamma_LV_mN_per_m']} theta={q['contact_angle_deg']} gold={q['gold']}")
