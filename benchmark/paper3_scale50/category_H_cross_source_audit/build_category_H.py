"""Category H -- cross-source audit: the SAME real compound (AOT), tested
under a completely different property/dataset than Round 3 already used
for it.

Round 3's R3-04 gave AOT's real raw tensiometry curve (Shah, Das &
Bhattarai 2025, Heliyon, PMC11835642) and asked for CMC/Gamma_max/
isotherm/deltaG_mic via the curve-orchestrator. This category gives the
SAME real compound's real wetting data (CA/adhesion tension, from the
SAME underlying paper's own Table 2, but a fully independent measurement
-- contact angle on glass, nothing to do with the tensiometry curve) and
asks for work of adhesion + spreading coefficient, exactly like Category
G. The point isn't a harder computation -- it's whether a full 42-
question benchmark run stays internally consistent when the same real
compound shows up more than once under a different property lens, which
is exactly the kind of thing a real research audit would check.

Same real back-derivation and verification discipline as Category G
(gamma = AT/cos(theta), verified self-consistent to 0.03%).

Usage:
  python build_category_H.py
"""
from __future__ import annotations
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(r"H:\CodeProjects\SurfactantKit\src")))
from surfactantkit.wetting import work_of_adhesion, spreading_coefficient

SOURCE = "Shah, Das & Bhattarai, Heliyon 2025, 11, e42XXX (open access, PMC11835642), Table 2 -- the SAME real paper whose raw tensiometry curve for this compound was already used in Round 3's R3-04 (a fully independent measurement: contact angle on glass, not tensiometry). gamma and theta back-derived from the paper's own cos(theta) and adhesion-tension columns, verified self-consistent (WA reproduces to 0.03%)."

ROWS = [
    ("H-01", "sodium bis(2-ethylhexyl) sulfosuccinate (AOT / Aerosol OT)", "water", 0.9109, 23.48),
    ("H-02", "sodium bis(2-ethylhexyl) sulfosuccinate (AOT / Aerosol OT)", "0.20 volume fraction 2-Propanol in water", 0.9376, 18.64),
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
        "cross_reference_note_INTERNAL_ONLY": "Same real compound as Round 3's R3-04 -- independent dataset (wetting, not tensiometry) from the same source paper. Not shown to the model.",
        "gold": {"work_of_adhesion_mJ_per_m2": wa, "spreading_coefficient_mJ_per_m2": s},
    })

if __name__ == "__main__":
    (HERE / "category_H_questions.json").write_text(json.dumps(QUESTIONS, indent=2), encoding="utf-8")
    print(f"Built {len(QUESTIONS)} Category H questions.")
    for q in QUESTIONS:
        print(f"  {q['id']}: gamma={q['gamma_LV_mN_per_m']} theta={q['contact_angle_deg']} gold={q['gold']}")
