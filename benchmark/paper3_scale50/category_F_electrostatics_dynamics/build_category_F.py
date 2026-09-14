"""Category F -- electrostatics and dynamics. Different real tool family
again: Debye screening, Henry-equation zeta potential, Stokes-Einstein
hydrodynamic radius. No curve, no counterion-factor convention question --
straightforward real-formula application, but on real literature data
including a genuine NEGATIVE control most benchmark question sets would
never think to include.

F-01/F-02: debye_length at two real ionic strengths, cross-checked against
two independent sources (a textbook reference point at I=0.1M, and a real
2026 paper's own explicit formula statement at I=0.15M, physiological
saline) -- both already validated in this project's own literature notes.

F-03: zeta_potential_henry, a real (mobility, zeta) PAIR from Serafini et
al. 2019's Supporting Information (mixed Triton X-100/DTAB micelles) --
genuinely hard-to-find raw data (most papers report only zeta, not the
underlying mobility), Smoluchowski regime confirmed (not assumed) by
back-solving f(kappa*a) from all 7 of the paper's own points and finding
them clustered tightly near the Smoluchowski limit.

F-04: hydrodynamic_radius_stokes_einstein, Milone et al.'s DOSY-NMR case
-- a genuine, real NEGATIVE CONTROL. Feeding their real D into
Stokes-Einstein correctly gives ~1.57 nm, nowhere near their OWN
separately-reported DLS hydrodynamic radius (~40 nm) for the same system
-- not a bug, a real physical-chemistry fact (DOSY-NMR on a surfactant
system near/above its aggregation point typically reports a fast-exchange
population-averaged D dominated by free monomer, not the large-aggregate D
DLS measures). Tests whether the model applies the formula correctly to
the D it's actually given, rather than being swayed toward "correcting"
its answer to match a different, unrelated technique's number it might
recall or expect.

F-05: hydrodynamic_radius_stokes_einstein, Sutherland et al. 2009's real
Taylor-dispersion SDS micelle diffusion data (mean of 4 real measured D
values) -- the positive counterpart to F-04, landing in the real,
independently-established SDS micelle R_h range (~2.0-2.8 nm).

Usage:
  python build_category_F.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(r"H:\CodeProjects\SurfactantKit\src")))
from surfactantkit.electrostatics import debye_length, zeta_potential_henry
from surfactantkit.dynamics import hydrodynamic_radius_stokes_einstein

QUESTIONS = [
    {
        "id": "F-01", "kind": "debye",
        "source": "\"Revisiting Colloid Theory for Biomedicine,\" ACS Nano Medicine (2026), PMC13455036; cross-checked against an independent textbook reference point",
        "ionic_strength_M": 0.1, "temperature_K": 298.15,
        "gold": {"debye_length_nm": debye_length(0.1, 298.15)},
    },
    {
        "id": "F-02", "kind": "debye",
        "source": "\"Revisiting Colloid Theory for Biomedicine,\" ACS Nano Medicine (2026), PMC13455036 (their own explicitly stated formula and worked value for physiological saline)",
        "ionic_strength_M": 0.15, "temperature_K": 298.15,
        "gold": {"debye_length_nm": debye_length(0.15, 298.15)},
    },
    {
        "id": "F-03", "kind": "zeta",
        "source": "Serafini, Fernandez-Leyes, Sanchez, Pereyra, Schulz, Durand, Schulz & Ritacco, Colloids Surf. A (2019), Supporting Information Table SI-III, mixed Triton X-100/DTAB micelles",
        "electrophoretic_mobility_um_cm_per_Vs": 0.585, "viscosity_mPas": 0.89, "regime": "smoluchowski",
        "gold": {"zeta_potential_mV": zeta_potential_henry(0.585, 0.89, "smoluchowski")},
        "literature_cross_check": {"paper_zeta_mV": 7.47},
    },
    {
        "id": "F-04", "kind": "hydrodynamic_radius",
        "source": "Milone et al., Int. J. Mol. Sci., PMC12786710, CTAB/pillararene host-guest nanoparticle system (DOSY-NMR)",
        "diffusion_coefficient_cm2_per_s": 1.56e-6, "viscosity_mPas": 0.89, "temperature_K": 298.15,
        "gold": {"hydrodynamic_radius_nm": hydrodynamic_radius_stokes_einstein(1.56e-6, 0.89, 298.15)},
    },
    {
        "id": "F-05", "kind": "hydrodynamic_radius",
        "source": "Sutherland, Mercer, Everist & Leaist, J. Chem. Eng. Data 54(2) (2009) 272-278, Table 2, aqueous SDS (NaDS), Taylor dispersion of a trace solubilized decanol tracer -- mean of 4 real measured micelle diffusion coefficients across their tested concentration range",
        "diffusion_coefficient_cm2_per_s": sum([0.095e-5, 0.094e-5, 0.092e-5, 0.091e-5]) / 4,
        "viscosity_mPas": 0.8903, "temperature_K": 298.15,
        "gold": {"hydrodynamic_radius_nm": hydrodynamic_radius_stokes_einstein(sum([0.095e-5, 0.094e-5, 0.092e-5, 0.091e-5]) / 4, 0.8903, 298.15)},
    },
]

if __name__ == "__main__":
    (HERE / "category_F_questions.json").write_text(json.dumps(QUESTIONS, indent=2), encoding="utf-8")
    print(f"Built {len(QUESTIONS)} Category F questions.")
    for q in QUESTIONS:
        print(f"  {q['id']} ({q['kind']}): gold={q['gold']}")
