"""Category L: unit conversion & dimensional reasoning (30 questions).

Deliberately NOT tool-mapped -- these are simple enough that a calculation
tool shouldn't be needed, and that's the point: this category exists to show
tool augmentation is neutral here (per taxonomy.md), not to test whether the
model can multiply by a constant. Conversion factors are exact SI/CGS
identities or definitions, not values requiring literature lookup, so this
category is generated programmatically like Tier 1 (same "gold answer
computed, not hand-typed" discipline) even though it lives in Tier 2.

Sources for every conversion factor (verified, not recalled from memory
alone): NIST Guide to the SI (NIST SP 811) and standard CGS/SI definitions.
"""

from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from schema import Question, next_id

# --- exact conversion factors, cited ---
DYN_CM_PER_MN_M = 1.0          # mN/m and dyn/cm are the SAME unit (both = mJ/m^2 = 1e-3 N/m); identity, not approximate
PA_S_PER_CP = 1.0e-3           # 1 cP = 1 mPa.s = 1e-3 Pa.s (NIST SP 811)
NM_PER_ANGSTROM = 0.1          # 1 Angstrom = 1e-10 m = 0.1 nm, exact by definition
CM2_PER_M2 = 1.0e4             # 1 m^2 = 1e4 cm^2, exact
J_PER_CAL = 4.184              # thermochemical calorie, exact by definition (NIST)
WATER_DENSITY_KG_PER_L_25C = 0.997  # g/mL at 25C, standard reference value (used only for M<->mol/kg at dilute conditions)

i = 0
def qid():
    global i
    i += 1
    return next_id("L", i)

QUESTIONS: list[Question] = []

# 1. mN/m <-> dyn/cm (identity -- tests whether the model knows these are literally the same unit)
for val, direction in [(72.0, "mN/m to dyn/cm"), (35.5, "mN/m to dyn/cm"),
                       (28.3, "dyn/cm to mN/m"), (45.0, "dyn/cm to mN/m")]:
    src, dst = direction.split(" to ")
    QUESTIONS.append(Question(
        id=qid(), category="L", subcategory="surface_tension_unit_identity",
        tools_required=[], difficulty="easy",
        question_text=f"Convert {val} {src} to {dst}. (Note: think carefully about "
                       f"what these two units actually represent physically before converting.)",
        given_data={"value": val, "from_unit": src, "to_unit": dst},
        gold_answer=round(val * DYN_CM_PER_MN_M, 4),
        tolerance={"rel": 0.001}, grading_method="numeric_tolerance",
        source_note="mN/m and dyn/cm are the identical CGS/SI-derived unit for surface tension (both equal mJ/m^2); conversion factor is exactly 1, not approximately 1 -- NIST SP 811.",
    ))

# 2. cP <-> Pa.s (viscosity, exactly the kind of unit trap that appears throughout Tier 1)
for val, direction in [(0.89, "cP to Pa.s"), (1.5, "cP to Pa.s"), (10.0, "cP to Pa.s"),
                       (0.001, "Pa.s to cP"), (0.0015, "Pa.s to cP")]:
    src, dst = direction.split(" to ")
    factor = PA_S_PER_CP if src == "cP" else 1.0 / PA_S_PER_CP
    QUESTIONS.append(Question(
        id=qid(), category="L", subcategory="viscosity_unit_conversion",
        tools_required=[], difficulty="easy",
        question_text=f"Convert {val} {src} to {dst}.",
        given_data={"value": val, "from_unit": src, "to_unit": dst},
        gold_answer=round(val * factor, 6),
        tolerance={"rel": 0.001}, grading_method="numeric_tolerance",
        source_note="1 cP = 1 mPa.s = 1e-3 Pa.s exactly (NIST SP 811); water at 25C is ~0.89 cP.",
    ))

# 3. mM <-> M (simple order of magnitude, but a real, common Tier-1-style trap)
for val, direction in [(15.0, "mM to M"), (0.008, "M to mM"), (250.0, "mM to M"), (0.0006, "M to mM")]:
    src, dst = direction.split(" to ")
    factor = 1e-3 if src == "mM" else 1e3
    QUESTIONS.append(Question(
        id=qid(), category="L", subcategory="concentration_unit_conversion",
        tools_required=[], difficulty="easy",
        question_text=f"Convert {val} {src} to {dst}.",
        given_data={"value": val, "from_unit": src, "to_unit": dst},
        gold_answer=val * factor,
        tolerance={"rel": 0.001}, grading_method="numeric_tolerance",
        source_note="mM to M is a direct 1000x SI prefix conversion, exact by definition.",
    ))

# 4. Angstrom <-> nm (molecular length scale, relevant to CPP/Tanford-style quantities)
for val, direction in [(16.68, "Angstrom to nm"), (2.3, "nm to Angstrom"),
                       (350.2, "cubic Angstrom to nm^3")]:
    if "cubic" in direction:
        src, dst = "A^3", "nm^3"
        factor = NM_PER_ANGSTROM ** 3
    else:
        src, dst = direction.split(" to ")
        factor = NM_PER_ANGSTROM if src == "Angstrom" else 1.0 / NM_PER_ANGSTROM
    QUESTIONS.append(Question(
        id=qid(), category="L", subcategory="length_unit_conversion",
        tools_required=[], difficulty="easy",
        question_text=f"Convert {val} {src} to {dst}.",
        given_data={"value": val, "from_unit": src, "to_unit": dst},
        gold_answer=round(val * factor, 8),
        tolerance={"rel": 0.001}, grading_method="numeric_tolerance",
        source_note="1 Angstrom = 0.1 nm exactly, by SI definition (1 Angstrom = 1e-10 m).",
    ))

# 5. diffusion coefficient cm^2/s <-> m^2/s (the classic Stokes-Einstein unit trap, tested standalone here)
for val, direction in [(1.0e-6, "cm^2/s to m^2/s"), (5.0e-10, "m^2/s to cm^2/s"), (2.0e-6, "cm^2/s to m^2/s")]:
    src, dst = direction.split(" to ")
    factor = 1e-4 if src == "cm^2/s" else 1e4
    QUESTIONS.append(Question(
        id=qid(), category="L", subcategory="diffusion_coefficient_unit_conversion",
        tools_required=[], difficulty="medium",
        question_text=f"Convert a diffusion coefficient of {val} {src} to {dst}.",
        given_data={"value": val, "from_unit": src, "to_unit": dst},
        gold_answer=val * factor,
        tolerance={"rel": 0.001}, grading_method="numeric_tolerance",
        source_note="1 cm^2/s = 1e-4 m^2/s exactly (length^2 scaling of the 100x cm-to-m factor); this exact unit trap appears in Category F (hydrodynamic_radius) as a live tool-input trap.",
    ))

# 6. mol/m^2 <-> molecules/nm^2 (surface excess concentration, ties to Gibbs adsorption)
AVOGADRO = 6.02214076e23
for val in [2.04e-6, 3.94e-6, 1.71e-6]:
    QUESTIONS.append(Question(
        id=qid(), category="L", subcategory="surface_excess_unit_conversion",
        tools_required=[], difficulty="medium",
        question_text=f"A surface excess concentration is {val} mol/m^2. Express this as "
                       f"molecules per nm^2.",
        given_data={"value_mol_per_m2": val},
        gold_answer=round(val * AVOGADRO * 1e-18, 5),
        tolerance={"rel": 0.005}, grading_method="numeric_tolerance",
        source_note="mol/m^2 -> molecules/nm^2 requires multiplying by Avogadro's number then converting m^2 to nm^2 (1 m^2 = 1e18 nm^2) -- a two-step unit conversion, not a single factor, testing whether the model chains both steps correctly.",
    ))

# 7. degrees <-> radians (contact angle context, ties to Young-Dupre questions)
import math
for val, direction in [(41.3, "deg to rad"), (1.221, "rad to deg"), (90.0, "deg to rad")]:
    src, dst = direction.split(" to ")
    QUESTIONS.append(Question(
        id=qid(), category="L", subcategory="angle_unit_conversion",
        tools_required=[], difficulty="easy",
        question_text=f"Convert a contact angle of {val} {src} to {dst}.",
        given_data={"value": val, "from_unit": src, "to_unit": dst},
        gold_answer=round(math.radians(val) if src == "deg" else math.degrees(val), 5),
        tolerance={"rel": 0.001}, grading_method="numeric_tolerance",
        source_note="Standard degrees-radians conversion (pi/180), relevant since Young-Dupre/spreading-coefficient formulas require radians internally.",
    ))

# 8. mg/L (or g/L) <-> mM, using a known surfactant molecular weight
# (MW values are standard, widely-tabulated formula weights, not literature
# measurements -- SDS = C12H25SO4Na = 288.38 g/mol; CTAB = C19H42BrN = 364.45 g/mol)
for mw, mw_label, val, direction in [
    (288.38, "SDS", 2.31, "g/L to mM"),
    (288.38, "SDS", 8.0, "mM to g/L"),
    (364.45, "CTAB", 0.328, "g/L to mM"),
    (364.45, "CTAB", 1.0, "mM to g/L"),
    (650.0, "a rhamnolipid (MW 650 g/mol)", 130.0, "g/L to mM"),
]:
    src, dst = direction.split(" to ")
    if src == "g/L":
        gold = round((val / mw) * 1000.0, 4)  # g/L -> mol/L -> mM
    else:
        gold = round((val / 1000.0) * mw, 4)  # mM -> mol/L -> g/L
    QUESTIONS.append(Question(
        id=qid(), category="L", subcategory="mass_to_molar_concentration",
        tools_required=[], difficulty="medium",
        question_text=f"{mw_label} has a molecular weight of {mw} g/mol. Convert a concentration "
                       f"of {val} {src} to {dst}.",
        given_data={"value": val, "from_unit": src, "to_unit": dst, "molecular_weight_g_per_mol": mw},
        gold_answer=gold, tolerance={"rel": 0.01}, grading_method="numeric_tolerance",
        source_note=f"{mw_label} molecular weight is a standard tabulated formula weight, not a "
                    f"literature-measured value subject to interpretation -- SDS/CTAB formulas and "
                    f"MW are unambiguous (C12H25SO4Na=288.38, C19H42BrN=364.45).",
    ))

if __name__ == "__main__":
    print(f"Category L: {len(QUESTIONS)} questions")
    from collections import Counter
    print(Counter(q.subcategory for q in QUESTIONS))
