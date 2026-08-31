"""Category E: electrostatics. debye_screening_length, zeta_potential."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from surfactantkit.electrostatics import ionic_strength, debye_length, zeta_potential_henry
from schema import Question, next_id

# (ionic_strength_M, temperature_K, label)
DEBYE_CASES = [
    (0.1, 298.15, "0.1 M 1:1 electrolyte (e.g. NaCl), 25 C"),
    (0.001, 298.15, "dilute, 1 mM ionic strength"),
    (0.01, 298.15, "10 mM ionic strength"),
    (0.5, 298.15, "concentrated brine-like ionic strength"),
    (0.1, 310.15, "0.1 M ionic strength at body temperature"),
    (0.05, 298.15, "50 mM ionic strength"),
]

# ionic-strength-from-concentrations cases (need to compute I first, a
# natural multi-step chain requiring z^2 weighting). Values are chosen
# so a naive "just sum the concentrations" shortcut gives a genuinely
# different (wrong) number than the correct charge-weighted formula --
# verified by direct computation, not assumed (a first draft using
# exact-stoichiometry CaCl2/Na2SO4 pairs accidentally made naive summing
# equal the correct answer for those specific ratios, which would have
# silently failed to test the thing it was meant to test).
IONIC_STRENGTH_CASES = [
    ({1: 0.1, -1: 0.1}, "0.1 M NaCl (1:1 salt)"),
    ({2: 0.1, -1: 0.15}, "Ca2+/Cl- mixed with excess background chloride"),
    ({1: 0.08, -2: 0.025}, "Na+/SO4(2-) mixed, non-stoichiometric background"),
]

# (mobility_um_cm_per_Vs, viscosity_mPas, regime, label)
ZETA_CASES = [
    (2.0, 0.89, "smoluchowski", "large colloidal particle, aqueous, 25 C"),
    (1.5, 0.89, "huckel", "small particle / low ionic strength"),
    (3.0, 1.0, "smoluchowski", "higher-mobility particle, slightly higher viscosity"),
    (0.8, 0.89, "huckel", "weakly charged small particle"),
]


def gen() -> list[Question]:
    qs: list[Question] = []
    i = 0

    for I, T, label in DEBYE_CASES * 2:
        i += 1
        ld = debye_length(I, T)
        use_unit_trap = i % 3 == 0
        text = (
            f"For an aqueous electrolyte solution with ionic strength I = {I} M at {T} K "
            f"({label}), compute the Debye screening length in nm."
        )
        qs.append(Question(
            id=next_id("E", i), category="E", subcategory="debye_screening_length",
            tools_required=["debye_screening_length"], difficulty="medium",
            question_text=text,
            given_data={"ionic_strength_M": I, "temperature_K": T},
            gold_answer=round(ld, 4), tolerance={"rel": 0.01}, source_note=label,
        ))

    for conc_dict, label in IONIC_STRENGTH_CASES * 3:
        i += 1
        I = ionic_strength(conc_dict)
        ld = debye_length(I, 298.15)
        conc_str = ", ".join(f"z={z}: {c} M" for z, c in conc_dict.items())
        qs.append(Question(
            id=next_id("E", i), category="E", subcategory="debye_screening_length",
            tools_required=["debye_screening_length"], difficulty="hard", trap_type="multi_tool_chain",
            question_text=(
                f"A solution contains the following ions ({label}): {conc_str}. First compute "
                f"the ionic strength I (mol/L), then the Debye screening length (nm) at 298.15 K."
            ),
            given_data={"ion_concentrations_M_by_charge": conc_dict, "temperature_K": 298.15},
            gold_answer={"ionic_strength_M": round(I, 4), "debye_length_nm": round(ld, 4)},
            tolerance={"rel": 0.01}, source_note=label,
        ))

    for mobility, eta, regime, label in ZETA_CASES * 3:
        i += 1
        zeta = zeta_potential_henry(mobility, eta, regime)
        qs.append(Question(
            id=next_id("E", i), category="E", subcategory="zeta_potential",
            tools_required=["zeta_potential"], difficulty="hard",
            question_text=(
                f"An electrophoretic mobility measurement gives {mobility} (um*cm)/(V*s) for a "
                f"system with solvent viscosity {eta} mPa.s ({label}). Given that this falls in "
                f"the {regime} regime, compute the zeta potential in mV via the Henry equation."
            ),
            given_data={"electrophoretic_mobility_um_cm_per_Vs": mobility, "viscosity_mPas": eta, "regime": regime},
            gold_answer=round(zeta, 3), tolerance={"rel": 0.01}, source_note=label,
        ))

    return qs


if __name__ == "__main__":
    from collections import Counter
    questions = gen()
    print(f"Generated {len(questions)} questions for Category E (electrostatics)")
    print(Counter(q.subcategory for q in questions))
