"""Category H: wetting / EOR. wetting_work_of_adhesion,
wetting_spreading_coefficient, eor_capillary_number."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from surfactantkit.wetting import work_of_adhesion, spreading_coefficient, capillary_number
from schema import Question, next_id

# (gamma_LV_mN_m, contact_angle_deg, label)
WETTING_CASES = [
    (72.0, 41.3, "gemini surfactant-treated bentonite surface (Murtaza et al.-style)"),
    (72.0, 62.0, "EAPB + polyelectrolyte treated quartz (Gaynanova et al.-style)"),
    (30.0, 20.0, "low-surface-tension surfactant solution, good wetting"),
    (72.0, 90.0, "neutral wetting, untreated hydrophobic surface"),
    (45.0, 110.0, "poor wetting, hydrophobic-treated surface"),
    (72.0, 0.0, "complete wetting limit"),
]

# (viscosity_mPas, velocity_m_per_s, IFT_mN_m, label)
CA_CASES = [
    (1.0, 1e-5, 20.0, "waterflooding, normal oil-water IFT"),
    (1.0, 1e-5, 0.01, "surfactant flooding, ultra-low IFT"),
    (5.0, 1e-5, 0.001, "polymer-surfactant flood, ultra-low IFT, elevated viscosity"),
    (1.0, 1e-4, 20.0, "higher injection velocity, normal IFT"),
    (1.0, 1e-5, 3.4e-4, "extreme ultra-low IFT (Ahmad Wazir et al.-style, ~3.4e-4 mN/m)"),
]


def gen() -> list[Question]:
    qs: list[Question] = []
    i = 0

    for gamma, theta, label in WETTING_CASES * 2:
        i += 1
        wa = work_of_adhesion(gamma, theta)
        qs.append(Question(
            id=next_id("H", i), category="H", subcategory="wetting_work_of_adhesion",
            tools_required=["wetting_work_of_adhesion"], difficulty="easy",
            question_text=(
                f"A liquid with surface tension {gamma} mN/m forms a contact angle of {theta} "
                f"degrees on a solid surface ({label}). Compute the Young-Dupre work of "
                f"adhesion, in mJ/m^2."
            ),
            given_data={"gamma_LV_mN_m": gamma, "contact_angle_deg": theta},
            gold_answer=round(wa, 3), tolerance={"rel": 0.01}, source_note=label,
        ))

    for gamma, theta, label in WETTING_CASES * 2:
        i += 1
        s = spreading_coefficient(gamma, theta)
        qs.append(Question(
            id=next_id("H", i), category="H", subcategory="wetting_spreading_coefficient",
            tools_required=["wetting_spreading_coefficient"], difficulty="easy",
            question_text=(
                f"For the same system (surface tension {gamma} mN/m, contact angle {theta} "
                f"degrees, {label}), compute the spreading coefficient S in mN/m."
            ),
            given_data={"gamma_LV_mN_m": gamma, "contact_angle_deg": theta},
            gold_answer=round(s, 3), tolerance={"rel": 0.01}, source_note=label,
        ))

    for eta, v, ift, label in CA_CASES * 3:
        i += 1
        ca = capillary_number(eta, v, ift)
        qs.append(Question(
            id=next_id("H", i), category="H", subcategory="eor_capillary_number",
            tools_required=["eor_capillary_number"], difficulty="medium",
            question_text=(
                f"In an EOR waterflood, the displacing fluid has viscosity {eta} mPa.s and "
                f"velocity {v} m/s; the oil-water interfacial tension is {ift} mN/m ({label}). "
                f"Compute the capillary number and state whether flow is capillary-dominated "
                f"(Ca < 1e-5) or not."
            ),
            given_data={"viscosity_mPas": eta, "velocity_m_per_s": v, "interfacial_tension_mN_m": ift},
            gold_answer={"capillary_number": ca, "regime": "capillary-dominated" if ca < 1e-5 else "transitional/viscous-dominated"},
            tolerance={"rel": 0.02}, source_note=label,
        ))

    return qs


if __name__ == "__main__":
    from collections import Counter
    questions = gen()
    print(f"Generated {len(questions)} questions for Category H (wetting/EOR)")
    print(Counter(q.subcategory for q in questions))
