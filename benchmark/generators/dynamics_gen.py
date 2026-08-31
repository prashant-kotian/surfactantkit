"""Category F: dynamics. hydrodynamic_radius."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from surfactantkit.dynamics import hydrodynamic_radius_stokes_einstein
from schema import Question, next_id

# (D_cm2_per_s, viscosity_mPas, temperature_K, label)
CASES = [
    (1.0e-6, 0.89, 298.15, "small micelle, water, 25 C"),
    (5.0e-7, 0.89, 298.15, "larger micelle/aggregate, water, 25 C"),
    (2.0e-6, 0.89, 298.15, "small, fast-diffusing aggregate"),
    (1.0e-6, 1.5, 298.15, "more viscous medium (e.g. glycerol-water mix)"),
    (1.0e-6, 0.89, 310.15, "physiological temperature"),
    (3.0e-7, 1.0, 298.15, "large vesicle-scale aggregate"),
]


def gen() -> list[Question]:
    qs: list[Question] = []
    i = 0
    for D, eta, T, label in CASES * 3:
        i += 1
        rh = hydrodynamic_radius_stokes_einstein(D, eta, T)
        use_unit_trap = i % 3 == 0
        if use_unit_trap:
            D_wrong_unit = D * 1e-4  # what it would look like in m^2/s -- a plausible trap value
            text = (
                f"A DLS experiment reports a diffusion coefficient of {D_wrong_unit:.3e} m^2/s "
                f"({label}), solvent viscosity {eta} mPa.s, at {T} K. Compute the hydrodynamic "
                f"radius in nm via the Stokes-Einstein equation."
            )
            given = {"diffusion_coefficient_m2_per_s_as_reported": D_wrong_unit, "viscosity_mPas": eta, "temperature_K": T}
        else:
            text = (
                f"A DLS experiment reports a diffusion coefficient of {D:.2e} cm^2/s ({label}), "
                f"solvent viscosity {eta} mPa.s, at {T} K. Compute the hydrodynamic radius in nm "
                f"via the Stokes-Einstein equation."
            )
            given = {"diffusion_coefficient_cm2_per_s": D, "viscosity_mPas": eta, "temperature_K": T}
        qs.append(Question(
            id=next_id("F", i), category="F", subcategory="hydrodynamic_radius",
            tools_required=["hydrodynamic_radius"], difficulty="medium",
            trap_type="unit_trap" if use_unit_trap else "none",
            question_text=text, given_data=given,
            gold_answer=round(rh, 4), tolerance={"rel": 0.01}, source_note=label,
        ))
    return qs


if __name__ == "__main__":
    from collections import Counter
    questions = gen()
    print(f"Generated {len(questions)} questions for Category F (dynamics)")
    print(Counter(q.trap_type for q in questions))
