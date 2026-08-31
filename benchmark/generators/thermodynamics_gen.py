"""Category G: thermodynamics. counterion_binding_degree,
gibbs_free_energy_micellization, vant_hoff_enthalpy, entropy_of_micellization."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from surfactantkit.thermodynamics import (
    cmc_to_mole_fraction, counterion_binding_degree,
    gibbs_free_energy_micellization, vant_hoff_enthalpy, entropy_micellization,
)
from schema import Question, next_id

# (slope_below, slope_above, label)
BETA_CASES = [
    (2.0, 1.0, "50% counterion binding"),
    (3.5, 1.2, "typical ionic surfactant conductometric titration"),
    (1.8, 1.5, "weak counterion binding"),
    (4.0, 1.0, "strong counterion binding, e.g. long-chain gemini"),
]

# (cmc_M, temperature_K, counterion_factor, label)
DG_CASES = [
    (0.008, 298.15, 1.6, "ionic surfactant, SDS-like CMC scale, beta~0.4"),
    (0.0006, 298.15, 1.0, "nonionic surfactant, low CMC"),
    (0.014, 298.15, 1.5, "DTAB-like ionic surfactant"),
    (0.00004, 298.15, 1.7, "gemini cationic surfactant, very low CMC"),
]

# (cmc1_M, T1_K, cmc2_M, T2_K, label)
VANTHOFF_CASES = [
    (0.008, 293.15, 0.0075, 313.15, "CMC decreasing slightly with temperature"),
    (0.014, 288.15, 0.0155, 308.15, "CMC increasing with temperature"),
    (0.00004, 293.15, 0.000038, 313.15, "gemini surfactant, small CMC-T dependence"),
]


def gen() -> list[Question]:
    qs: list[Question] = []
    i = 0

    for s_below, s_above, label in BETA_CASES * 3:
        i += 1
        beta = counterion_binding_degree(s_below, s_above)
        qs.append(Question(
            id=next_id("G", i), category="G", subcategory="counterion_binding_degree",
            tools_required=["counterion_binding_degree"], difficulty="easy",
            question_text=(
                f"In a conductometric titration ({label}), the slope of conductivity vs "
                f"concentration is {s_below} below the CMC and {s_above} above the CMC. "
                f"Compute the degree of counterion binding, beta, using the slope-ratio method."
            ),
            given_data={"slope_below_cmc": s_below, "slope_above_cmc": s_above},
            gold_answer=round(beta, 4), tolerance={"rel": 0.01}, source_note=label,
        ))

    for cmc_M, T, factor, label in DG_CASES * 3:
        i += 1
        x_cmc = cmc_to_mole_fraction(cmc_M)
        dg = gibbs_free_energy_micellization(x_cmc, T, factor)
        qs.append(Question(
            id=next_id("G", i), category="G", subcategory="gibbs_free_energy_micellization",
            tools_required=["gibbs_free_energy_micellization"], difficulty="hard",
            trap_type="multi_tool_chain",
            question_text=(
                f"A surfactant has CMC = {cmc_M*1000:.3f} mM in water at {T} K ({label}). Using "
                f"a counterion factor of {factor} (i.e. 2-beta for this ionic surfactant), first "
                f"convert the CMC to mole fraction (using the standard ~55.5 mol/L water "
                f"molarity approximation), then compute the standard Gibbs free energy of "
                f"micellization, deltaG_mic, in kJ/mol."
            ),
            given_data={"cmc_M": cmc_M, "temperature_K": T, "counterion_factor": factor},
            gold_answer=round(dg, 3), tolerance={"rel": 0.02}, source_note=label,
        ))

    for cmc1, T1, cmc2, T2, label in VANTHOFF_CASES * 4:
        i += 1
        x1 = cmc_to_mole_fraction(cmc1)
        x2 = cmc_to_mole_fraction(cmc2)
        dh = vant_hoff_enthalpy(x1, T1, x2, T2)
        qs.append(Question(
            id=next_id("G", i), category="G", subcategory="vant_hoff_enthalpy",
            tools_required=["vant_hoff_enthalpy"], difficulty="hard", trap_type="multi_tool_chain",
            question_text=(
                f"A surfactant's CMC is {cmc1*1000:.4f} mM at {T1} K and {cmc2*1000:.4f} mM at "
                f"{T2} K ({label}). Using the van't Hoff relation, compute the enthalpy of "
                f"micellization, deltaH_mic, in kJ/mol."
            ),
            given_data={"cmc1_M": cmc1, "temperature1_K": T1, "cmc2_M": cmc2, "temperature2_K": T2},
            gold_answer=round(dh, 3), tolerance={"rel": 0.03}, source_note=label,
        ))

    # entropy triad-completion questions, chaining the two previous results.
    # Full cross-product of DG_CASES x VANTHOFF_CASES (4x3=12) rather than a
    # length-limited zip, which was silently capping this at 4 questions --
    # the thinnest tool in the first draft of the bank.
    import itertools
    for (cmc_M, T, factor, dg_label), (cmc1, T1, cmc2, T2, dh_label) in itertools.product(DG_CASES, VANTHOFF_CASES):
        i += 1
        x_cmc = cmc_to_mole_fraction(cmc_M)
        dg = gibbs_free_energy_micellization(x_cmc, T, factor)
        x1, x2 = cmc_to_mole_fraction(cmc1), cmc_to_mole_fraction(cmc2)
        dh = vant_hoff_enthalpy(x1, T1, x2, T2)
        ds = entropy_micellization(dg, dh, T)
        qs.append(Question(
            id=next_id("G", i), category="G", subcategory="entropy_of_micellization",
            tools_required=["gibbs_free_energy_micellization", "vant_hoff_enthalpy", "entropy_of_micellization"],
            difficulty="hard", trap_type="multi_tool_chain",
            question_text=(
                f"For a surfactant with deltaG_mic = {dg:.3f} kJ/mol and deltaH_mic = {dh:.3f} "
                f"kJ/mol at {T} K, compute the entropy of micellization, deltaS_mic, in J/(mol.K)."
            ),
            given_data={"deltaG_mic_kJ_per_mol": round(dg, 3), "deltaH_mic_kJ_per_mol": round(dh, 3), "temperature_K": T},
            gold_answer=round(ds, 3), tolerance={"rel": 0.02},
            source_note=f"{dg_label} / {dh_label}",
        ))

    return qs


if __name__ == "__main__":
    from collections import Counter
    questions = gen()
    print(f"Generated {len(questions)} questions for Category G (thermodynamics)")
    print(Counter(q.subcategory for q in questions))
