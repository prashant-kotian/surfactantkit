"""Category B -- mixed-surfactant (Clint/Rubingh) entanglement.

Different real tool family than Category A (chain-length QSPR) and Round
3 (single-curve orchestrator): binary regular-solution theory. 8 real
published systems, reused directly from this project's own
literature_validation_notes.md / tests/test_mixed_micelle.py -- every
system there was already independently checked against a real paper
before this benchmark existed, so gold here is naturally double-sourced
(tool output computed fresh below, cross-referenced against each paper's
own reported value, with the real, sometimes-imperfect agreement
disclosed exactly as literature_validation_notes.md already documents it
-- never silently rounded to "matches").

Real entangled variety across the 8 (not the same question shape repeated):
  B-01 Clint-only + sign classification (DTAB-SDS, antagonistic case)
  B-02 x1+beta, NO Clint comparison available (Muherei TX-100/SDS -- the
       paper's own ideal-CMC column is a known, disclosed dead end, see
       literature_validation_notes.md; asking for it here would just
       reward guessing against an unverifiable number)
  B-03 full triad, exact/close match (Azum gemini G6/TX-114)
  B-04 x1+beta, Clint computable but NOT independently literature-checked
       for this system (Lee TTAB/Tween-20) -- tests whether the model
       (dis)honestly flags which of its own numbers are literature-cross-
       checked vs. tool-only
  B-05 CONVENTION-FLIP trap (McLachlan gemini/zwitterionic): the paper
       reports the mole fraction of component 2, not component 1 -- the
       question asks explicitly for component 1's, so a model that
       doesn't track which component the tool's x1 actually refers to
       will silently report the wrong one
  B-06 full triad, real biosurfactant system (Liu TX-100/rhamnolipid)
  B-07 full triad, the strongest literature control case (cholate/SDS --
       beta matches the paper almost exactly, not just in sign)
  B-08 ROOT-FINDING-ONLY trap (Bales DHPC/SDS, the real disclosed table
       outlier point): the model is given raw (alpha1, cmc_mix, pure
       CMCs) and must SOLVE for x1 itself -- no literature x1 value is
       given to copy. Gold is the tool's real solved root (~0.41), which
       this project's own validation work concluded is more likely
       correct than the source paper's own printed table value (0.29,
       likely a transcription error -- it breaks an otherwise smooth,
       monotonic trend across the paper's 5 neighboring compositions).
       Tests genuine equation-solving, not literature recall.

Usage:
  python build_category_B.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(r"H:\CodeProjects\SurfactantKit\src")))
from surfactantkit.mixed_micelle import clint_ideal_cmc, solve_rubingh_x, rubingh_beta


def _rub(alpha1, cmc_mix, cmc1, cmc2):
    x1 = solve_rubingh_x(alpha1, cmc_mix, cmc1, cmc2)
    beta = rubingh_beta(x1, alpha1, cmc_mix, cmc1) if x1 is not None else None
    return x1, beta


QUESTIONS = []

# --- B-01: DTAB(1)-SDS(2), antagonistic case, PMC6554738 ---
cmc1, cmc2 = 14.80, 8.00
alpha1, cmc_mix = 0.75, 13.00
clint = clint_ideal_cmc(alpha1, cmc1, cmc2)
QUESTIONS.append({
    "id": "B-01", "system": "DTAB (1) - SDS (2), DTAB-rich", "source": "Rodriguez et al., PMC6554738",
    "component1": {"name": "dodecyltrimethylammonium bromide (DTAB)", "pure_cmc_mM": cmc1},
    "component2": {"name": "sodium dodecyl sulfate (SDS)", "pure_cmc_mM": cmc2},
    "alpha1": alpha1, "cmc_mix_mM": cmc_mix,
    "ask": ["clint_ideal_cmc_mM", "synergy_classification"],
    "gold": {"clint_ideal_cmc_mM": clint, "rubingh_x1": None, "rubingh_beta": None,
             "synergy_classification": "antagonistic"},
    "literature_cross_check": {"clint_ideal_cmc_mM": {"paper": 12.21, "tool": clint, "pct_diff": 100*abs(clint-12.21)/12.21}},
    "prompt_note": None,
    "design_notes": "Real mixed CMC (13.00 mM) exceeds the Clint ideal (12.21 mM), so it classifies as antagonistic (real mixed CMC is HIGHER than the Clint ideal prediction here) -- classify synergy from that comparison directly, not from a Rubingh beta (not requested here). NOT shown to the model -- classifying synergy from its own computed Clint ideal vs. the given real mixed CMC is exactly the task; stating the comparison outcome in the prompt would hand over the answer. Gold field itself must stay the bare classification word ('antagonistic'), not this parenthetical explanation -- a grader doing exact-match categorical comparison (as grade_all.py does) would otherwise fail every model that correctly answers just 'antagonistic', which is exactly what happened until this was caught (2026-09-14, all 4 models' real transcripts independently converged on the same fix-revealing symptom).",
})

# --- B-02: Muherei TX-100(1)-SDS(2) ---
cmc1, cmc2 = 0.387, 3.468
alpha1, cmc_mix = 0.47, 0.547
x1, beta = _rub(alpha1, cmc_mix, cmc1, cmc2)
QUESTIONS.append({
    "id": "B-02", "system": "Triton X-100 (1) - SDS (2)", "source": "Muherei & Junin 2009, Asian J. Appl. Sci. 2(2), 115-127",
    "component1": {"name": "Triton X-100", "pure_cmc_mM": cmc1},
    "component2": {"name": "sodium dodecyl sulfate (SDS)", "pure_cmc_mM": cmc2},
    "alpha1": alpha1, "cmc_mix_mM": cmc_mix,
    "ask": ["rubingh_x1", "rubingh_beta", "synergy_classification"],
    "gold": {"clint_ideal_cmc_mM": None, "rubingh_x1": x1, "rubingh_beta": beta,
             "synergy_classification": "synergistic" if beta < 0 else "antagonistic"},
    "literature_cross_check": {"rubingh_x1": {"paper": 0.7501, "tool": x1, "pct_diff": 100*abs(x1-0.7501)/0.7501},
                                "rubingh_beta": {"paper": -1.888, "tool": beta, "pct_diff": 100*abs(beta-(-1.888))/abs(-1.888)}},
    "prompt_note": "Note: a Clint ideal CMC is not requested for this question -- only x1, beta, and synergy classification.",
    "design_notes": "Clint ideal CMC deliberately excluded from 'ask': this project's own literature validation found the paper's own stated 'ideal CMC' column does not reproduce from its own stated pure-component CMCs via the standard Clint formula (a real, disclosed dead end in the source paper itself, not a tool gap) -- asking for it here would only reward guessing against an already-known-unreliable number. The prompt_note only scopes the task (don't compute it), without explaining why, so as not to hint at a data-quality issue the model should not be assumed to know about.",
})

# --- B-03: Azum gemini G6(1)-TX114(2) ---
cmc1, cmc2 = 0.041, 0.263
alpha1, cmc_mix = 0.48, 0.061
clint = clint_ideal_cmc(alpha1, cmc1, cmc2)
x1, beta = _rub(alpha1, cmc_mix, cmc1, cmc2)
QUESTIONS.append({
    "id": "B-03", "system": "Gemini G6 (1) - Triton X-114 (2)", "source": "Azum et al. 2022, Biointerface Res. Appl. Chem. 12(6), 7416-7428",
    "component1": {"name": "gemini surfactant G6 (16-6-16)", "pure_cmc_mM": cmc1},
    "component2": {"name": "Triton X-114", "pure_cmc_mM": cmc2},
    "alpha1": alpha1, "cmc_mix_mM": cmc_mix,
    "ask": ["clint_ideal_cmc_mM", "rubingh_x1", "rubingh_beta", "synergy_classification"],
    "gold": {"clint_ideal_cmc_mM": clint, "rubingh_x1": x1, "rubingh_beta": beta, "synergy_classification": "synergistic"},
    "literature_cross_check": {
        "clint_ideal_cmc_mM": {"paper": 0.073, "tool": clint, "pct_diff": 100*abs(clint-0.073)/0.073},
        "rubingh_x1": {"paper": 0.764, "tool": x1, "pct_diff": 100*abs(x1-0.764)/0.764},
        "rubingh_beta": {"paper": -1.211, "tool": beta, "pct_diff": 100*abs(beta-(-1.211))/abs(-1.211)},
    },
    "prompt_note": None,
    "design_notes": "Strong agreement case: Clint ideal exact match, x1 close match, beta correct sign (~4% magnitude difference from the paper's own regression-fit value -- a real, expected pointwise-vs-regression gap, not an error).",
})

# --- B-04: Lee TTAB(1)-Tween20(2) ---
cmc1, cmc2 = 2.25, 0.41
alpha1, cmc_mix = 0.6, 0.70
clint = clint_ideal_cmc(alpha1, cmc1, cmc2)
x1, beta = _rub(alpha1, cmc_mix, cmc1, cmc2)
QUESTIONS.append({
    "id": "B-04", "system": "TTAB (1) - Tween-20 (2)", "source": "Lee & Lee 2012, J. Korean Chem. Soc. 56(5), 556-562",
    "component1": {"name": "tetradecyltrimethylammonium bromide (TTAB)", "pure_cmc_mM": cmc1},
    "component2": {"name": "Tween-20 (polysorbate 20)", "pure_cmc_mM": cmc2},
    "alpha1": alpha1, "cmc_mix_mM": cmc_mix,
    "ask": ["clint_ideal_cmc_mM", "rubingh_x1", "rubingh_beta", "synergy_classification"],
    "gold": {"clint_ideal_cmc_mM": clint, "rubingh_x1": x1, "rubingh_beta": beta, "synergy_classification": "synergistic"},
    "literature_cross_check": {"rubingh_x1": {"paper": 0.28, "tool": x1, "pct_diff": 100*abs(x1-0.28)/0.28}},
    "prompt_note": None,
    "design_notes": "The source paper does not report a Clint ideal CMC at all -- clint_ideal_cmc_mM is real, computable, and included in gold as a tool-only value (no independent literature check exists for it in THIS system, unlike x1/beta which do). Originally also asked the model to self-report which fields were literature-checkable, but dropped that sub-ask 2026-09-14: it isn't fairly answerable from the given data alone (it's a question about what the SOURCE PAPER reports, not something derivable from the numbers given, and telling the model would itself be a leak).",
})

# --- B-05: McLachlan gemini 12-4-12(1)-ZW3-12(2), CONVENTION FLIP ---
cmc1, cmc2 = 1.10, 2.63
alpha1, cmc_mix = 0.5, 1.52
clint = clint_ideal_cmc(alpha1, cmc1, cmc2)
x1, beta = _rub(alpha1, cmc_mix, cmc1, cmc2)
QUESTIONS.append({
    "id": "B-05", "system": "gemini 12-4-12 (1) - zwitterionic ZW3-12 (2)", "source": "McLachlan et al. 2020, RSC Adv. 10(6), 3221-3232",
    "component1": {"name": "gemini 12-4-12", "pure_cmc_mM": cmc1},
    "component2": {"name": "zwitterionic ZW3-12", "pure_cmc_mM": cmc2},
    "alpha1": alpha1, "cmc_mix_mM": cmc_mix,
    "ask": ["clint_ideal_cmc_mM", "rubingh_x1_component1"],
    "gold": {"clint_ideal_cmc_mM": clint, "rubingh_x1": x1, "rubingh_beta": None, "synergy_classification": None},
    "literature_cross_check": {
        "clint_ideal_cmc_mM": {"paper": 1.55, "tool": clint, "pct_diff": 100*abs(clint-1.55)/1.55},
        "rubingh_x1_component1": {"paper_reports_component2_x": 0.301, "component1_x_derived": 1-0.301,
                                   "tool_component1_x1": x1, "pct_diff": 100*abs(x1-(1-0.301))/(1-0.301)},
    },
    "prompt_note": None,
    "design_notes": "CONVENTION TRAP: the source paper reports the MICELLAR MOLE FRACTION OF COMPONENT 2 (ZW3-12), 0.301 -- not component 1's. The question explicitly asks for component 1's (gemini 12-4-12) micellar mole fraction. Correct gold is x1 (component 1) = 1 - 0.301 = 0.699, matching the tool's own x1 output directly (the tool always solves for component 1 as passed). A model that reports 0.301 as if it were x1 for component 1 has made a real, checkable convention error. No hint given in the prompt -- the question is fully self-contained (alpha1, cmc_mix, both pure CMCs given), so no literature lookup is needed to answer correctly regardless; this mainly tests whether a model that DOES recall or half-recall this real paper gets tripped up by its component-2 convention, and separately (augmented condition) whether the tool's own x1 output gets correctly attributed to component 1 in the model's final report.",
})

# --- B-06: Liu TX100(1)-rhamnolipid(2) ---
cmc1, cmc2 = 0.309, 0.134
alpha1, cmc_mix = 0.888, 0.253
clint = clint_ideal_cmc(alpha1, cmc1, cmc2)
x1, beta = _rub(alpha1, cmc_mix, cmc1, cmc2)
QUESTIONS.append({
    "id": "B-06", "system": "Triton X-100 (1) - rhamnolipid biosurfactant (2)", "source": "Liu et al. 2020, Molecules 25(18), 4327, Table 9",
    "component1": {"name": "Triton X-100", "pure_cmc_mM": cmc1},
    "component2": {"name": "rhamnolipid biosurfactant", "pure_cmc_mM": cmc2},
    "alpha1": alpha1, "cmc_mix_mM": cmc_mix,
    "ask": ["clint_ideal_cmc_mM", "rubingh_x1", "rubingh_beta", "synergy_classification"],
    "gold": {"clint_ideal_cmc_mM": clint, "rubingh_x1": x1, "rubingh_beta": beta, "synergy_classification": "synergistic"},
    "literature_cross_check": {
        "clint_ideal_cmc_mM": {"paper": 0.270, "tool": clint, "pct_diff": 100*abs(clint-0.270)/0.270},
        "rubingh_x1": {"paper": 0.744, "tool": x1, "pct_diff": 100*abs(x1-0.744)/0.744},
        "rubingh_beta": {"paper": -0.379, "tool": beta, "pct_diff": 100*abs(beta-(-0.379))/abs(-0.379)},
    },
    "prompt_note": None,
    "design_notes": "Real synthetic/biosurfactant mixture -- same discipline applies to a less \"textbook\" real-world system.",
})

# --- B-07: Cholate(1)-SDS(2), strongest control case ---
cmc1, cmc2 = 11.50, 11.98
alpha1, cmc_mix = 0.5, 4.07
clint = clint_ideal_cmc(alpha1, cmc1, cmc2)
x1, beta = _rub(alpha1, cmc_mix, cmc1, cmc2)
QUESTIONS.append({
    "id": "B-07", "system": "sodium cholate (1) - SDS (2), 1:1", "source": "Kang, Bahadur et al., PMC4087020",
    "component1": {"name": "sodium cholate (NaCA)", "pure_cmc_mM": cmc1},
    "component2": {"name": "sodium dodecyl sulfate (SDS)", "pure_cmc_mM": cmc2},
    "alpha1": alpha1, "cmc_mix_mM": cmc_mix,
    "ask": ["clint_ideal_cmc_mM", "rubingh_x1", "rubingh_beta", "synergy_classification"],
    "gold": {"clint_ideal_cmc_mM": clint, "rubingh_x1": x1, "rubingh_beta": beta, "synergy_classification": "synergistic"},
    "literature_cross_check": {
        "clint_ideal_cmc_mM": {"paper": 11.74, "tool": clint, "pct_diff": 100*abs(clint-11.74)/11.74},
        "rubingh_x1": {"paper": 0.503, "tool": x1, "pct_diff": 100*abs(x1-0.503)/0.503},
        "rubingh_beta": {"paper": -4.23, "tool": beta, "pct_diff": 100*abs(beta-(-4.23))/abs(-4.23)},
    },
    "prompt_note": None,
    "design_notes": "The strongest single validation case in this project's own survey: a FIXED 1:1 mixture study (no multi-point regression possible for the paper's own beta), so beta itself matches almost exactly here (not just its sign) -- use this as the case where a model's answer SHOULD land very close on every field, a real control. No hint given in the prompt -- this is deliberately presented as an ordinary question, not flagged as an 'easy' one.",
})

# --- B-08: Bales DHPC(2)-SDS(1) real disclosed table-outlier point ---
cmc_sds, cmc_dhpc = 8.3, 1.8
alpha1, cmc_mix = 0.60, 1.4  # the real disclosed outlier composition
x1, beta = _rub(alpha1, cmc_mix, cmc_sds, cmc_dhpc)
QUESTIONS.append({
    "id": "B-08", "system": "SDS (1) - DHPC (2)", "source": "Vautier-Giongo, Bakshi, Singh, Ranganathan, Hajdu & Bales, J. Colloid Interface Sci. 282(1) (2005) 149-155",
    "component1": {"name": "sodium dodecyl sulfate (SDS)", "pure_cmc_mM": cmc_sds},
    "component2": {"name": "1,2-diheptanoyl-sn-glycero-3-phosphocholine (DHPC)", "pure_cmc_mM": cmc_dhpc},
    "alpha1": alpha1, "cmc_mix_mM": cmc_mix,
    "ask": ["rubingh_x1"],
    "gold": {"clint_ideal_cmc_mM": None, "rubingh_x1": x1, "rubingh_beta": None, "synergy_classification": None},
    "literature_cross_check": {
        "rubingh_x1": {"paper_table_value": 0.29, "tool": x1, "pct_diff_vs_paper_table": 100*abs(x1-0.29)/0.29,
                        "caveat": "This project's own literature validation work concluded the paper's own printed table value (0.29) is LIKELY A TRANSCRIPTION ERROR, not the true solved root -- it breaks an otherwise smooth, monotonically-decreasing x1 sequence across this paper's 5 neighboring bulk compositions (0.63, 0.56, 0.48, [0.29 printed / 0.41 solved], 0.24), while the tool's own solved value (~0.41) fits that trend. Gold here is the real, honestly-computed root, not the paper's own printed (probably erroneous) number."},
    },
    "prompt_note": None,
    "design_notes": "ROOT-FINDING TRAP: only raw (alpha1, cmc_mix, pure CMCs) are given -- no literature x1 value is supplied anywhere in the prompt for the model to copy. It must actually solve Rubingh's implicit equation. This also tests whether the model, if it happens to recall or look up this specific real paper's printed table value (0.29), correctly prefers its own genuine calculation over a recalled number that a real published table itself is likely to have gotten wrong. No hint about this given in the prompt -- coaching 'trust your own calculation over a recalled number' would itself bias the test.",
})

if __name__ == "__main__":
    (HERE / "category_B_questions.json").write_text(json.dumps(QUESTIONS, indent=2), encoding="utf-8")
    print(f"Built {len(QUESTIONS)} Category B questions.")
    for q in QUESTIONS:
        print(f"  {q['id']}: {q['system']} -- gold={q['gold']}")
