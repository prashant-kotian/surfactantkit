"""General grader for the Paper 3 50-question scale-up (Categories A-H).

Extends the real per-field tolerance rule already locked down in
paper3_pilot_round3/grade_round3.py (20% relative tolerance on numeric
fields, exact match on categorical fields, strict null-vs-value mismatch
-- no partial credit for a caveated-but-wrong FINAL_JSON value) to all 8
new categories, each with its own real gold schema (see each category's
own build_category_*.py for how gold was computed and double-sourced).

Designed to be run INCREMENTALLY as transcripts arrive -- categories/
models/conditions with no transcript file yet are skipped and reported as
"pending", not treated as failures. Point this at a transcripts folder
per (model, condition) containing one .txt file per category (or a
combined file per category -- the splitter finds every "<ID>" section
regardless of how many categories' worth of text are in one file), using
the same real header-splitting convention already proven across Round 3
and Categories A-H's own manual-run files: a real header line starting
with the question's own ID (e.g. "A-01", "B-05", "H-02"), optionally
prefixed with "#" or followed by more text on the same line.

Usage:
  python grade_all.py --transcripts-dir <dir with one .txt per (model,condition)>
  python grade_all.py  # uses the default layout under manual_runs/, see MODEL_RUNS below
"""
from __future__ import annotations
import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "harness"))
from extract import extract_final
from grading import _as_float, _num_ok

REL_TOLERANCE = 0.20  # same convention as grade_round3.py, stated explicitly

CATEGORIES = {
    "A": {"dir": "category_A_chain_length", "questions_file": "category_A_questions.json"},
    "B": {"dir": "category_B_mixed_surfactant", "questions_file": "category_B_questions.json"},
    "C": {"dir": "category_C_vanthoff_thermo", "questions_file": "category_C_questions.json"},
    "D": {"dir": "category_D_counterion_deltaG", "questions_file": "category_D_questions.json"},
    "E": {"dir": "category_E_structure_geometry", "questions_file": "category_E_questions.json"},
    "F": {"dir": "category_F_electrostatics_dynamics", "questions_file": "category_F_questions.json"},
    "G": {"dir": "category_G_wetting", "questions_file": "category_G_questions.json"},
    "H": {"dir": "category_H_cross_source_audit", "questions_file": "category_H_questions.json"},
}

# field name -> "num" (tolerance-checked) or "cat" (exact match, case-insensitive)
# or "skip" (reported/extracted but not graded -- process-disclosure fields,
# not a fixed-answer field).
FIELD_SPECS = {
    "A": {"predicted_cmc_mM": "num", "slope": "num", "intercept": "num", "r_squared": "num",
          "in_applicability_domain": "cat", "loo_cv_mean_pct_error": "num"},
    "B": {"clint_ideal_cmc_mM": "num", "rubingh_x1": "num", "rubingh_beta": "num",
          "synergy_classification": "cat"},
    "C": {"delta_h_endpoint_kJ_per_mol": "num", "self_consistent": "cat"},
    "D": {"counterion_factor": "num", "delta_g_mic_kJ_per_mol": "num"},
    "E_cpp_chain": {"tail_volume_A3": "num", "critical_length_A": "num", "cpp": "num",
                     "predicted_morphology": "cat"},
    "E_hlb_griffin": {"mw_hydrophilic_g_per_mol": "num", "mw_total_g_per_mol": "num", "hlb": "num"},
    "F_debye": {"debye_length_nm": "num"},
    "F_zeta": {"zeta_potential_mV": "num"},
    "F_hydrodynamic_radius": {"hydrodynamic_radius_nm": "num"},
    "G": {"work_of_adhesion_mJ_per_m2": "num", "spreading_coefficient_mJ_per_m2": "num"},
    "H": {"work_of_adhesion_mJ_per_m2": "num", "spreading_coefficient_mJ_per_m2": "num"},
}


def field_specs_for(category: str, question: dict) -> dict:
    key = category
    if category in ("E", "F"):
        key = f"{category}_{question.get('kind')}"
    return FIELD_SPECS[key]


_SECTION_RE = re.compile(r"^#?\s*([A-H]-\d{2})\b", re.MULTILINE)


def split_transcript(text: str) -> dict[str, str]:
    matches = list(_SECTION_RE.finditer(text))
    out = {}
    for i, m in enumerate(matches):
        qid = m.group(1)
        start = m.end()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        out[qid] = text[start:end].strip()
    return out


def grade_field(kind: str, gold_v, model_v) -> bool:
    gold_is_null = gold_v is None
    model_is_null = model_v is None or (isinstance(model_v, str) and model_v.strip().lower() in ("null", "none", "n/a", ""))
    if kind == "cat":
        if gold_is_null:
            return model_is_null
        if model_is_null:
            return False
        return str(model_v).strip().lower() == str(gold_v).strip().lower()
    # numeric
    if gold_is_null:
        return model_is_null
    if model_is_null:
        return False
    try:
        return _num_ok(model_v, float(gold_v), REL_TOLERANCE, None)
    except (TypeError, ValueError):
        return False


def load_all_questions() -> dict[str, dict]:
    """{question_id: question_dict}, question_dict augmented with 'category'."""
    out = {}
    for cat, spec in CATEGORIES.items():
        path = HERE / spec["dir"] / spec["questions_file"]
        qs = json.loads(path.read_text(encoding="utf-8"))
        for q in qs:
            q = dict(q, category=cat)
            out[q["id"]] = q
    return out


def grade_transcript(text: str, questions: dict) -> tuple[dict, dict]:
    """Returns (per_question_field_results, extracted_answers)."""
    sections = split_transcript(text)
    results, extracted = {}, {}
    for qid, q in questions.items():
        resp = sections.get(qid)
        ans = extract_final(resp) if resp else None
        extracted[qid] = ans
        specs = field_specs_for(q["category"], q)
        gold = q["gold"]
        field_results = {}
        for field, kind in specs.items():
            model_v = ans.get(field) if ans else None
            field_results[field] = grade_field(kind, gold.get(field), model_v)
        results[qid] = field_results
    return results, extracted


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--transcript", action="append", nargs=2, metavar=("LABEL", "PATH"),
                     help="a (label, path) pair, e.g. --transcript claude/unaugmented path/to/file.txt -- "
                          "repeat for every (model, condition) transcript you have so far")
    args = ap.parse_args()

    questions = load_all_questions()
    print(f"Loaded {len(questions)} real questions across {len(CATEGORIES)} categories.")

    if not args.transcript:
        print("No --transcript given. Usage:\n"
              "  python grade_all.py --transcript claude/unaugmented path/to/claude_unaug.txt "
              "--transcript gemini/augmented path/to/gemini_aug.txt ...")
        return

    all_results = {}
    for label, path in args.transcript:
        text = Path(path).read_text(encoding="utf-8")
        results, extracted = grade_transcript(text, questions)
        found = sum(1 for qid in questions if split_transcript(text).get(qid))
        totals = {"score": 0, "max": 0}
        cat_totals = {c: {"score": 0, "max": 0, "n_found": 0} for c in CATEGORIES}
        for qid, q in questions.items():
            fr = results[qid]
            n_pass = sum(fr.values())
            n_fields = len(fr)
            totals["score"] += n_pass
            totals["max"] += n_fields
            cat_totals[q["category"]]["score"] += n_pass
            cat_totals[q["category"]]["max"] += n_fields
            if extracted[qid] is not None:
                cat_totals[q["category"]]["n_found"] += 1

        print(f"\n=== {label} ({path}) ===")
        print(f"Questions with a parsed FINAL_JSON: {found}/{len(questions)}")
        pct = 100 * totals["score"] / totals["max"] if totals["max"] else 0.0
        print(f"Overall: {totals['score']}/{totals['max']} = {pct:.1f}%")
        for c in CATEGORIES:
            ct = cat_totals[c]
            cpct = 100 * ct["score"] / ct["max"] if ct["max"] else 0.0
            print(f"  Category {c}: {ct['score']}/{ct['max']} = {cpct:.1f}% ({ct['n_found']} questions found)")

        all_results[label] = {
            "totals": totals, "category_totals": cat_totals,
            "per_question": results, "extracted": extracted,
        }

    out_path = HERE / "grade_all_results.json"
    out_path.write_text(json.dumps(all_results, indent=2, default=str))
    print(f"\nFull results written to {out_path}")


if __name__ == "__main__":
    main()
