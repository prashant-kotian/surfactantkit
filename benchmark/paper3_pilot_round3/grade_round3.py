"""Real grading script for the Round 3 manual unaugmented pilot -- replaces
the hand-graded eyeball pass from 2026-09-14 with a locked-down, repeatable
tolerance rule, using this project's own real grading primitives
(harness/grading.py's _num_ok, harness/extract.py's extract_final) rather
than reimplementing numeric-tolerance logic from scratch.

Schema note: Round 3's gold answers are a 7-field structured dict per
question (charge_type, cmc_mM, gamma_max_mol_per_m2, a_min_nm2,
isotherm_model, frumkin_a, delta_g_mic_kJ_per_mol), not the single
gold_answer + tolerance shape harness/grading.py's grade() expects for
Tier 1/2 -- so this is a dedicated per-field grader, not a reuse of grade()
itself.

Per-field rule, stated explicitly (matches the 2026-09-14 hand-grading
exactly, now locked down instead of eyeballed):
  - charge_type, isotherm_model: exact string match (both null -> pass).
  - cmc_mM, gamma_max_mol_per_m2, a_min_nm2, frumkin_a, delta_g_mic_kJ_per_mol:
    if gold is null, PASS only if the answer is also null/None (a non-null
    answer where gold requires refusal is a fabrication, not partial credit).
    If gold is a real number, PASS only if the answer is also a real number
    within REL_TOLERANCE (default 20%) of gold (an answer of null where gold
    has a real number is a missed-determinable-value fail, not neutral).

Usage:
  python grade_round3.py
"""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "harness"))
from extract import extract_final
from grading import _as_float, _num_ok

REL_TOLERANCE = 0.20  # 20% relative tolerance for numeric fields, stated explicitly, not hidden

NUMERIC_FIELDS = ["cmc_mM", "gamma_max_mol_per_m2", "a_min_nm2", "frumkin_a", "delta_g_mic_kJ_per_mol"]
CATEGORICAL_FIELDS = ["charge_type", "isotherm_model"]
ALL_FIELDS = CATEGORICAL_FIELDS + NUMERIC_FIELDS

_SECTION_RE = re.compile(r"^#?\s*(R3-0[1-8])\s*$", re.MULTILINE)


def split_transcript(text: str) -> dict[str, str]:
    """Split a transcript into {question_id: response_text} using the real
    header lines (standalone 'R3-0N' or '# R3-0N' lines) -- deliberately
    anchored to a FULL LINE match so inline mentions elsewhere in the prose
    (e.g. Claude's own "Contrast R3-01 (same SDS...)") are not mistaken for
    section boundaries."""
    matches = list(_SECTION_RE.finditer(text))
    if not matches:
        raise ValueError(f"No 'R3-0N' section headers found in transcript -- check the format.")
    out = {}
    for i, m in enumerate(matches):
        qid = m.group(1)
        start = m.end()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        out[qid] = text[start:end].strip()
    return out


def grade_field(field: str, gold_v, model_v) -> bool:
    gold_is_null = gold_v is None
    model_is_null = model_v is None or (isinstance(model_v, str) and model_v.strip().lower() in ("null", "none", "n/a", ""))

    if field in CATEGORICAL_FIELDS:
        if gold_is_null:
            return model_is_null
        if model_is_null:
            return False
        return str(model_v).strip().lower() == str(gold_v).strip().lower()

    # numeric field
    if gold_is_null:
        return model_is_null
    if model_is_null:
        return False
    return _num_ok(model_v, float(gold_v), REL_TOLERANCE, None)


def main():
    questions = {q["id"]: q for q in json.loads((HERE / "round3_questions.json").read_text(encoding="utf-8"))}

    transcripts = {
        "chatgpt": HERE / "unagumented manual" / "R3-01_to_R3-08_answers chatgpt.txt",
        "claude": HERE / "unagumented manual" / "R3_surfactant_analyses claude.txt",
    }

    results = {}  # model -> qid -> field -> bool
    extracted = {}  # model -> qid -> parsed FINAL_JSON dict (or None)

    for model, path in transcripts.items():
        text = path.read_text(encoding="utf-8")
        sections = split_transcript(text)
        results[model] = {}
        extracted[model] = {}
        for qid, q in questions.items():
            resp = sections.get(qid)
            ans = extract_final(resp) if resp else None
            extracted[model][qid] = ans
            gold = q["gold"]
            field_results = {}
            for field in ALL_FIELDS:
                model_v = ans.get(field) if ans else None
                field_results[field] = grade_field(field, gold.get(field), model_v)
            results[model][qid] = field_results

    # ---- per-question x per-field table ----
    print("=" * 100)
    print("ROUND 3 GRADING -- per-question, per-field (P=pass, F=fail)")
    print(f"Numeric tolerance: {REL_TOLERANCE*100:.0f}% relative. Null-vs-value mismatch always fails "
          f"(no partial credit for a caveated-but-wrong FINAL_JSON value).")
    print("=" * 100)
    header = f"{'QID':6s} {'Model':8s} " + " ".join(f"{f[:10]:10s}" for f in ALL_FIELDS) + "  Score"
    print(header)
    print("-" * len(header))
    totals = {m: 0 for m in transcripts}
    max_total = {m: 0 for m in transcripts}
    field_totals = {m: {f: 0 for f in ALL_FIELDS} for m in transcripts}

    for qid in questions:
        for model in transcripts:
            fr = results[model][qid]
            row = " ".join(f"{'PASS' if fr[f] else 'FAIL':10s}" for f in ALL_FIELDS)
            n_pass = sum(fr.values())
            totals[model] += n_pass
            max_total[model] += len(ALL_FIELDS)
            for f in ALL_FIELDS:
                field_totals[model][f] += int(fr[f])
            print(f"{qid:6s} {model:8s} {row}  {n_pass}/{len(ALL_FIELDS)}")
        print()

    # ---- summary ----
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)
    for model in transcripts:
        pct = 100 * totals[model] / max_total[model]
        print(f"{model:10s}: {totals[model]}/{max_total[model]} = {pct:.1f}%")
    print()
    print(f"{'Field':22s} " + " ".join(f"{m:10s}" for m in transcripts))
    for f in ALL_FIELDS:
        row = " ".join(f"{field_totals[m][f]}/8={100*field_totals[m][f]/8:5.1f}%" for m in transcripts)
        print(f"{f:22s} {row}")

    out = {
        "rel_tolerance": REL_TOLERANCE,
        "totals": {m: {"score": totals[m], "max": max_total[m], "pct": 100*totals[m]/max_total[m]} for m in transcripts},
        "field_totals": {m: {f: field_totals[m][f] for f in ALL_FIELDS} for m in transcripts},
        "per_question": {m: {qid: results[m][qid] for qid in questions} for m in transcripts},
        "extracted_answers": extracted,
    }
    out_path = HERE / "round3_grading_results.json"
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nFull results written to {out_path}")


if __name__ == "__main__":
    main()
