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

_SECTION_RE = re.compile(r"^#?\s*(R3-0[1-8])\b", re.MULTILINE)
# Anchored to LINE START + word boundary after the id (not full-line-only
# anymore, per the 2026-09-14 augmented-condition transcripts which have
# trailing text on the same header line, e.g.
# "R3-02  (SMILES: CCCCCCCCCCCCN+(C)C.[Br-], 296.35 K)"). Still safe against
# inline mid-paragraph mentions (e.g. "Contrast R3-01 (same SDS...)") since
# those never start a line.

# 2026-09-14: the SMILES pasted into Claude for the augmented run had its
# [N+] bracket notation stripped on these 4 questions (confirmed by directly
# comparing against Round3_Questions_Augmented.txt's real content -- not the
# same corruption in the ChatGPT transcript, which has all 8 SMILES intact).
# The tool correctly refused to parse the malformed SMILES on all 4 ("charge_
# type": "unparseable") rather than silently guessing -- a real, disciplined,
# correct response to a BAD INPUT, not a reasoning failure. Graded normally
# against gold below (so the aggregate score reflects reality), but flagged
# separately so this data-corruption artifact isn't silently blended into
# Claude's real augmented-condition performance without explanation.
CORRUPTED_INPUT = {("claude", "augmented", qid) for qid in ("R3-02", "R3-03", "R3-06", "R3-08")}


def split_transcript(text: str) -> dict[str, str]:
    """Split a transcript into {question_id: response_text} using the real
    header lines -- deliberately anchored to LINE START so inline mentions
    elsewhere in the prose (e.g. Claude's own "Contrast R3-01 (same SDS...)")
    are not mistaken for section boundaries."""
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


RUNS = [
    ("chatgpt", "unaugmented", HERE / "unagumented manual" / "R3-01_to_R3-08_answers chatgpt.txt"),
    ("claude", "unaugmented", HERE / "unagumented manual" / "R3_surfactant_analyses claude.txt"),
    ("chatgpt", "augmented", HERE / "aggumented manual" / "R3-01_to_R3-08_generated_outputs gpt.txt"),
    ("claude", "augmented", HERE / "aggumented manual" / "R3_answers claude.txt"),
]


def main():
    questions = {q["id"]: q for q in json.loads((HERE / "round3_questions.json").read_text(encoding="utf-8"))}

    results = {}     # (model,cond) -> qid -> field -> bool
    extracted = {}   # (model,cond) -> qid -> parsed FINAL_JSON dict (or None)

    for model, cond, path in RUNS:
        key = (model, cond)
        text = path.read_text(encoding="utf-8")
        sections = split_transcript(text)
        results[key] = {}
        extracted[key] = {}
        for qid, q in questions.items():
            resp = sections.get(qid)
            ans = extract_final(resp) if resp else None
            extracted[key][qid] = ans
            gold = q["gold"]
            field_results = {}
            for field in ALL_FIELDS:
                model_v = ans.get(field) if ans else None
                field_results[field] = grade_field(field, gold.get(field), model_v)
            results[key][qid] = field_results

    # ---- per-question x per-field table ----
    print("=" * 110)
    print("ROUND 3 GRADING -- per-question, per-field (P=pass, F=fail)")
    print(f"Numeric tolerance: {REL_TOLERANCE*100:.0f}% relative. Null-vs-value mismatch always fails "
          f"(no partial credit for a caveated-but-wrong FINAL_JSON value).")
    print("* = corrupted-SMILES input (see CORRUPTED_INPUT) -- graded normally but flag before trusting.")
    print("=" * 110)
    header = f"{'QID':6s} {'Model':8s} {'Cond':12s} " + " ".join(f"{f[:10]:10s}" for f in ALL_FIELDS) + "  Score"
    print(header)
    print("-" * len(header))
    totals = {k: 0 for k in results}
    max_total = {k: 0 for k in results}
    field_totals = {k: {f: 0 for f in ALL_FIELDS} for k in results}
    totals_clean = {k: 0 for k in results}       # excluding corrupted-input rows
    max_total_clean = {k: 0 for k in results}

    for qid in questions:
        for model, cond, _ in RUNS:
            key = (model, cond)
            fr = results[key][qid]
            corrupted = key + (qid,) in CORRUPTED_INPUT
            flag = "*" if corrupted else " "
            row = " ".join(f"{'PASS' if fr[f] else 'FAIL':10s}" for f in ALL_FIELDS)
            n_pass = sum(fr.values())
            totals[key] += n_pass
            max_total[key] += len(ALL_FIELDS)
            if not corrupted:
                totals_clean[key] += n_pass
                max_total_clean[key] += len(ALL_FIELDS)
            for f in ALL_FIELDS:
                field_totals[key][f] += int(fr[f])
            print(f"{qid:5s}{flag} {model:8s} {cond:12s} {row}  {n_pass}/{len(ALL_FIELDS)}")
        print()

    # ---- summary ----
    print("=" * 110)
    print("SUMMARY (literal grade, corrupted-input rows included)")
    print("=" * 110)
    for key in results:
        pct = 100 * totals[key] / max_total[key]
        print(f"{key[0]:10s} {key[1]:12s}: {totals[key]}/{max_total[key]} = {pct:.1f}%")

    print()
    print("SUMMARY excluding the 4 corrupted-SMILES rows (claude/augmented only affected)")
    print("-" * 80)
    for key in results:
        if max_total_clean[key] != max_total[key]:
            pct = 100 * totals_clean[key] / max_total_clean[key]
            print(f"{key[0]:10s} {key[1]:12s}: {totals_clean[key]}/{max_total_clean[key]} = {pct:.1f}% "
                  f"(literal, corrupted included: {100*totals[key]/max_total[key]:.1f}%)")

    print()
    print(f"{'Field':22s} " + " ".join(f"{k[0]}/{k[1][:4]:6s}" for k in results))
    for f in ALL_FIELDS:
        row = " ".join(f"{field_totals[k][f]}/8={100*field_totals[k][f]/8:5.1f}%" for k in results)
        print(f"{f:22s} {row}")

    out = {
        "rel_tolerance": REL_TOLERANCE,
        "corrupted_input_rows": sorted(f"{m}/{c}/{q}" for (m, c, q) in CORRUPTED_INPUT),
        "totals": {f"{k[0]}/{k[1]}": {"score": totals[k], "max": max_total[k], "pct": 100*totals[k]/max_total[k]} for k in results},
        "totals_excluding_corrupted": {
            f"{k[0]}/{k[1]}": {"score": totals_clean[k], "max": max_total_clean[k],
                                "pct": (100*totals_clean[k]/max_total_clean[k]) if max_total_clean[k] else None}
            for k in results
        },
        "field_totals": {f"{k[0]}/{k[1]}": {f: field_totals[k][f] for f in ALL_FIELDS} for k in results},
        "per_question": {f"{k[0]}/{k[1]}": {qid: results[k][qid] for qid in questions} for k in results},
        "extracted_answers": {f"{k[0]}/{k[1]}": extracted[k] for k in results},
    }
    out_path = HERE / "round3_grading_results.json"
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nFull results written to {out_path}")


if __name__ == "__main__":
    main()
