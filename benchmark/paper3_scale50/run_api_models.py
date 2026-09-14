"""Runs the real Gemini and Qwen flagship models against all 42 Category A-H
questions, both conditions (unaugmented / surfmcp-augmented), via the harness's
real provider adapters (benchmark/harness/providers.py) -- the same tool
dispatch (benchmark/harness/tools.py, 26 tools) that Claude Desktop and ChatGPT
Codex use through MCP, just called directly through each provider's native
function-calling API instead of through the MCP protocol.

Companion to the manual runs the user is doing by hand in the Claude and
ChatGPT chat apps (same 42 questions, same two conditions, pasted from the
CategoryX_Unaugmented.txt / CategoryX_Augmented.txt files in each category
folder) -- this script exists so all 4 models in Paper 3's comparison get
graded against the exact same question set under the exact same rule.

Reuses the real per-question prompt text already rendered into
CategoryX_Unaugmented.txt / CategoryX_Augmented.txt (the same text a human
pastes for Claude/GPT) -- extracted between the existing
">>> PASTE BELOW >>>" / "<<< PASTE ABOVE <<<" markers so the API models see
byte-identical task text to what the manual runs see. The
"[reference only, do not paste ...]" line above each marker (which can name
the held-out gold compound) is deliberately excluded from what's sent.

Output: one plain-text transcript per (model, condition), written to
api_runs/<model>_<condition>.txt, with "<ID>\n<response>\n\n" sections in the
exact format grade_all.py's splitter already expects -- so grading is just
  python grade_all.py --transcript gemini/augmented api_runs/gemini_augmented.txt ...
Also writes a parallel .json per (model, condition) with the full adapter
output (tool_calls, thinking, raw) for methodology/tool-usage analysis later.

RESUMABLE BY DESIGN: this machine has been observed killing long-running
background bash tasks unrelated to anything in this script (see memory:
feedback_background_task_reliability.md) -- a single ~45min "run everything,
write at the end" job loses all progress if that happens. Instead, each
question's real API response is written to its own file the moment it comes
back (api_runs/raw/<model>_<condition>/<qid>.json) -- a kill mid-run loses at
most the one in-flight question, and re-running the same command skips every
question whose raw file already exists rather than re-spending API calls on
it. The consolidated .txt/.json are rebuilt from whatever raw files exist
every time the script runs, so they are always gradable even mid-run.

Usage:
  python run_api_models.py                      # both models, both conditions, all 42
  python run_api_models.py --models gemini       # just one model
  python run_api_models.py --conditions unaugmented
  python run_api_models.py --categories A,B,C    # restrict to these categories (for
                                                  # chunking a run into pieces that fit
                                                  # a foreground command's time budget)
  python run_api_models.py --limit 2             # smoke test, 2 questions per category
  python run_api_models.py --force               # ignore existing raw files, re-run all
"""
from __future__ import annotations
import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
HARNESS = HERE.parent / "harness"
sys.path.insert(0, str(HARNESS))


def _load_env():
    envf = HARNESS / ".env"
    if not envf.exists():
        return
    for line in envf.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        os.environ.setdefault(k.strip(), v.strip().strip('"').strip("'"))


_load_env()

import providers
import prompts

OUT_DIR = HERE / "api_runs"

CATEGORY_FILES = {
    "A": "category_A_chain_length/CategoryA",
    "B": "category_B_mixed_surfactant/CategoryB",
    "C": "category_C_vanthoff_thermo/CategoryC",
    "D": "category_D_counterion_deltaG/CategoryD",
    "E": "category_E_structure_geometry/CategoryE",
    "F": "category_F_electrostatics_dynamics/CategoryF",
    "G": "category_G_wetting/CategoryG",
    "H": "category_H_cross_source_audit/CategoryH",
}

MODELS = {
    "gemini": {"provider": "gemini", "model_id": "gemini-3.1-pro-preview", "kw": {"thinking_budget": 8000}},
    "qwen": {"provider": "qwen", "model_id": "qwen-plus-2025-09-11", "kw": {"thinking_budget": 4000}},
}

_ID_RE = re.compile(r"^([A-H]-\d{2})\s*$", re.MULTILINE)
_BLOCK_RE = re.compile(r">>> PASTE BELOW >>>\s*\n(.*?)\n<<< PASTE ABOVE <<<", re.DOTALL)


def extract_questions(path: Path) -> dict[str, str]:
    """{question_id: paste-block text only}, no reference/gold-adjacent lines."""
    text = path.read_text(encoding="utf-8")
    ids = [(m.group(1), m.start()) for m in _ID_RE.finditer(text)]
    out = {}
    for i, (qid, start) in enumerate(ids):
        end = ids[i + 1][1] if i + 1 < len(ids) else len(text)
        chunk = text[start:end]
        bm = _BLOCK_RE.search(chunk)
        if not bm:
            print(f"  WARNING: no paste block found for {qid} in {path.name}, skipping")
            continue
        out[qid] = bm.group(1).strip()
    return out


def load_all(condition_suffix: str) -> dict[str, str]:
    """condition_suffix: 'Unaugmented' or 'Augmented'. Returns {qid: prompt_text} across all 8 categories."""
    out = {}
    for cat, prefix in CATEGORY_FILES.items():
        path = HERE / f"{prefix}_{condition_suffix}.txt"
        out.update(extract_questions(path))
    return out


def rebuild_consolidated(model_key: str, cond_label: str, raw_dir: Path, all_ids: list[str]):
    """Rebuild the .txt/.json consolidated files from whatever raw per-question
    files exist right now, in question-ID order, skipping ones not yet run."""
    transcript_lines = []
    detail = {}
    for qid in all_ids:
        raw_path = raw_dir / f"{qid}.json"
        if not raw_path.exists():
            continue
        out = json.loads(raw_path.read_text(encoding="utf-8"))
        transcript_lines.append(f"{qid}\n{out.get('final_text') or ''}\n")
        detail[qid] = out
    OUT_DIR.mkdir(exist_ok=True)
    txt_path = OUT_DIR / f"{model_key}_{cond_label}.txt"
    txt_path.write_text("\n".join(transcript_lines), encoding="utf-8")
    json_path = OUT_DIR / f"{model_key}_{cond_label}.json"
    json_path.write_text(json.dumps(detail, indent=2, default=str), encoding="utf-8")
    return txt_path, len(detail)


def run_all(model_key: str, condition: str, limit: int | None, categories: list[str] | None, force: bool):
    spec = MODELS[model_key]
    suffix = "Augmented" if condition == "surfmcp" else "Unaugmented"
    all_questions = load_all(suffix)
    all_ids = sorted(all_questions)  # full 42 -- used to keep consolidated output in order
    ids = all_ids
    if categories:
        ids = [qid for qid in ids if qid[0] in categories]
    if limit:
        # take `limit` from each category letter, preserving order
        seen = {}
        kept = []
        for qid in ids:
            cat = qid[0]
            seen.setdefault(cat, 0)
            if seen[cat] < limit:
                kept.append(qid)
                seen[cat] += 1
        ids = kept

    system = prompts.build_system(condition)
    cond_label = "augmented" if condition == "surfmcp" else "unaugmented"
    raw_dir = OUT_DIR / "raw" / f"{model_key}_{cond_label}"
    raw_dir.mkdir(parents=True, exist_ok=True)

    def _already_succeeded(qid: str) -> bool:
        p = raw_dir / f"{qid}.json"
        if not p.exists():
            return False
        try:
            return json.loads(p.read_text(encoding="utf-8")).get("error") is None
        except (json.JSONDecodeError, OSError):
            return False  # a corrupt/partial raw file is not a success -- retry it

    todo = [qid for qid in ids if force or not _already_succeeded(qid)]
    skipped = len(ids) - len(todo)
    print(f"\n=== {model_key} / {cond_label} : {len(ids)} selected"
          f"{f', {skipped} already succeeded (skipping)' if skipped else ''} ===")

    for qid in todo:
        user = all_questions[qid]
        t0 = time.time()
        try:
            out = providers.run(spec["provider"], spec["model_id"], system, user, condition, **spec["kw"])
        except Exception as e:
            out = {"final_text": "", "thinking": "", "tool_calls": [], "raw": {},
                   "error": f"{type(e).__name__}: {e}"}
        dt = time.time() - t0
        n_tools = len(out.get("tool_calls", []))
        err = out.get("error")
        flag = "ERR" if err else "OK "
        print(f"  [{flag}] {qid}  tools={n_tools}  {dt:.1f}s" + (f"  | {err}" if err else ""), flush=True)
        (raw_dir / f"{qid}.json").write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")

    txt_path, n_total = rebuild_consolidated(model_key, cond_label, raw_dir, all_ids)
    print(f"  wrote {txt_path} ({n_total}/{len(all_ids)} questions answered so far)")
    return txt_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="gemini,qwen")
    ap.add_argument("--conditions", default="unaugmented,surfmcp")
    ap.add_argument("--limit", type=int, default=None,
                     help="cap questions PER CATEGORY LETTER (smoke test); omit for the full 42")
    ap.add_argument("--categories", default=None,
                     help="comma list of category letters to restrict to, e.g. A,B,C "
                          "(chunk a full run into pieces that fit a foreground time budget)")
    ap.add_argument("--force", action="store_true",
                     help="re-run questions even if a raw response file already exists for them")
    args = ap.parse_args()

    models = [m.strip() for m in args.models.split(",") if m.strip()]
    conditions = [c.strip() for c in args.conditions.split(",") if c.strip()]
    categories = [c.strip().upper() for c in args.categories.split(",")] if args.categories else None
    for m in models:
        if m not in MODELS:
            raise SystemExit(f"unknown model '{m}'. known: {list(MODELS)}")

    written = []
    for m in models:
        for c in conditions:
            written.append(run_all(m, c, args.limit, categories, args.force))

    print("\nAll runs done. Grade with:")
    grade_args = " ".join(f"--transcript {p.stem.replace('_', '/')} {p}" for p in written)
    print(f"  python grade_all.py {grade_args}")


if __name__ == "__main__":
    main()
