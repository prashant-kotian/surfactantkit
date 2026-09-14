"""Category A -- chain-length series entanglement, real held-out compound.

Different entangled skill chain than Round 3's curve-orchestrator: no
tensiometry curve at all here. The model is given N-1 REAL literature CMC
values across a real homologous series (each with its own real citation,
sourced from this PhD project's own already-validated SurfQSPR
SEED_DATASET -- see H:\\CodeProjects\\SurfQSPR\\src\\surfqspr\\dataset.py)
and must fit the real Klevens/Stauff-Klevens log-linear chain-length
relationship, then predict the CMC of ONE held-out real compound (given
only by name/SMILES/chain length, not its real CMC).

Gold is DOUBLE-SOURCED, per the explicit 2026-09-14 requirement: (1) the
real tool's own fit-and-predict output (SurfQSPR-MCP's
fit_and_predict_cmc_from_points, which fits ONLY on the given points, no
internal dataset access -- see mcp_server.py's own docstring for why
predict_cmc() would leak the held-out answer), AND (2) the real held-out
compound's own independently-published literature CMC (already sitting in
SEED_DATASET with its own citation) -- both included in gold so grading
can be checked against either the tool's prediction or the literature's
own real value, and the real percent difference between them is reported
honestly rather than hidden.

Design choices made specifically to avoid an easy, ungrounded pass (per
explicit 2026-09-14 user instruction -- "don't keep loose ends... turn our
expected unaugmented output towards more than 60%"):
  - loo_cv_mean_pct_error is REQUIRED in the schema. This needs N real
    leave-one-out sub-fits (refit N times, each holding out one more
    point) -- not something a model can do accurately by hand across 7-8
    points without either a tool or a lot of careful arithmetic; Round 3
    already found models happily approximate fit-quality numbers rather
    than compute them for real.
  - slope/intercept/r_squared are graded to real fitted precision (20%
    relative tolerance, same convention as Round 3), not just "is the
    predicted CMC roughly right" -- closes off the shortcut of eyeballing
    a rule-of-thumb ("each 2 carbons roughly halves CMC") and getting the
    point prediction close by luck while the fit parameters themselves are
    wrong.
  - 3 of 6 questions deliberately remove a SERIES BOUNDARY point (not an
    interior one), turning the query into a genuine extrapolation --
    testing whether the model recognizes and flags degraded confidence,
    not just whether it can also compute a number there.
  - Each given point carries its real confidence tier ("validated" vs
    "textbook", copied verbatim from SEED_DATASET) -- present as
    real research input, not filtered/cleaned first.

Usage:
  python build_category_A.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(r"H:\CodeProjects\SurfQSPR\src")))
from surfqspr.predict import fit_and_predict_cmc_from_points
from surfqspr.dataset import SEED_DATASET

_SULFATE = {d.n_carbons: d for d in SEED_DATASET if "sulfate" in d.name and "cholate" not in d.name}
_TAB = {d.n_carbons: d for d in SEED_DATASET if "trimethylammonium bromide" in d.name}

QUESTIONS_SPEC = [
    # (id, series_dict, series_label, holdout_n, held_out_smiles_name)
    ("A-01", _SULFATE, "sodium alkyl sulfate (anionic)", 11),
    ("A-02", _SULFATE, "sodium alkyl sulfate (anionic)", 16),
    ("A-03", _SULFATE, "sodium alkyl sulfate (anionic)", 8),
    ("A-04", _TAB, "alkyltrimethylammonium bromide (cationic)", 12),
    ("A-05", _TAB, "alkyltrimethylammonium bromide (cationic)", 6),
    ("A-06", _TAB, "alkyltrimethylammonium bromide (cationic)", 14),
]


def build_question(qid, series, label, holdout_n):
    given_points = [(n, d.cmc_mM) for n, d in series.items() if n != holdout_n]
    held_out_d = series[holdout_n]
    r = fit_and_predict_cmc_from_points(given_points, holdout_n)
    real_cmc = held_out_d.cmc_mM
    pct_diff_tool_vs_lit = 100 * abs(r.predicted_cmc_mM - real_cmc) / real_cmc

    given_rows = [
        {"n_carbons": n, "name": series[n].name, "cmc_mM": series[n].cmc_mM,
         "confidence": series[n].confidence, "source": series[n].source}
        for n in sorted(series) if n != holdout_n
    ]

    return {
        "id": qid,
        "series_label": label,
        "holdout_n_carbons": holdout_n,
        "held_out_compound": {"name": held_out_d.name, "smiles": held_out_d.smiles,
                               "real_literature_cmc_mM": real_cmc, "confidence": held_out_d.confidence,
                               "source": held_out_d.source},
        "given_points": given_rows,
        "gold": {
            "predicted_cmc_mM": r.predicted_cmc_mM,
            "slope": r.slope,
            "intercept": r.intercept,
            "r_squared": r.r_squared,
            "in_applicability_domain": r.in_applicability_domain,
            "training_range": list(r.training_range),
            "loo_cv_mean_pct_error": r.loo_cv_mean_pct_error,
        },
        "gold_cross_check": {
            "real_literature_cmc_mM": real_cmc,
            "tool_predicted_cmc_mM": r.predicted_cmc_mM,
            "pct_difference_tool_vs_literature": pct_diff_tool_vs_lit,
            "note": ("Extrapolation (query outside the given points' range)" if not r.in_applicability_domain
                     else "Interpolation (query within the given points' range)")
                    + f" -- tool prediction and the real independently-published literature CMC "
                      f"agree to within {pct_diff_tool_vs_lit:.1f}%.",
        },
    }


UNAUG_TEMPLATE = """A-{n:02d}
[reference only, do not paste -- real held-out compound: {held_out_name}, real literature CMC {held_out_cmc} mM ({held_out_conf}); withheld from the prompt below, used only for grading]

>>> PASTE BELOW >>>
Below is a real set of published critical micelle concentration (CMC) values for a homologous {series_label} surfactant series (same headgroup class, varying alkyl chain length, all at 25C), each with its own real literature source. Using ONLY these points -- do not recall or assume any CMC value for any compound in this series from your own training, including the one asked about below -- fit the real relationship between chain length and CMC for this series, and use that fit to predict the CMC of the following compound:

Query compound: {query_name}
SMILES: {query_smiles}
Chain length: {query_n} carbons

Given real literature data points (n_carbons, CMC in mM, confidence tier, source):
{points_block}

State explicitly what fitting method you actually used (e.g. real least-squares linear regression on log10(CMC) vs. chain length, vs. an approximate/eyeballed estimate) -- this is part of what's being assessed, not just the final number. Report a leave-one-out cross-validation (LOO-CV) mean percent error for this fit if you can genuinely compute one (refit the line N times, each time holding out one more of the given points and testing on it) -- do not report a plausible-looking LOO-CV number you have not actually computed; report it as not computed if you didn't. State explicitly whether the query chain length falls within or outside the given points' own chain-length range, and if outside, treat your prediction as a genuine extrapolation, not equally confident as an interpolation would be.

End your response with a block in exactly this form (use null for anything not genuinely determined):

FINAL_JSON:
{{
  "predicted_cmc_mM": <number or null>,
  "slope": <number or null, the real fitted log10(CMC) vs n_carbons slope>,
  "intercept": <number or null>,
  "r_squared": <number or null>,
  "in_applicability_domain": true or false,
  "training_range": [<min_n_carbons>, <max_n_carbons>],
  "loo_cv_mean_pct_error": <number or null, ONLY if you genuinely computed real leave-one-out sub-fits>,
  "fitting_method": "real_least_squares_regression" or "approximate_estimate" or other short description,
  "undeterminable": ["short phrase + reason", ...]
}}
<<< PASTE ABOVE <<<
"""

AUG_TEMPLATE = """A-{n:02d}
[reference only, do not paste -- real held-out compound: {held_out_name}, real literature CMC {held_out_cmc} mM ({held_out_conf}); withheld from the prompt below, used only for grading]

>>> PASTE BELOW >>>
You have access to the surfqspr-mcp server's fit_and_predict_cmc_from_points tool -- call it to analyze the real data below. Do not compute anything by hand; use the tool, and report exactly what it returns. Do not recall or assume a CMC value for the query compound from your own training -- pass it only as query_n_carbons to the tool.

Below is a real set of published critical micelle concentration (CMC) values for a homologous {series_label} surfactant series (same headgroup class, varying alkyl chain length, all at 25C), each with its own real literature source.

Query compound: {query_name}
SMILES: {query_smiles}
Chain length: {query_n} carbons

Given real literature data points (n_carbons, CMC in mM, confidence tier, source):
{points_block}

Call the tool with points=[[n_carbons, cmc_mM], ...] for exactly the {n_points} given points above (not the query compound), and query_n_carbons={query_n}. Report the tool's own returned values faithfully -- do not override, "correct", or second-guess a value the tool returns, and do not fill in a value the tool itself left null.

End your response with a block in exactly this form:

FINAL_JSON:
{{
  "predicted_cmc_mM": <number or null>,
  "slope": <number or null>,
  "intercept": <number or null>,
  "r_squared": <number or null>,
  "in_applicability_domain": true or false,
  "training_range": [<min_n_carbons>, <max_n_carbons>],
  "loo_cv_mean_pct_error": <number or null>,
  "fitting_method": "tool_call",
  "undeterminable": ["short phrase + reason", ...]
}}
<<< PASTE ABOVE <<<
"""


def render(template, i, q):
    points_block = "\n".join(
        f"  n_carbons={p['n_carbons']}, CMC={p['cmc_mM']} mM, [{p['confidence']}] {p['source']}"
        for p in q["given_points"]
    )
    return template.format(
        n=i, held_out_name=q["held_out_compound"]["name"],
        held_out_cmc=q["held_out_compound"]["real_literature_cmc_mM"],
        held_out_conf=q["held_out_compound"]["confidence"],
        series_label=q["series_label"],
        query_name=q["held_out_compound"]["name"], query_smiles=q["held_out_compound"]["smiles"],
        query_n=q["holdout_n_carbons"], points_block=points_block, n_points=len(q["given_points"]),
    )


def main():
    questions = [build_question(*spec) for spec in QUESTIONS_SPEC]
    (HERE / "category_A_questions.json").write_text(json.dumps(questions, indent=2), encoding="utf-8")

    unaug = "\n\n".join(render(UNAUG_TEMPLATE, i, q) for i, q in enumerate(questions, 1))
    (HERE / "CategoryA_Unaugmented.txt").write_text(unaug, encoding="utf-8")

    aug = "\n\n".join(render(AUG_TEMPLATE, i, q) for i, q in enumerate(questions, 1))
    (HERE / "CategoryA_Augmented.txt").write_text(aug, encoding="utf-8")

    print(f"Built {len(questions)} Category A questions.")
    for q in questions:
        print(f"  {q['id']}: holdout n={q['holdout_n_carbons']} ({q['held_out_compound']['name']}), "
              f"tool_pred={q['gold']['predicted_cmc_mM']:.3f} mM, real_lit={q['held_out_compound']['real_literature_cmc_mM']} mM, "
              f"agree={q['gold_cross_check']['pct_difference_tool_vs_literature']:.1f}%, "
              f"in_domain={q['gold']['in_applicability_domain']}")


if __name__ == "__main__":
    main()
