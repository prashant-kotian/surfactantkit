# SurfBench

The benchmark half of the SurfactantKit paper: questions designed to measure how well LLMs (unaugmented vs. tool-augmented via SurfMCP) perform on surfactant and interfacial-science calculations.

See `taxonomy.md` for the full category/tool mapping and design principles.

## Status

**Tier 1 (tool-mapped, 315 questions): done, all 25 SurfactantKit tools covered, every tool at 8+ questions.** Generated programmatically -- every gold answer comes from calling the real library function, not hand-derivation, so gold-answer correctness inherits directly from SurfactantKit's own 106-test, 8-literature-system validation.

**Tier 2 (non-tool control categories, 206 questions): done**, all 6 categories at their exact taxonomy target (J:40, K:40, L:30, M:40, N:30, O:26). Hand-curated (not generator scripts, since there's no formula to compute a gold answer from for classification/recall/conceptual questions), except Category L (unit conversion) which *is* script-generated since conversion factors are exact, cited SI/CGS identities, not literature-dependent values. Category K (property recall) draws every gold value directly from `literature_validation_notes.md` -- i.e. numbers this project already independently verified against real source papers -- rather than fresh recall, since that's the one category where sourcing quality directly determines correctness. See `tier2/build_tier2_bank.py` to regenerate; writes `question_bank_tier2.json`.

Combined Tier 1 + Tier 2: 315 + 206 = **521 total questions** (the original taxonomy's 294+206=500 target was superseded when Tier 1 was rebalanced upward earlier in the project; see the rebalancing log above).

**Real bug caught building Tier 2, worth flagging for anyone extending this bank**: the grading method `category_match` was originally built for exactly one purpose (Tier 1's no-solution traps -- checking whether a model's answer contains no-solution language). Reusing the same `grading_method` string for Tier 2's classification/keyword questions (checking whether the answer contains a specific short phrase like "anionic" or "synergistic") silently broke, because `grade()` used to route on `grading_method == "category_match"` alone and always ran the no-solution check -- meaning every Tier 2 question would have graded as wrong regardless of what the model actually said. Caught by the same offline synthetic-correct/synthetic-wrong test used to validate Tier 1's grading originally; fixed by routing on `trap_type == "no_solution"` specifically, with a separate fuzzy-keyword-match path for genuine classification questions. `benchmark/harness/grading.py`'s `grade()` now takes a `trap_type` argument for this reason -- don't drop it when calling `grade()` from anywhere new.

## Rebalancing pass (first draft -> current)

- `entropy_of_micellization`: 4 -> 12 questions. Root cause was a `zip()` over two case lists that silently capped iteration at the shorter list's length; switched to the full cross-product (`itertools.product`), which is also more question variety, not just more count.
- `gibbs_area_per_molecule`: 6 -> 12 questions.
- Unit-trap questions: 11 -> 16 (added a viscosity Pa.s-vs-mPa.s trap to the capillary number questions).
- `no_solution` traps deliberately left at 7 rather than padded to hit a target ratio -- the existing ones (Rubingh's narrow valid-root range, Davies' unverified HLB groups, MSR below the CMC) are genuinely subtle "looks solvable but isn't" cases; adding more via obviously-degenerate inputs (e.g. two identical temperatures) would lower average trap quality just to move a count.

## Known gaps in Tier 1, still open

- `tanford_chain_geometry`'s reported count (35) looks high but is not a real imbalance -- it's a *prerequisite* tool for the CPP and aggregation-number chain questions too (15 direct + 12 CPP-chain + 8 agg-chain = 35), not 35 standalone questions. Flagging this so it isn't "fixed" by mistake in a future pass.
- One real bug caught and fixed during generation: two of the electrostatics "compute ionic strength from ion concentrations" test cases used exact-stoichiometry salt ratios (CaCl2, Na2SO4) where naive concentration-summing accidentally produces the *same* number as the correct charge-weighted formula -- meaning they would have silently failed to test the thing they were designed to test. Caught by directly computing both ways before finalizing, not assumed correct. Replaced with non-stoichiometric ion mixes verified to actually discriminate.
- Two unit-trap questions originally included in-question hints that gave away the trap ("remember the formula needs Kelvin", "note: convert first") -- removed, since a hinted trap doesn't test anything.
- ~~The newer 17 tools (everything past the original Rubingh/Clint work) have not had the same depth of independent literature cross-validation as the original 8-system sweep~~ -- resolved as of the 2026-09-02 literature validation pass (rounds 2-4, see `literature_validation_notes.md`): CPP/Tanford, adsorption (Gibbs/Szyszkowski), wetting, electrostatics, dynamics, and thermodynamics (ΔG, ΔH) all now have real-paper cross-checks, most within 0.5% agreement. Two genuine open gaps remain, not swept under the rug: `molar_solubilization_ratio` has no confirmed real-paper match yet (candidates found so far use a different experimental regime than the function's convention), and the Szyszkowski `K` constant has no cleanly-sourced literature value (one candidate's extracted digits were internally inconsistent, not used). See the main repo README's Validation section for the current per-tool breakdown.
- Contamination check: `question_text` must never contain author/paper attribution ("X et al."), since that hands the model a search query instead of a computation. Audited all 9 generators -- `wetting_gen.py` (3 cases: Murtaza, Gaynanova, Ahmad Wazir) and `solubilization_gen.py` (2 cases: Patel, Yadav) were interpolating attributed labels straight into the prompt text; stripped to generic descriptions. `mixed_micelle_gen.py`'s literature attributions (Muherei, Azum, Liu, McLachlan, etc.) were already clean -- they only ever went into `source_note`, never `question_text`. This check needs re-running on any future generator before it ships. Separately, for the planned real-published-data validation pilot (feeding actual paper CMC/alpha values through a live AI+MCP session), citation-stripping isn't sufficient on its own -- the exact numeric combination is itself a fingerprint a model could recall or look up. That test needs its own controls: no web-search tool on any condition, jittered (not exact) numeric inputs with freshly-computed gold answers, and cross-checking model answers against the *paper's own* (sometimes formula-inconsistent, per `literature_validation_notes.md`) reported values as a memorization tripwire.

## Model IDs for the scored run (pinned 2026-09-02)

Confirmed live against each provider's own API (models.list where available, or a
direct minimal call), not guessed or copied from docs:

- **Claude**: `claude-opus-4-8` — fixed id, no separate dated alias exists.
- **GPT**: `gpt-5.6-sol` — fixed id, exact match to the model originally chosen; unlike its `gpt-5`/`gpt-5.1`.../`gpt-5.5` siblings it has no date-suffixed variant, which appears to be this tier's own naming convention rather than an oversight.
- **Qwen**: `qwen-plus-2025-09-11` — a real dated snapshot, confirmed responsive. Note: `qwen3-plus`, mentioned earlier in the project as the target model, does **not** exist as an actual DashScope model (`InvalidParameter` on call) -- that was an imprecise recollection, not a real id. The correct target has always been Qwen Plus; this is its current dated snapshot.
- **Gemini**: `gemini-2.5-pro` / `gemini-2.5-flash` — **no dated snapshot alias is exposed on this account** (checked via `models.list`). Gemini is rolling-only here. This is a genuine, disclosed reproducibility limitation for the methods section, not something further digging is expected to resolve -- state plainly that the Gemini results reflect whatever `gemini-2.5-pro` resolved to on the run date.

## Regenerating

```bash
cd benchmark
python build_question_bank.py
```

Rerunning is idempotent and safe -- deterministic inputs, no randomness, so gold answers never drift unless a generator or the underlying library changes.
