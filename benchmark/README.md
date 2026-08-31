# SurfBench

The benchmark half of the SurfactantKit paper: questions designed to measure how well LLMs (unaugmented vs. tool-augmented via SurfMCP) perform on surfactant and interfacial-science calculations.

See `taxonomy.md` for the full category/tool mapping and design principles.

## Status

**Tier 1 (tool-mapped, 315 questions): done, all 25 SurfactantKit tools covered, every tool at 8+ questions.** Generated programmatically -- every gold answer comes from calling the real library function, not hand-derivation, so gold-answer correctness inherits directly from SurfactantKit's own 106-test, 8-literature-system validation.

**Tier 2 (non-tool control categories, ~206 questions planned): not started.** These need hand-curation from literature/domain knowledge rather than a generator script, since there's no formula to compute a gold answer from for classification/recall/conceptual questions.

## Rebalancing pass (first draft -> current)

- `entropy_of_micellization`: 4 -> 12 questions. Root cause was a `zip()` over two case lists that silently capped iteration at the shorter list's length; switched to the full cross-product (`itertools.product`), which is also more question variety, not just more count.
- `gibbs_area_per_molecule`: 6 -> 12 questions.
- Unit-trap questions: 11 -> 16 (added a viscosity Pa.s-vs-mPa.s trap to the capillary number questions).
- `no_solution` traps deliberately left at 7 rather than padded to hit a target ratio -- the existing ones (Rubingh's narrow valid-root range, Davies' unverified HLB groups, MSR below the CMC) are genuinely subtle "looks solvable but isn't" cases; adding more via obviously-degenerate inputs (e.g. two identical temperatures) would lower average trap quality just to move a count.

## Known gaps in Tier 1, still open

- `tanford_chain_geometry`'s reported count (35) looks high but is not a real imbalance -- it's a *prerequisite* tool for the CPP and aggregation-number chain questions too (15 direct + 12 CPP-chain + 8 agg-chain = 35), not 35 standalone questions. Flagging this so it isn't "fixed" by mistake in a future pass.
- One real bug caught and fixed during generation: two of the electrostatics "compute ionic strength from ion concentrations" test cases used exact-stoichiometry salt ratios (CaCl2, Na2SO4) where naive concentration-summing accidentally produces the *same* number as the correct charge-weighted formula -- meaning they would have silently failed to test the thing they were designed to test. Caught by directly computing both ways before finalizing, not assumed correct. Replaced with non-stoichiometric ion mixes verified to actually discriminate.
- Two unit-trap questions originally included in-question hints that gave away the trap ("remember the formula needs Kelvin", "note: convert first") -- removed, since a hinted trap doesn't test anything.
- The newer 17 tools (everything past the original Rubingh/Clint work) have not had the same depth of independent literature cross-validation as the original 8-system sweep -- most rely on textbook reference values, literature ranges, or mathematical round-trip tests rather than an external paper match. See the main repo README's Validation section for the current breakdown by tool.

## Regenerating

```bash
cd benchmark
python build_question_bank.py
```

Rerunning is idempotent and safe -- deterministic inputs, no randomness, so gold answers never drift unless a generator or the underlying library changes.
