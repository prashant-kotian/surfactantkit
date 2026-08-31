# SurfBench

The benchmark half of the SurfactantKit paper: questions designed to measure how well LLMs (unaugmented vs. tool-augmented via SurfMCP) perform on surfactant and interfacial-science calculations.

See `taxonomy.md` for the full category/tool mapping and design principles.

## Status

**Tier 1 (tool-mapped, 301 questions): done, all 25 SurfactantKit tools covered.** Generated programmatically -- every gold answer comes from calling the real library function, not hand-derivation, so gold-answer correctness inherits directly from SurfactantKit's own 106-test, 8-literature-system validation.

**Tier 2 (non-tool control categories, ~206 questions planned): not started.** These need hand-curation from literature/domain knowledge rather than a generator script, since there's no formula to compute a gold answer from for classification/recall/conceptual questions.

## Known gaps in Tier 1, flagged for a refinement pass

- **Per-tool coverage is uneven** (4 to 35 questions per tool, vs. the taxonomy's target of ~12 each). `entropy_of_micellization` (4) and `gibbs_area_per_molecule` (6) are thinnest; `tanford_chain_geometry` (35) and `debye_screening_length` (21) are thickest. Not wrong, but worth rebalancing before this is called final.
- **Unit-trap questions are underrepresented**: 11/301 (~3.7%) vs. the taxonomy's ~15-20% target. `no_solution` traps are similarly thin (7/301, ~2.3%).
- One real bug already caught and fixed during generation: two of the electrostatics "compute ionic strength from ion concentrations" test cases used exact-stoichiometry salt ratios (CaCl2, Na2SO4) where naive concentration-summing accidentally produces the *same* number as the correct charge-weighted formula -- meaning they would have silently failed to test the thing they were designed to test. Caught by directly computing both ways before finalizing, not assumed correct. Replaced with non-stoichiometric ion mixes that were verified to actually discriminate.
- Two unit-trap questions originally included in-question hints that gave away the trap ("remember the formula needs Kelvin", "note: convert first") -- removed, since a hinted trap doesn't test anything.

## Regenerating

```bash
cd benchmark
python build_question_bank.py
```

Rerunning is idempotent and safe -- deterministic inputs, no randomness, so gold answers never drift unless a generator or the underlying library changes.
