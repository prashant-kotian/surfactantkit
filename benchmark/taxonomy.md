# SurfBench question taxonomy

500 questions total, split so every one of SurfactantKit's 25 MCP tools gets systematic, multi-question coverage (Tier 1), plus a control group of non-tool categories (Tier 2) that measures whether tool augmentation is neutral where it should be -- pure recall/classification isn't something a calculation tool can help with, and showing that plainly is part of the paper's honesty.

## Tier 1 -- tool-mapped categories (294 questions)

Each category maps to one SurfactantKit module. Question count per category is roughly proportional to the number of distinct tools in that module (~12 questions per tool), so no tool is represented by only one or two questions (too easy for a fluke pass/fail to dominate that tool's apparent score).

| Category | Tools tested | Module | Count |
|---|---|---|---|
| A. Mixed micelle theory | `clint_ideal_cmc`, `rubingh_solve`, `rubingh_activity_coefficients`, `excess_free_energy`, `rosen_monolayer_solve`, `corrin_harkins_predict` | mixed_micelle | 60 |
| B. Adsorption | `gibbs_surface_excess`, `gibbs_area_per_molecule`, `szyszkowski_predict_surface_tension` | adsorption | 36 |
| C. HLB | `hlb_from_mw`, `hlb_from_groups` | hlb | 24 |
| D. Molecular geometry | `tanford_chain_geometry`, `critical_packing_parameter`, `aggregation_number` | cpp | 36 |
| E. Electrostatics | `debye_screening_length`, `zeta_potential` | electrostatics | 24 |
| F. Dynamics | `hydrodynamic_radius` | dynamics | 15 |
| G. Thermodynamics | `counterion_binding_degree`, `gibbs_free_energy_micellization`, `vant_hoff_enthalpy`, `entropy_of_micellization` | thermodynamics | 48 |
| H. Wetting / EOR | `wetting_work_of_adhesion`, `wetting_spreading_coefficient`, `eor_capillary_number` | wetting | 36 |
| I. Solubilization | `molar_solubilization_ratio` | solubilization | 15 |
| **Subtotal** | | | **294** |

## Tier 2 -- non-tool control categories (206 questions)

No SurfactantKit tool applies to these by design. They exist to show tool augmentation is *neutral* here (doesn't help, shouldn't hurt) -- if augmentation appeared to help uniformly across every category including these, that would itself be suspicious and worth investigating rather than reporting at face value.

| Category | What it tests | Count |
|---|---|---|
| J. Classification & nomenclature | Gemini/bolaform/amidoamine/zwitterionic identification from structure or description | 40 |
| K. Property recall | Known literature CMC/HLB/etc. values -- direct hallucination-rate measurement | 40 |
| L. Unit conversion & dimensional reasoning | mN/m <-> dyn/cm, mM <-> mol/kg, cP <-> Pa.s, etc. -- simple enough that a tool shouldn't be needed | 30 |
| M. Conceptual theory | "Why does X happen" / mechanism explanations -- Clint/Rubingh/Gibbs/CPP interpretation, not calculation | 40 |
| N. Gemini-specific structural reasoning | Spacer length/type effects, architecture comparisons -- ties to the PhD's own domain | 30 |
| O. Real-world application | Formulation/EOR/solubilization reasoning without a clean numeric answer | 26 |
| **Subtotal** | | | **206** |

**Total: 500**

## Question design principles (apply across Tier 1)

1. **Gold answers are computed by SurfactantKit itself, not hand-derived.** Each Tier 1 question is generated from a template + parameter set; a generator script calls the actual library function to produce the gold answer. This guarantees gold-answer correctness is exactly as good as the library's own test suite (106 passing tests, 8 literature-validated systems) and makes the question bank trivial to regenerate if a formula is ever revised.
2. **Realistic parameter ranges.** Where possible, parameters are drawn from or near the 8 literature systems in `literature_validation_notes.md` (real pure CMCs, real chain lengths, real ionic strengths) rather than arbitrary numbers, so questions read like real problems, not textbook plug-ins.
3. **Deliberate unit traps, a meaningful minority (~15-20%) of each category.** E.g. giving a diffusion coefficient in m²/s when the tool expects cm²/s, or a concentration in M when mM is expected. This directly targets the failure mode the whole project is about -- catching whether a model (with and without the tool) notices and converts correctly, rather than only testing "can you plug numbers into a formula."
4. **Multi-tool chains, a meaningful minority (~15-20%) of each category** where it's natural (e.g. Tanford geometry -> CPP -> morphology; Gibbs slope -> Szyszkowski prediction; counterion binding -> Gibbs free energy). Tests genuine tool orchestration, not single-call lookup.
5. **"No valid answer" trap questions, a small minority (~5-10%).** E.g. Rubingh inputs with no real root, a Davies group not in the verified table, an out-of-range packing parameter. Gold answer is "correctly refuses / reports no solution," not a number -- directly tests whether a model (and the tool) resist producing a plausible-looking wrong answer, which is the whole thesis of this project.
6. **Grading**: numeric-tolerance match (typically relative 1-2%, wider for quantities with inherent approximation like aggregation number) for standard questions; exact string/category match for trap and classification questions.

## Directory layout

```
benchmark/
  taxonomy.md              (this file)
  schema.py                (Question dataclass, shared by all generators)
  generators/
    mixed_micelle_gen.py    (Category A)
    adsorption_gen.py        (Category B)
    hlb_gen.py                (Category C)
    geometry_gen.py            (Category D)
    electrostatics_gen.py       (Category E)
    dynamics_gen.py               (Category F)
    thermodynamics_gen.py          (Category G)
    wetting_gen.py                   (Category H)
    solubilization_gen.py             (Category I)
  build_question_bank.py    (runs all generators, writes the full bank to JSON)
  question_bank_pilot.json  (first pilot batch, for review before scaling)
```

Tier 2 (non-tool) categories are written by hand/curated from literature (not generator scripts, since there's no formula to compute a gold answer from), except Category L (unit conversion), which is script-generated since conversion factors are exact SI/CGS identities. Done -- see `tier2/` for the 6 category modules and `tier2/build_tier2_bank.py` for the assembly script, writing `question_bank_tier2.json`.
