# SurfactantKit Capability Audit -- for Paper 3 ground-zero redesign

Purpose: before redesigning the LLM benchmark, establish exactly what
SurfactantKit can and cannot do today, and which of its real capabilities are
genuinely tool-dependent (an unaugmented LLM should plausibly fail or
hallucinate) versus arithmetic an LLM can already fake convincingly by hand
(already empirically confirmed in the old 42-question benchmark: chain-length
regression, Debye length, HLB/Griffin, plain wetting work-of-adhesion, Tanford
CPP all scored 83-100% unaugmented across all 4 models tested).

All function names, formulas, and docstring claims below are taken directly
from reading the real source in `src/surfactantkit/` on 2026-09-14, not
recalled from memory. 61 functions are exposed via `mcp_server.py`
(essentially the full library surface -- not a curated subset).

## Trap taxonomy used below

- **GOOD TRAP** -- a real, tool-dependent judgment call or computation an
  unaugmented LLM is likely to get wrong: a genuine model-selection decision
  among competing theories, a real per-system empirical parameter not in
  general training data, real iterative numerical solving, real
  cheminformatics parsing, or a formula obscure/complex enough that
  chain-of-thought reproduction is unreliable.
- **BAD TRAP** -- straightforward formula/arithmetic an LLM with decent
  chain-of-thought already does reliably (confirmed empirically this
  session).
- **NOT YET BUILDABLE** -- capability the toolkit does not have at all.

## Capability catalog by module

### classify.py / orchestrate.py -- structure-driven autonomous pipeline

| Capability | Needs | Trap |
|---|---|---|
| `classify_surfactant_charge_type` | real SMILES | **GOOD** -- real RDKit SMARTS substructure matching (sulfate/sulfonate/phosphate ester, carboxylate, quaternary N, free amines); correctly returns `zwitterionic`/`ambiguous_pH_dependent`/`unparseable` rather than guessing. An LLM without a real parser is pattern-matching on memorized SMILES fragments, not parsing -- error-prone on unfamiliar/branched structures, and it cannot correctly refuse the pH-ambiguous case the way this function does. |
| `derive_all_properties_from_smiles_and_curve` | SMILES + raw pre/post-CMC tensiometry curve | **GOOD** -- chains classification into real CMC break-point detection, Langmuir-vs-Frumkin BIC model selection, Gamma_max/A_min, and deltaG_mic, propagating every genuinely-missing input (electrolyte condition, counterion alpha) as a named gap instead of guessing. This is the single highest-value "give raw SMILES + raw data, no formula told" question shape already built, but the old benchmark never actually asked it this way -- it always pre-fit the curve and handed the slope/CMC to the model. |
| Structural-family detection beyond charge type (gemini/dimeric, biosurfactant/glycolipid, etc.) | -- | **NOT YET BUILDABLE** -- confirmed absent. `classify.py` only returns charge type. |

### curve_analysis.py -- raw-data extraction methods

| Capability | Needs | Trap |
|---|---|---|
| `cmc_from_surface_tension_curve` | raw (C, gamma) points | **GOOD** -- real 2-segment vs. 3-segment (baseline+decline+plateau) BIC model selection, with a documented real bug history (naive `min_points_per_segment` made 3-segment unreachable for realistic 6-10-point curves, silently corrupting CMC/A_min). Also computes C20/pC20 surfactant efficiency, only when the decline genuinely reaches a 20 mN/m drop. An LLM eyeballing a plotted/tabulated curve is very unlikely to reproduce a real BIC-driven segment choice or correctly withhold C20 when undefined. |
| `cmc_from_conductivity_curve` | raw conductivity-vs-concentration | **GOOD** -- same break-point-finding problem, different real experimental technique, different noise/curvature signature. |
| `aggregation_number_from_quenching_curve`, `_from_sls_debye_plot`, `_from_dls`, `_from_svedberg_equation` | 4 distinct real experimental techniques (fluorescence quenching, static light scattering/Debye plot, DLS, sedimentation/Svedberg) | **GOOD** -- an LLM must know WHICH formula applies to which raw data shape; these are not interchangeable, and confusing them (e.g. applying Svedberg's sedimentation-diffusion formula to SLS Debye-plot data) produces a plausible-looking but wrong number. |

### mixed_micelle.py -- 9 real competing mixture theories, no auto-selector

| Model | Trap |
|---|---|
| Clint ideal | BAD alone (simple harmonic mean) but the REAL trap is knowing it's only valid as an ideal-mixing baseline, not a fit |
| Rubingh regular-solution (point + regression) | **GOOD** -- real grid-scan+bisection root solve for x1, non-trivial residual equation |
| Asymmetric Margules (w12/w21) | **GOOD** -- 2-parameter activity-coefficient model, genuinely different structural assumption from Rubingh's single beta |
| EOMMM global fit | **GOOD** -- 2D grid search + golden-section refinement over (w12, w21), real infeasibility-margin handling |
| Rodenas (point + series, quadratic Lagrange derivative) | **GOOD** -- requires a numerical derivative of a fitted curve, a genuinely different computational path from Rubingh |
| Maeda free energy of mixing | **GOOD** -- distinct thermodynamic convention |
| Rosen monolayer | **GOOD** -- different physical picture (surface excess/monolayer, not micellar mole fraction) |
| Motomura ideal composition | BAD alone (closed form) but real as a baseline-comparison input |
| Corrin-Harkins (counterion-effect CMC prediction) | **GOOD** -- real counterion-concentration-dependent extrapolation |

**Confirmed gap:** no orchestrator runs these against one raw mixture dataset
and concludes which model actually fits (no mixture analogue of
`select_isotherm_model`'s BIC comparison). This is the single highest-payoff
build: the raw models already exist and are tested; only the
comparison/selection layer is missing.

### thermodynamics.py -- real convention/model choices, not just plug-in

| Capability | Trap |
|---|---|
| `gibbs_free_energy_micellization` counterion_factor convention | **GOOD** -- explicit `(2-beta)`/`(1+alpha)` convention that must come from independent measurement (e.g. a real EPR/conductivity value like Bales et al.'s alpha=0.272 for SDS), never guessed; the function raises rather than defaulting. |
| `mass_action_free_energy_of_micellization` vs. the pseudo-phase-separation formula above | **GOOD** -- a genuine model CHOICE: mass-action's correction terms are non-negligible below n~50-100; an LLM must recognize when the simpler pseudo-phase formula is NOT a safe approximation. |
| `critical_micellization_degree`, two formulas (`d2alpha_dc2` vs `d2alpha_dlnc2`) | **GOOD** -- both legitimate, answer a differently-posed question; picking wrong is a real, easy-to-miss error, not obviously wrong-looking. |
| `vant_hoff_enthalpy` (2-point) | **GOOD, with history** -- this exact function had a real sign bug (fixed 2026-09-04); Claude's real transcript this session reproduced the identical sign-flip failure mode unaugmented, on the CORRECT convention (prefactor=1), proving this is a genuine, reproducible trap, not a grading artifact. |
| `vant_hoff_multi_point_fit` (deltaH, deltaS, deltaCp simultaneously) | **GOOD** -- real 3x3 linear system solve (Gaussian elimination) on a reparameterized Gibbs-Helmholtz+Cp form; a 2-point estimate cannot substitute for this when Cp is non-negligible, and knowing WHEN 5+ points are needed instead of 2 is itself a judgment call. |
| `counterion_binding_degree` (slope-ratio) | BAD arithmetic once slopes are given, but real per-system data (conductivity slopes) an LLM cannot know without the raw curve. |

### adsorption.py -- isotherm model selection is the whole point

| Capability | Trap |
|---|---|
| `select_isotherm_model` (Langmuir vs Frumkin via BIC) | **GOOD, already proven** -- real competing-model selection; Frumkin always fits at least as well by construction (more parameters), so BIC penalization is the entire correctness question. Already used by `orchestrate.py`. |
| `frumkin_fit_K_and_a` | **GOOD** -- 2D grid+golden-section joint fit, real 2D-condensation regime (`|a|>=2`) disclosed as a genuine physical multi-valued-root case, not a solver bug. |
| `N_FACTOR_BY_SYSTEM_TYPE` (n=1 vs n=2 Gibbs prefactor) | **GOOD, already proven** -- this exact parameter is the one documented case where unaugmented Claude Opus 4.8 got it wrong 5/6 times AND tool-augmented Claude got it wrong 0/6 (a caller silently defaulting a plausible-looking number through). |

### hlb.py -- 3 competing methods, explicit refusal on missing groups

| Capability | Trap |
|---|---|
| Griffin (`hlb_griffin`, mass-ratio) vs. Davies (`hlb_davies`, group-additive) vs. Guo/Rong/Ying ECL refinement (`hlb_davies_guo_ecl`, nonlinear effective chain length) | **GOOD** -- three genuinely different methods with different data requirements and different accuracy regimes (Guo's ECL corrects Davies' known nonlinearity failure for polyethoxylated nonionics specifically); knowing which applies to a given structure is a real judgment call. |
| `hlb_davies` raising `KeyError` on an unlisted group (imidazoline deliberately unresolved; amide/sulfonate/betaine/phosphate group numbers are this project's own DERIVED, not-literature values, gated behind `allow_derived_groups=True`) | **GOOD** -- tests whether an LLM invents a plausible-looking group number instead of correctly reporting "not determinable," exactly the discipline this project enforces on itself. |

### hld.py -- salinity-scan fitting, per-headgroup-class empirical constants

| Capability | Trap |
|---|---|
| `hld_ionic`/`hld_nonionic` (log-linear vs linear-in-salinity, different physical mechanism) | **GOOD** -- picking the wrong salinity dependence (ionic vs nonionic mechanism) is a real, easy, physically-motivated error. |
| `fit_k_and_cc_from_salinity_scan` | **GOOD** -- real linear regression on `ln(S*)` vs EACN from a raw multi-oil scan; k and Cc are genuinely surfactant-class-specific (`CATIONIC_QUAT_K=0.7`, per-compound Cc table) and cannot be recalled from general knowledge. |
| `hlb_from_cc_cationic_quat` | **GOOD** -- an empirical bridge valid ONLY for one headgroup class; applying it elsewhere is a real, plausible-looking error. |

### cpp.py -- Tanford vs Nagarajan (tail length matters or doesn't)

| Capability | Trap |
|---|---|
| `tanford_tail_volume`/`tanford_critical_length` -> `critical_packing_parameter` -> `classify_aggregate_morphology` with a FIXED head area | **BAD TRAP, already confirmed** -- this is exactly what the old benchmark tested and all 4 models did near-perfectly by hand. |
| `nagarajan_equilibrium_area_ionic` (chain-length-DEPENDENT headgroup area via Debye screening) | **GOOD** -- a real, sourced alternative showing the fixed-a0 assumption above is wrong for ionic surfactants; requires chaining a Debye-Huckel kappa calc into an equilibrium-area formula, then INTO the packing parameter -- multi-step, and the "which model applies" judgment (ionic vs nonionic, tail-length-sensitive or not) is the actual test. |

### electrostatics.py / dynamics.py -- regime selection + a genuinely hard formula

| Capability | Trap |
|---|---|
| `henry_function` regime selection (huckel/smoluchowski/ohshima/swan_furst/qin) | **GOOD** -- 5 real regimes, each with a stated validity domain by kappa*a; real micelles (kappa*a ~0.2-4.2) land in the poorly-approximated intermediate zone where the two textbook limits (huckel/smoluchowski) are both wrong -- picking one of those two anyway is the exact failure mode this function's docstring warns against. |
| `electrophoretic_mobility_relaxation_corrected` (Ohshima-Healy-White 1983) | **GOOD, strong candidate** -- infinite-series terms, empirically-fit correction factor, non-monotonic mobility-vs-zeta behavior above ~100mV; essentially impossible for an LLM to reproduce from memory rather than real computation. |
| `perrin_friction_factor` + `hydrodynamic_radius_perrin_corrected` | **GOOD** -- directly explains this session's own F-04 "DOSY vs DLS radius mismatch" puzzle; plain Stokes-Einstein on a rodlike micelle gives an inflated apparent radius unless shape-corrected, and axial_ratio must come from an independent CPP-morphology judgment, not be guessed. |

### wetting.py -- Owens-Wendt vs van Oss-Chaudhury-Good (a great trap)

| Capability | Trap |
|---|---|
| `work_of_adhesion`/`spreading_coefficient` (single-liquid Young-Dupre) | **BAD TRAP, already confirmed** -- Category G/H of the old benchmark, 100% unaugmented across every model. |
| `owens_wendt_solid_surface_energy` (2-component dispersive/polar) | **GOOD** -- needs 2+ liquids with known dispersive/polar surface-tension components, a real least-squares/linear solve. |
| `van_oss_chaudhury_good_solid_surface_energy` (3-component LW/acid/base) | **GOOD, excellent candidate** -- requires EXACTLY 3 independent test liquids and a real 3x3 linear solve whose cross-terms are genuinely counter-intuitive (solid's acid component pairs with each liquid's BASE component, not acid-with-acid) -- a very easy, very plausible-looking mistake for an LLM to make by naive pattern-matching, and a real, checkable model requiring literature-tabulated per-liquid LW/acid/base values that don't exist in a bare contact-angle number. |

### solubilization.py -- straightforward once inputs are real

`molar_solubilization_ratio`, `micelle_water_partition_coefficient`: simple
formulas, but both require an intrinsic water solubility that is a real,
compound-specific measured quantity, not derivable from SMILES/structure
alone -- **GOOD as a "what's missing" refusal test**, BAD as arithmetic once
supplied.

### Confirmed total gaps (checked via full-tree grep, not assumed)

- **No mixture model-selection orchestrator.** 9 real theories exist; nothing
  runs them together and concludes which fits.
- **No structural-family classifier beyond charge type.** No gemini/dimeric,
  biosurfactant/glycolipid, or other structural-family detection from SMILES.
- **Zero graph-digitization capability.** Nothing extracts (x,y) values from
  a plot image anywhere in the codebase.
- **Zero output packaging.** No plotting, no Excel/CSV export, no
  auto-generated narrative conclusions anywhere in `src/`.

## Build-priority ranking for Paper 3's redesigned benchmark

1. **Mixture model-selection orchestrator (highest payoff, lowest cost).**
   The 9 models and their math already exist and are tested; this is
   orchestration, not new science -- directly analogous to
   `select_isotherm_model`'s already-proven BIC pattern. Unlocks an entire
   new category of "give raw mixture SMILES+data, no model named" questions
   that cannot be faked by arithmetic, because the correct MODEL CHOICE is
   the actual test, not the arithmetic once a model is chosen.
2. **Structural-family classifier.** Meaningful but narrower payoff than #1
   for LLM-differentiation purposes -- extends what orchestrate.py can
   automatically route to, but charge-type classification alone already
   supports strong ground-zero questions today.
3. **Output packaging (graphs/Excel/conclusions).** High value for making
   SurfactantKit useful to real researchers (the tool's actual mission), but
   secondary for Paper 3's specific LLM-differentiation goal -- an LLM's
   ability to describe/interpret a result is a different capability than
   whether it can derive the result correctly, and conflating the two would
   muddy the benchmark's construct validity.
4. **Graph digitization.** Real, substantial, standalone computer-vision
   build (closer to a WebPlotDigitizer clone than an extension of existing
   code) -- highest implementation cost of the four, and orthogonal to
   chemistry judgment specifically (it's an image-processing capability, not
   a chemistry-reasoning one). Lowest priority for Paper 3 specifically,
   though real value for the tool's broader research-user mission.

## What a ground-zero-redesigned question looks like

**Old style (already built, confirmed too easy -- 83-100% unaugmented):**
> "Here is CMC data at 3 chain lengths [pre-selected points], here is the
> regression slope/intercept [already implied], predict CMC at n=11."

The model only has to do algebra on numbers the question has effectively
already organized for it.

**New style, pure surfactant (uses `derive_all_properties_from_smiles_and_curve`):**
> "Here is [SMILES], and here is a raw surface-tension-vs-concentration
> table [8-10 raw (C, gamma) pairs, no CMC pre-marked, no ionic/nonionic
> label]. Determine: (a) the surfactant's charge type from structure alone,
> (b) the CMC, (c) whether the pre-CMC adsorption is better described by
> Langmuir or Frumkin, and (d) Gamma_max, A_min. State explicitly anything
> that cannot be determined from what's given (e.g. deltaG_mic requires a
> counterion-binding measurement this data doesn't include) rather than
> guessing it."
No formula is named. No model is pre-chosen. No "ionic" label is given --
the model must derive it from the SMILES the way `classify.py` does. This
directly tests the exact judgment chain `orchestrate.py` already performs and
that no LLM has been asked to replicate unaugmented yet.

**New style, mixture (once the model-selection orchestrator exists):**
> "Here are 2 real surfactant SMILES, their pure CMCs, and a raw series of
> (bulk mole fraction, mixed CMC) points [no model named]. Determine the
> most appropriate mixture theory for this system and report its
> interaction parameter, with your reasoning for why that model (not one of
> the other 8) fits best."
This cannot be faked by arithmetic -- Rubingh, Margules, EOMMM, and Rodenas
all consume the same raw inputs and produce DIFFERENT, individually
plausible-looking numbers; only real model-comparison (the orchestrator's
BIC-style logic) distinguishes a correct answer from a plausible wrong one.

**New style, wetting (uses van Oss-Chaudhury-Good):**
> "Here are contact angles for 3 test liquids [with real literature LW/acid/
> base surface-tension components given] on the same solid. Decompose the
> solid's surface energy into its Lifshitz-van der Waals and Lewis acid/base
> components."
The acid-pairs-with-base cross-term is a natural, plausible-looking mistake
for an LLM to make; getting the RIGHT number requires actually solving the
3x3 system with the correct pairing, not narrating a formula.
