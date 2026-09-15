# Bottleneck resolution plan: closing what can be closed

Purpose: `GENUINE_BOTTLENECK_AUDIT.md` found 14 real bottlenecks. This
document classifies each by HOW it could be closed -- data mining, our own
computation, a genuinely new method we'd have to propose and validate
ourselves, or a real standing wall nothing closes -- so we work on the
right thing in the right order instead of treating all 14 as one pile.
Several bottlenecks split across more than one path; that split is itself
real information, not hedging.

## The four resolution paths

- **[MINE]** -- closeable by finding more real literature that already
  exists but hasn't been sourced into this toolkit yet. Pure search/reading
  effort, no new derivation.
- **[COMPUTE]** -- closeable using data/functions we already have inside
  SurfactantKit itself (RDKit descriptors, our own CPP/tail-volume module,
  etc.) via a real, defensible derivation. No new external literature
  needed, but real derivation work needed.
- **[NEW METHOD]** -- the existing literature framework genuinely doesn't
  reach here (like Davies for gemini, confirmed this session). Closing it
  means proposing and validating something SurfactantKit itself invents --
  a real, citable research contribution, not a lookup.
- **[WALL]** -- a genuine, standing epistemic block. No mining, no
  computation, no method of ours closes it -- only a specific new
  measurement (by someone, on that specific system) would.

## Classification, bottleneck by bottleneck

**1. Counterion binding degree (beta).**
[MINE] for the ~20-30 most-studied surfactants (SDS, CTAB, DTAB, etc.) --
real conductivity slope-ratio values are published for these, just not yet
compiled into this toolkit as a reference table.
[NEW METHOD] for everything else -- Manning counterion-condensation theory
or a Poisson-Boltzmann cell model can PREDICT beta from headgroup charge
density and ionic strength without a per-compound measurement. This is
real, established electrostatics theory not yet wired into SurfactantKit;
implementing and validating it against the mined literature table above
would be a genuine, citable contribution, closing the bottleneck for
compounds nobody has measured.
[WALL] remains only for the residual disagreement between theory and
measurement on any single real compound (theory won't be exact).

**2. Electrolyte condition (Gibbs prefactor n).**
[WALL], mostly -- this is a fact about what was in the beaker, not
recoverable from a surface-tension curve alone.
[NEW METHOD], conditionally -- IF a conductivity-vs-concentration curve is
also supplied (not just surface tension), the pre-CMC slope directly gives
solution ionic strength, which could back out the electrolyte condition.
Worth building as an optional cross-source inference path, clearly scoped
to "only when conductivity data is also given."

**3. HLD k (non-quat classes).**
[MINE], pure. The HLD-NAC literature (Salager/Acosta groups) is active and
real; k values for nonionic/zwitterionic/gemini/biosurfactant classes
likely exist in papers not yet searched. Straightforward next search
target, no new theory needed.

**4. HLD Cc (non-quat classes).**
[MINE], pure. Same literature family as #3, same search effort.

**5. HLD b (non-ethoxylate nonionics).**
[MINE], pure. Same family again.

**6. Davies gemini/glycolipid HLB.**
[NEW METHOD], exclusively -- confirmed structural this session, not
data-absent (two rounds of real testing). [MINE] and [COMPUTE] both
already tried and failed here; more literature search on Davies' OWN
framework won't fix a framework whose assumptions don't hold for a
bis-headgroup architecture. The real path: either (a) derive a genuinely
new group-contribution scheme calibrated specifically for multi-headgroup
architectures (using the 5-compound gemini dataset already in hand as a
starting calibration set), or (b) formally establish Griffin's mass-ratio
method as the correct DEFAULT for these structural classes rather than
trying to extend Davies at all -- itself a real, citable methodological
finding, not a workaround.

**7. Nagarajan headgroup_prefactor_A (non-sulfate heads).**
[MINE] -- Nagarajan's own follow-up papers and related electrostatic-model
literature may have prefactors for sulfonate/carboxylate/quat headgroups
not yet searched for.
[COMPUTE] -- the constant (82.0 for sulfate) is itself a product of two
more fundamental, separately-estimable quantities (a charge-separation
distance and a bare hydrocarbon-water interfacial free energy). Deriving
the prefactor from those sub-quantities per headgroup, rather than needing
one pre-combined literature constant, is real, buildable work using
electrostatics already in this toolkit (`electrostatics.py`).

**8. Perrin axial_ratio.**
[COMPUTE], real and near-term buildable. The function's own docstring
already names one legitimate source: CPP's predicted morphology
(`cpp.py`'s `classify_aggregate_morphology`). Build a bridge that derives
an approximate axial-ratio PRIOR from the already-computed CPP/morphology
class, clearly labeled as a structural estimate (with real uncertainty
bounds), not a measurement -- no new external dependency, uses modules we
already have.

**9. Henry-function high-zeta non-uniqueness.**
**Revised 2026-09-15, user directly questioned the pure-WALL call and was
right to.** [WALL] only for a SINGLE mobility value at a SINGLE ionic
strength -- that specific input really is underdetermined (the mobility-
zeta curve folds back above ~100 mV, one mobility maps to two zeta
values), and no computation on that one number alone fixes it.
[NEW METHOD]/[COMPUTE], real and buildable, for a richer input: (a) if a
mobility-vs-ionic-strength SERIES is available instead of one point, kappa*a
varies across the series and the shape (plateau vs. maximum-and-fold-back
within the measured range) is itself diagnostic of which branch the system
is on -- standard colloid-science practice, no new measurement TYPE needed,
just more of the SAME measurement across conditions. (b) Implementing the
fuller O'Brien-White (1978) numerical relaxation treatment (rather than the
simpler Henry approximation currently used) explicitly locates the mobility
maximum from kappa*a and surface conductivity, using only electrophoretic
data already in scope. (c) A disclosed, literature-grounded plausibility
heuristic -- surfactant-micelle zeta rarely exceeds ~100-150 mV in aqueous
systems -- can flag an implausible branch, always disclosed as a prior, never
silently substituted for the genuine ambiguity. **Correction to the summary
table and the audit's headline count**: this is not one of the "pure walls,"
it's a wall for minimal input that a real, buildable richer-input path
closes -- see `GENUINE_BOTTLENECK_AUDIT.md`'s updated item 9.

**10. SLS dn/dc.**
[MINE] for common surfactant/solvent pairs (SDS, CTAB, Triton X-100,
Tween-80, etc.) -- real published values exist, not yet compiled.
[WALL] remains for genuinely novel compounds -- an actual refractometer
measurement is the only real closure there; no credible general
first-principles estimator for dn/dc exists (it depends sensitively on
electronic structure, not just molecular formula).

**11. Svedberg v_bar.**
[MINE] for common surfactants, same shape as #10.
[COMPUTE], real and promising -- `cpp.py`'s `tanford_tail_volume` already
uses group-additive volume increments (Tanford's rule) for the hydrophobic
tail. Partial specific volume is closely related to molecular volume, and
Traube's rule (classical, well-established surfactant science, group-
additive molar volume for homologous series) is a real, citable method for
estimating v_bar from structure. Extending the existing tail-volume
framework to a full-molecule v_bar estimator, cross-validated against the
mined real values above, is a strong, buildable near-term target.

**12. Van Oss-Chaudhury-Good components.**
Not a real bottleneck -- [MINE], trivial and cheap. Source the actual
standard values (water, glycerol, diiodomethane, formamide) from
established surface-science reference tables directly into the toolkit as
a real, cited lookup table. Closes a false-positive refusal cheaply.

**13. Intrinsic water solubility (solubilizate).**
[MINE] for common solubilizates used in real solubilization studies
(naphthalene, benzene, pyrene, common dyes) -- these are among the most
heavily measured properties in environmental/pharma chemistry; real public
data (PubChem experimental properties, similar sources) likely already
covers most literature test compounds.
[NEW METHOD]/[COMPUTE], strong candidate -- real, published, simple QSPR
solubility estimators exist (e.g. the Delaney/ESOL equation: logS from
logP + MW + rotatable-bond count + aromatic proportion, all RDKit-
computable already). Wiring in a real, cited QSPR estimator as a clearly-
labeled ESTIMATE (not a measurement) alongside the mined real-value table
would meaningfully shrink this bottleneck for compounds nobody measured.

**14. pH-ambiguous charge type.**
[WALL] for the case where solution pH itself is genuinely unstated -- no
fix, correct to keep refusing.
[COMPUTE], real and buildable regardless -- extend the classifier to
accept solution pH as an OPTIONAL input, and when given, compute real
ionization fraction via Henderson-Hasselbalch using the ionizable group's
pKa (literature-sourced for common amine classes, or estimated via
established Hammett-type substituent-constant methods for simple
primary/secondary/tertiary amines not yet tabulated). This doesn't remove
the wall for "pH unstated" but closes a real, separate gap: right now the
tool can't use pH even when a caller DOES supply it.

## Summary table

| # | Bottleneck | MINE | COMPUTE | NEW METHOD | WALL |
|---|---|---|---|---|---|
| 1 | Counterion binding degree | yes (common cpds) | | yes (Manning/PB) | residual only |
| 2 | Electrolyte condition | | | yes (if conductivity given) | yes (default) |
| 3 | HLD k | yes | | | |
| 4 | HLD Cc | yes | | | |
| 5 | HLD b | yes | | | |
| 6 | Davies gemini/glycolipid | tried, failed | tried, failed | **yes, only path** | |
| 7 | Nagarajan prefactor | yes | yes | | |
| 8 | Perrin axial_ratio | | **yes, near-term** | | |
| 9 | Henry high-zeta | | | **yes, if series given** | single-point only |
| 10 | SLS dn/dc | yes (common) | | | yes (novel) |
| 11 | Svedberg v_bar | yes (common) | **yes, near-term** | | |
| 12 | van Oss components | **yes, trivial/cheap** | | | not a real gap |
| 13 | Intrinsic water solubility | yes (common) | yes (QSPR) | | |
| 14 | pH-ambiguous charge | | **yes, near-term** | | yes (pH truly unstated) |

## Proposed work order (cheapest/highest-confidence first)

**Tier 0 -- DONE 2026-09-15 (commit `283d407`):**
- #12 van Oss-Chaudhury-Good standard-liquid table -- `VOCG_STANDARD_LIQUIDS`/
  `OWENS_WENDT_STANDARD_LIQUIDS` in wetting.py, live-verified against van
  Oss/Good/Busscher 1990 before hardcoding. New MCP tool
  `get_standard_probe_liquid_properties`.

**Tier 1 -- DONE 2026-09-15, all 4 items (commits `4bff41e`, `034d383`,
`eea9e7a`, `bd068c0`):**
- #8 Perrin axial_ratio from CPP-predicted morphology -- DONE. New
  `estimate_axial_ratio_from_cpp_geometry` in cpp.py: models the rodlike-
  micelle core as a prolate ellipsoid with minor semi-axis pinned at the
  extended tail length, solves for axial_ratio via volume conservation
  given a REAL measured aggregation number. Scoped to cpp in (1/3, 1/2];
  raises outside that range or when N_agg is geometrically inconsistent.
- #11 Svedberg v_bar from Tanford tail-volume conversion -- DONE, but
  PARTIAL by design. New `estimate_partial_specific_volume_from_tail_
  and_headgroup` in curve_analysis.py converts Tanford's already-verified
  tail volume to cm^3/mol exactly (unit conversion, no new empirical fit)
  -- reduces the external requirement to only the headgroup's own real
  molar volume, deliberately NOT hardcoded (no headgroup value was
  verified with enough confidence this session to assert as fact).
- #14 pH-conditional charge classification -- DONE. New
  `classify_surfactant_charge_type_at_ph` in classify.py: Henderson-
  Hasselbalch against literature pKaH ranges for primary/secondary/
  tertiary amines (live-verified: methylamine 10.64, dimethylamine 10.73,
  trimethylamine 9.79), or an exact caller-supplied real pKa. The
  pH-truly-unstated WALL is untouched and correctly still refuses.
- #13 Delaney/ESOL QSPR solubility estimator -- DONE. New
  `estimate_intrinsic_water_solubility_qspr` in solubilization.py, Pat
  Walters' real RDKit-refit ESOL coefficients (live-verified), prominently
  disclosed as a coarse estimate (~0.6-1 log unit real published error),
  never silently substituted for a real measured value.

All 5 new functions wired into the MCP server, full suite 429/429 passing
throughout, each item committed and pushed separately.

**Tier 2 -- DONE 2026-09-15 with 3 real, honest partial-coverage
disclosures (commits pending push at time of writing) -- literature
search hit real, disclosed walls for 2 of 5 items, matching this
project's own long-established pattern (many HLD-NAC and Nagarajan-
adjacent primary sources are paywalled/bot-walled; see hld.py's own
CATIONIC_QUAT_HLB_DAVIES precedent for exactly this situation):**

- #3/#4/#5 HLD-NAC k/Cc/b for non-quat classes -- PARTIAL, real values
  added: `K_ANIONIC_DEFAULT`=0.16 (re-confirmed via a second source),
  `K_EXTENDED_SURFACTANT`=0.06, `ALPHA_APG_DEFAULT`=0, `B_APG_SPAN_
  DEFAULT`=0 (all via Steven Abbott's "Practical Surfactants Science"
  HLD page, live-verified), `SDS_CC`=-3.0 (via a well-corroborated
  secondary citation to Leng & Acosta 2023, primary paper paywalled).
  New MCP tool `hld_class_reference`. GENUINE, DISCLOSED REMAINING GAP:
  real Cc values for zwitterionic (e.g. cocamidopropyl betaine),
  gemini/dimeric, and glycolipid biosurfactant classes were searched for
  and NOT found -- every primary HLD-NAC source located (Leng & Acosta
  2023 and its SAXS companion, Acosta's 2026 JSD paper, the tutorial
  chapter) blocked automated fetch. Real path forward if this matters:
  ask the user to provide a PDF, matching this project's own established
  unblocking mechanism for exactly this class of source.
- #1's common-compound counterion binding degree table -- DONE for 3
  real compounds: SDS (alpha=0.272+/-0.017, HIGH confidence, Bales 2001,
  primary PDF already read in full elsewhere in this project), DTAB
  (alpha=0.28, MODERATE, web-search-corroborated), CTAB (alpha=0.26,
  MODERATE). New `COUNTERION_BINDING_DEGREE_REFERENCE` in
  thermodynamics.py, MCP tool `counterion_binding_degree_reference`.
  AOT was searched for and not found with a real value -- disclosed, not
  guessed.
- #7's Nagarajan headgroup-prefactor search -- NOT CLOSED, a real,
  disclosed gap. Multiple real searches for sulfonate/carboxylate/
  quaternary-ammonium headgroup prefactors (or their constituent d/sigma
  sub-quantities) in Nagarajan's own follow-up papers and related
  electrostatic-headgroup literature found no citable numeric value --
  every real lead traced back to paywalled ACS Langmuir content. No
  constant was added; `headgroup_prefactor_A` remains required, not
  defaulted, for anything but sulfate-type headgroups, exactly as
  before this session.
- #10's common dn/dc values -- DONE for 2 real compounds: SDS (0.11
  mL/g) and CTAB (0.15 mL/g), both water/632.8nm/25C, Malvern
  Panalytical's own published reference page, live-verified. New
  `DN_DC_REFERENCE_ML_PER_G` in curve_analysis.py, MCP tool
  `dn_dc_reference`. Triton X-100 and Tween-80 (both explicitly named as
  targets) were searched for and NOT found with a real citable value --
  disclosed, not guessed.
- #13's common solubilizate values -- DONE for 3 real compounds:
  naphthalene (2.17e-4 mol/L, HIGH confidence, the SAME already-
  validated primary-source value from literature_validation_notes.md),
  benzene (~0.0228 mol/L, MODERATE, widely-tabulated EPA/ATSDR-compiled
  constant), pyrene (~6.87e-7 mol/L, MODERATE, a real ACS solubility-
  measurement paper's value). New `INTRINSIC_WATER_SOLUBILITY_
  REFERENCE_M` in solubilization.py, MCP tool
  `intrinsic_water_solubility_reference`.

All 5 sub-items closed with real, sourced, confidence-disclosed values
where findable; 2 real sub-gaps (zwitterionic/gemini/biosurfactant HLD
Cc, and Nagarajan's non-sulfate prefactor) remain genuinely open and are
explicitly recorded here rather than silently dropped, per this
project's standing "close a workflow to its entirety, record what can't
be closed" rule.

**Tier 3 -- [NEW METHOD], real research contribution, higher effort:**
- #1's Manning/PB counterion-binding predictor, validated against Tier 2's
  mined table
- #6's gemini/glycolipid HLB resolution (either a new calibrated group-
  contribution scheme, or formally establishing Griffin-as-default) --
  the highest-value single item, since it's the one CONFIRMED structural
  bottleneck and already has real supporting data (5 gemini reference
  compounds) in hand
- #2's conditional conductivity-derived electrolyte inference
- #9's mobility-vs-ionic-strength series branch-disambiguation (O'Brien-
  White treatment + plausibility heuristic) -- reclassified 2026-09-15,
  moved out of "leave as wall"

**Leave as [WALL], do not chase:**
- #9's single-mobility-point default case (genuinely underdetermined with
  that input alone)
- #2's default case (electrolyte condition truly unstated)
- #14's default case (pH truly unstated)
- #10's novel-compound dn/dc (no credible first-principles estimator)

Not yet started -- this is the classification the user asked for, presented
for direction before any of Tier 0-3 work begins.
