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
disclosures (commit `5f8e20c`) -- literature search hit real, disclosed
walls for 2 of 5 items, matching this project's own long-established
pattern (many HLD-NAC and Nagarajan-
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

**2026-09-15, same day -- most of the above genuine gaps CLOSED with real
primary-source PDFs the user provided directly** (per the user's own
instruction after the Tier 2 report: "when you dont find appropriate
paper, dont flag the tool as not fixed, instaed hold it and give me list
papers you need... ill do both and give you"). 17 real papers provided,
identified and renamed to `Author_Year_Topic.pdf` (see
`~/.claude/projects/.../memory/feedback_papers_folder_rename_convention.md`),
13 directly relevant:

- **#7 Nagarajan headgroup prefactor -- SUBSTANTIALLY CLOSED**, not via
  a lookup, but via the ACTUAL primary source: Nagarajan & Ruckenstein,
  *Langmuir* 7 (1991) 2934-2969 ("Theory of Surfactant Self-Assembly: A
  Predictive Molecular Thermodynamic Approach") -- the real, earlier,
  more general founding theory the 2002 paper's `headgroup_prefactor_A`
  shortcut was itself built from. Its own Table I gives real molecular
  constants (a_p, a_o, delta, d) for sodium sulfate, sodium sulfonate
  (a genuinely NEW headgroup class), and N-betaine (ZWITTERIONIC,
  another genuinely new class), plus the paper's own complete ionic
  (eqs. 70-73), dipole (eqs. 67-68), and steric (eq. 66) headgroup free-
  energy formulas -- all implemented as 3 new real functions in cpp.py
  (`nagarajan_ruckenstein_ionic_headgroup_free_energy`,
  `_dipole_headgroup_free_energy`, `_steric_headgroup_free_energy`) plus
  the real constants table (`NAGARAJAN_RUCKENSTEIN_HEADGROUP_
  CONSTANTS`). Validated three ways: a real independent physical
  constant (the Bjerrum length in water, ~7.0-7.1 A, that both new
  functions' CGS unit handling must reduce to exactly -- caught a real,
  confusing unit-conversion mistake in the first draft, rewritten
  cleanly before shipping); a direct cross-check against this module's
  own already-paper-verified `nagarajan_debye_huckel_kappa_inverse`;
  and real physical-limit/sign checks. Cationic quaternary ammonium and
  carboxylate headgroups are STILL not in this paper's own table either
  -- a real, disclosed, narrower remaining gap than before, not fully
  closed.
- **HLD Cc gap -- SUBSTANTIALLY EXTENDED**: the real primary PDF (Leng &
  Acosta, *J. Surfactants Deterg.* 26(3) (2023) 287-301, both its
  published and ChemRxiv-preprint form) gives real, literature-cross-
  checked Cc values for SDS, SDHS, SLES, C10PO4S, C16DPODS, **AOT**
  (this project's single most-used compound), BCl, and DPCl -- 8 real
  compounds, replacing the earlier single unverified SDS_CC=-3.0
  placeholder with a real, dual-sourced (this-work + independent
  literature) value for SDS (-2.63 vs. literature -2.5) and 7 more. New
  `CC_REFERENCE_ANIONIC_CATIONIC` dict + MCP tool
  `hld_characteristic_curvature_reference` in hld.py. Also closed: a
  real BIOSURFACTANT Cc (-1.41, rhamnolipid, via Nguyen & Sabatini 2008
  as cited in Hellweg, Oberdisse & Sottmann, *Front. Soft Matter*
  3:1260211 (2023), both read in full) -- `RHAMNOLIPID_CC`. A real,
  disclosed DISCREPANCY was found and documented, not silently resolved:
  the Hellweg review states Cc(AOT)=-0.92, contradicting Leng & Acosta's
  own directly-read value (~2.4-3.5, matching independent literature
  2.5) -- -0.92 is suspiciously exactly SDHS's own value from the SAME
  paper, almost certainly a citation error in the review, not a second
  real AOT measurement; the directly-read primary value is what's
  shipped, with the discrepancy documented in both the code and here.
  ZWITTERIONIC and GEMINI Cc specifically remain genuinely open -- none
  of the 17 papers studied either class's own Cc directly (though the
  zwitterionic case now has a real, arguably better, mechanistic
  alternative: the dipole free-energy function above, using real
  N-betaine constants, rather than a bare Cc number).
- **#1 AOT counterion binding degree -- CLOSED**: Thapa, Ray, Dey,
  Sultana, Aswal & Ismail, *RSC Adv.* 5 (2015) 45956-45964, Table 1
  (Corrin-Harkins plot, NaCl medium, below the real critical salt
  concentration c*) gives a real beta=0.39 (alpha=0.61), HIGH
  confidence, added to `COUNTERION_BINDING_DEGREE_REFERENCE`.
- **#10 Triton X-100 dn/dc -- CLOSED**: Stubicar, Matejas, Zipper &
  Wilfing, in Mittal (ed.), *Surfactants in Solution*, Plenum Press
  (1989) 181-193, states it directly in their own instrument-
  calibration section: dn/dc=0.140+/-0.005 mL/g (546 nm, 20 C), stable
  across water and KCl/KBr/KI electrolyte solutions alike -- added to
  `DN_DC_REFERENCE_ML_PER_G`. Tween-80 remains genuinely open -- not in
  this batch either.

Remaining genuinely open after this batch: zwitterionic/gemini HLD Cc
(a real Cc number specifically); Nagarajan cationic-quat/carboxylate
prefactor; Tween-80 dn/dc. All explicitly recorded, not silently
dropped. Full suite: 486/486 passing after this round (one real bug
found and fixed along the way: a stale test assertion still expecting
the old SDS_CC=-3.0 placeholder after the real value was substituted).

**2026-09-15, same day -- SECOND real-paper batch, closing essentially
everything still listed as open just above.** User asked "give me
paper/book keyword to search for all open issues"; a grouped list
(Groups A-F) was given, and the user provided 10 more real PDFs into a
`New folder` subdirectory, immediately renamed per the established
convention. Every remaining item in this section was closed:

- **Nagarajan cationic-quat/carboxylate headgroup prefactor -- CLOSED**:
  Nagarajan's own later, more comprehensive review chapter, "Theory of
  Micelle Formation," Ch. 1 in *Structure-Performance Relationships in
  Surfactants*, 2nd ed., Taylor & Francis (2003), Table 3 ("Molecular
  Constants for Surfactant Headgroups") -- the user pointed directly at
  the exact table via the filename itself. This independent, later
  source reports IDENTICAL values for every headgroup already sourced
  from the 1991 paper (a real cross-confirmation, not just one paper's
  self-consistency), AND adds real constants for trimethyl ammonium
  bromide and pyridinium bromide (CATIONIC, both real gaps now closed),
  sodium and potassium carboxylate (the CARBOXYLATE gap now closed),
  four more nonionic classes, and a SECOND real zwitterionic class
  (lecithin -- genuinely needing BOTH the ionic delta AND the dipole d,
  a real structural complexity the table discloses rather than
  resolves one way). `NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS` in
  cpp.py grew from 4 to 13 real, sourced entries.
- **Zwitterionic HLD Cc -- CLOSED**, comprehensively: Acosta, Harwell &
  Sabatini (eds.), *Surfactant Formulation Engineering Using HLD and
  NAC*, Elsevier (2021), Table 1.1 -- Acosta's own definitive compiled
  database. Real Cc for FIVE zwitterionic surfactants: lecithin (5.5,
  matching Nouraei & Acosta 2017 -- also independently provided this
  round and read in full -- EXACTLY, a genuine cross-confirmation
  within the same real dataset), Epikuron 200 (5.1), C4-mPC (3.0),
  dodecylsulfobetaine/lauryl sultaine (-0.7), lauramine oxide (-4.0),
  and **cocamidopropyl betaine (CAPB)** -- the exact compound
  repeatedly named as the target example throughout this whole
  investigation -- with two real estimates (-5.2 and -2.1). New
  `CC_REFERENCE_ZWITTERIONIC` in hld.py. The AOT=-0.92 discrepancy
  flagged in the previous batch is now further confirmed as a citation
  error: this book's own text states AOT's real bi~0.32, consistent
  with the strongly positive Cc already sourced, nowhere near -0.92.
- **Gemini HLD Cc -- CLOSED**: same Table 1.1, a real gemini/dimeric
  surfactant ("Gemini benzene sulfonate C16 [Ph](SO3Na)O[Ph](SO3Na)"),
  k=0.16 (matching `K_ANIONIC_DEFAULT` exactly -- another real cross-
  consistency check), Cc=-7.4 (this work), cross-checked against an
  independent comparison value of -6.6. New `GEMINI_BENZENE_SULFONATE_
  CC`/`_LITERATURE` in hld.py.
- **Tween-80 dn/dc -- CLOSED**: Jiang, F., PhD dissertation, Virginia
  Tech (2011), "Effects of the Non-ionic Surfactant Tween 80 on the
  Enzymatic Hydrolysis of Model Cellulose and Lignocellulosic
  Substrates" -- its own SPR methods section states directly: dn/dc =
  0.132 mL/g (690 nm, 25 C, Wyatt Optilab rEX). Added to
  `DN_DC_REFERENCE_ML_PER_G`.
- **Real, identified resources for Tier 3 item 1 (counterion-binding
  predictor), NOT yet implemented -- deliberately deferred, not a
  smaller task than it looks.** The user also provided the REAL Manning
  1969 paper (*J. Chem. Phys.* 51, 924) AND, more importantly, a full
  3-paper Blankschtein-group series specifically calibrated for
  SURFACTANT micelles (not generic polyelectrolytes like Manning's own
  theory): Srinivasan & Blankschtein, *Langmuir* 19 (2003) 9932-9945 and
  9946-9961 (molecular-thermodynamic theory of counterion binding,
  validated against real SDS+NaCl, alkali dodecyl sulfates, multivalent
  counterions, and organic counterions), and Goldpise & Blankschtein,
  *Langmuir* 21 (2005) 9850-9865 (extending it to ionic-nonionic and
  ionic-zwitterionic MIXTURES). This is a real, better-targeted
  alternative to plain Manning theory for this specific goal -- but
  implementing the full free-energy-minimization model is a
  substantially larger undertaking than a lookup-table extension, on
  the same scale as reimplementing a piece of Nagarajan's own
  micellization theory. Deliberately scoped OUT of this session's work;
  flagged here as the real, identified starting point for whenever Tier
  3 item 1 is actually taken on.
- **Real, supporting (not load-bearing) evidence for Tier 3 item 6**
  (gemini/glycolipid Davies HLB): Liao et al., *Arabian J. Chem.* 16
  (2023) 105111, reports real Griffin-method HLB values (6.22-7.97) for
  four real NONIONIC gemini surfactants, cross-validated against their
  own observed W/O vs. O/W emulsion behavior -- real, published,
  recent (2023) confirmation that Griffin's method is already being
  used successfully as the practical default for gemini surfactants,
  supporting this project's own proposed resolution path (formally
  establish Griffin-as-default rather than extending Davies). No new
  code needed -- `hlb_griffin` already exists; this is documentation-
  grade supporting evidence, not a capability gap.

**Genuinely still open after this second batch**: only the Blankschtein
counterion-binding theory's actual implementation (Tier 3 item 1,
scoped as real future work, not blocked on sourcing) and the Davies-
gemini-HLB resolution's own code decision (Tier 3 item 6, now with even
more real supporting data in hand). Full suite: 494/494 passing.

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
