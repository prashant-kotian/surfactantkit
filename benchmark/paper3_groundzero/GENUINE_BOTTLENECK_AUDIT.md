# Genuine bottlenecks: where even the literature has no data

Purpose: a complete audit of every place across SurfactantKit's 61 MCP tools
where a real computation is gated on information that genuinely doesn't
exist in the literature for most compounds -- not a missing MCP tool, not a
solvable "if only we built X" gap, but a real, disclosed epistemic wall this
project's own "never guess" discipline already refuses to cross. Compiled
2026-09-15 by grepping every module's own docstrings for this project's
standing gap-disclosure language ("do not guess," "no default," "must be
independently measured," etc.) and reading each hit's real context.

14 distinct bottlenecks found, in 5 categories by KIND (not all "missing
data" is the same kind of problem -- distinguishing this matters for what,
if anything, could ever close each one).

## Category A: counterion/ionic-character facts (the most pervasive family)

**1. Counterion binding degree / dissociation alpha (beta).** Needed for
deltaG_mic on any ionic surfactant (`gibbs_free_energy_micellization`,
threaded through `orchestrate.py`'s whole pipeline, and through any van't
Hoff-derived deltaG). Real, independent measurement required (conductivity
slope-ratio method, EPR, or a literature value for that SPECIFIC compound)
-- genuinely absent for most real surfactants unless someone measured it.
This is the single most-tested bottleneck in this project's benchmark work
(GZ-01, GZ-05, Category D's entire design).

**2. Electrolyte condition (excess salt vs. none) for the Gibbs adsorption
prefactor n.** Appears in `gibbs_gamma_max`, `szyszkowski_surface_tension`,
`szyszkowski_fit_K`, `frumkin_fit_K_and_a` -- all four raise/refuse without
it. A real fact about the SOLUTION being characterized, not the compound --
no amount of SMILES parsing or curve-fitting recovers it; it must be stated
or refused. (Real nuance found this session, GZ-11: this ONLY blocks the
Gibbs-prefactor-dependent thermodynamics path -- the geometry/CPP path via
Nagarajan's method does NOT need it, since it derives screening length from
the CMC itself.)

## Category B: mixture-theory calibration constants (HLD-NAC)

**3. HLD's k (EACN-scaling constant).** Real, sourced values exist for only
two classes: `CATIONIC_QUAT_K=0.7` and a general anionic default (~0.16, via
the source paper's own citation of Abbott 2017). For any other surfactant
class -- nonionic beyond the one worked example, zwitterionic, gemini,
biosurfactant -- this needs real HLD-titration data specific to that class,
which was not found/sourced in this project's literature mining so far.

**4. HLD's Cc (characteristic curvature).** Same shape of gap: only
`CATIONIC_QUAT_CC` is a real sourced constant. System-specific, needs a real
HLD-titration/salinity-scan for that exact compound -- genuinely absent for
most real surfactants.

**5. HLD's b (nonionic salinity-scaling constant).** Only one real value
sourced (~0.13 dL/g, alkyl ethoxylates). Unsourced for other nonionic
structural classes.

## Category C: HLB group-contribution gaps

**6. Davies group numbers for untabulated headgroups.** Originally looked
like a "missing number" gap, cheaply fixable by extending
`derive_davies_group_number_from_griffin`'s already-proven cross-calibration
method. **Tested twice against real data this session (see
LLM_REFUSAL_RESEARCH.md) and confirmed NOT fixable that way** -- for a
gemini's two headgroups, neither naive doubling (unphysical HLB=32.95,
outside the 0-20 scale) nor a single combined group number (36% spread
across 5 real reference compounds, driven monotonically by spacer length)
converges. This is the one bottleneck in this whole audit that is
**structural, not just data-absent**: Davies' single-headgroup linear model
does not generalize to a bis-headgroup architecture, and no literature value
could ever fix that, because the problem isn't a missing number, it's that
the model's own assumptions don't hold. The glycolipid case (ring-methyl
partition ambiguity) is a genuine, separate, not-yet-fully-tested gap in the
same category.

## Category D: geometry/electrostatics/dynamics-specific physical constants

**7. Nagarajan's headgroup_prefactor_A.** The source paper (Nagarajan 2002)
only gives a real, validated COMBINED value (82.0, itself a product of two
separate real physical quantities -- a charge-separation distance and a
bare hydrocarbon-water interfacial free energy) for ONE headgroup class:
sodium-alkyl-sulfate-type. For any other ionic headgroup (sulfonate,
phosphate ester, carboxylate, quaternary ammonium), this real constant has
not been sourced -- using 82.0 for those would be exactly the kind of
"treat a system-specific number as universal" error this project's own
discipline exists to prevent.

**8. Perrin shape-correction's axial_ratio.** `hydrodynamic_radius_perrin_
corrected` needs a real, independently-known aggregate aspect ratio (from
CPP-predicted morphology, or a separate SANS/SAXS/cryo-TEM measurement) --
explicitly NOT inferable from the diffusion coefficient alone (the inverse
problem is underdetermined without an independent volume estimate).

**9. Henry-function mobility inversion at high zeta -- a DIFFERENT kind of
gap, not missing data, and NARROWER than first stated (revised 2026-09-15
after the user pushed on this directly).** `zeta_potential_relaxation_
corrected` can genuinely fail to find a unique zeta from a SINGLE mobility
value at a SINGLE ionic strength, because the real mobility-vs-zeta
relationship folds back (non-monotonic) above roughly 100 mV -- that part
is a real, physical non-uniqueness, not missing data, and no amount of
computation on that one number resolves it. **But it is not a pure wall**:
feed the tool a better dataset and it stops being underdetermined. Two real
paths, neither requiring a new measurement TYPE: (a) a mobility-vs-ionic-
strength SERIES instead of one point -- kappa*a changes with ionic
strength, so the shape of that series (plateau vs. a maximum-and-fold-back
within the measured range) is itself diagnostic of which branch the system
sits on, standard colloid-science practice, not a new technique; (b) the
fuller O'Brien-White (1978) numerical relaxation treatment in place of the
simpler Henry approximation, which explicitly locates the mobility maximum
from kappa*a and surface conductivity using only electrophoretic data. A
disclosed, literature-grounded heuristic (surfactant-micelle zeta rarely
exceeds ~100-150 mV in aqueous systems) can additionally flag implausible
branches -- always disclosed as a plausibility prior, never a silent pick.
See `BOTTLENECK_RESOLUTION_PLAN.md` item 9 for the reclassification.

## Category E: light-scattering / sedimentation instrument constants

**10. SLS's dn/dc (refractive index increment).** Real, solvent-wavelength-
system-specific constant needed for `aggregation_number_from_sls_debye_
plot` -- genuinely varies even within one macromolecule class (the
function's own docstring cites Zhao, Brown & Schuck 2011 on exactly this
variability), must be independently measured or cited per real system.

**11. Svedberg's v_bar (partial specific volume).** Same shape of gap for
`aggregation_number_from_svedberg_equation` -- real, system-specific,
must be independently known, not derivable from the compound's formula
alone.

## Category F: a case that LOOKS like this category but genuinely isn't

**12. Van Oss-Chaudhury-Good's per-liquid LW/acid/base surface-tension
components.** `wetting_van_oss_chaudhury_good` needs real values for its 3
test liquids -- but unlike everything else in this audit, these ARE
genuinely tabulated in standard surface-science literature for the handful
of liquids actually used for this method (water, glycerol, diiodomethane,
formamide, etc.). This is not a "literature has no data" bottleneck at all
for the standard case -- it would only become one for a genuinely uncommon
test liquid. Included here specifically to mark the boundary: not every
"caller must supply this" parameter in this codebase is a real literature
gap, and this project's own discipline (require it explicitly, never
default it) doesn't mean the discipline believes it's unknowable.

## Category G: solubilization

**13. Intrinsic water solubility (of the SOLUBILIZATE, not the
surfactant).** `molar_solubilization_ratio` / `micelle_water_partition_
coefficient` need this real, compound-specific measured value for whatever
is being solubilized -- genuinely absent for most solute/surfactant
combinations unless a specific uptake study measured it (this project's own
literature-validation notes document two real, named dead-end search
attempts before finding one real usable case, Paria & Yuet 2006).

## Category H: structural ambiguity that traces back to a real, unstated fact

**14. pH-ambiguous charge type (free amine, no stated solution pH).**
`classify_surfactant_charge_type` correctly refuses to resolve this rather
than defaulting to cationic. Different in kind from the others: the
compound's own pKa may well be literature-derivable in principle, but the
ACTUAL solution pH being characterized is a real, experiment-specific fact
no structure-only computation can ever recover -- structurally the same
category as electrolyte condition (item 2), a fact about the experiment,
not the compound.

## Summary table

| # | Bottleneck | Kind | Any path to closing it? |
|---|---|---|---|
| 1 | Counterion binding degree | missing measurement | only per-compound, from real data |
| 2 | Electrolyte condition | missing experimental fact | only if told |
| 3 | HLD k (non-quat classes) | missing calibration | real HLD-titration data, not yet sourced |
| 4 | HLD Cc (non-quat classes) | missing calibration | same |
| 5 | HLD b (non-ethoxylate nonionics) | missing calibration | same |
| 6 | Davies gemini/glycolipid groups | **structural, not data** | **no** -- confirmed via direct testing |
| 7 | Nagarajan prefactor (non-sulfate heads) | missing calibration | real per-headgroup-class literature source |
| 8 | Perrin axial_ratio | missing independent measurement | real SANS/SAXS/cryo-TEM or CPP-inferred |
| 9 | Henry mobility inversion, high zeta | **physical non-uniqueness, single point only** | **partial** -- a mobility-vs-ionic-strength series or O'Brien-White breaks it; single point does not |
| 10 | SLS dn/dc | missing instrument constant | real per-system measurement |
| 11 | Svedberg v_bar | missing instrument constant | real per-system measurement |
| 12 | van Oss LW/acid/base | **not actually a gap** | already tabulated for standard liquids |
| 13 | Intrinsic water solubility | missing measurement | real per-solute uptake study |
| 14 | pH-dependent charge | missing experimental fact | only if told |

**11 of 14 are real, standing epistemic walls in their default (minimal-
data) form** -- the toolkit's correct, durable contribution there is
exactly what it already does: refuse rather than guess. **1 (item 6) is
confirmed structural**, not fixable by any amount of future literature
mining. **1 (item 12) is not actually a bottleneck** once sourced correctly
for standard liquids. **1 (item 9) is a wall only for a single-point input**
-- richer input (a real series, not just more literature) breaks it, a
distinction worth building rather than assuming away (see
`BOTTLENECK_RESOLUTION_PLAN.md`). None of these are "go build an MCP tool"
opportunities in the naive sense -- they're refusal points to keep testing
and keep reinforcing UNLESS a real richer-dataset or first-principles path
exists, which is exactly what the resolution plan sorts out.
