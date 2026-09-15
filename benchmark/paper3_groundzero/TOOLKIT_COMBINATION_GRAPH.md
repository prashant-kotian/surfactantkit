# SurfactantKit ecosystem: full technique combination graph

Purpose: before designing the next batch of ground-zero questions, map every
real primary -> secondary -> tertiary chain across BOTH repos (SurfactantKit's
61 MCP tools + surfactantkit-reports' 10 report pipelines and their property
cascades), so new questions can chain real, already-implemented dependencies
instead of inventing ad hoc combinations. Built 2026-09-15 by reading every
cascade_*.py file in surfactantkit-reports (previously unexplored this
session) plus TOOLKIT_CAPABILITY_AUDIT.md's earlier module-by-module read of
SurfactantKit core.

## Entry points (primary, from raw data alone)

Each of these takes ONLY raw experimental data (+ SMILES where relevant) and
needs no other tool's output first:

1. `cmc_from_surface_tension_curve` (raw C-vs-gamma tensiometry)
2. `cmc_from_conductivity_curve` (raw C-vs-conductivity)
3. `aggregation_number_from_quenching_curve` (SSFQ)
4. `aggregation_number_from_sls_debye_plot` (SLS)
5. `aggregation_number_from_dls`, `_from_svedberg_equation` (2 more real techniques, not yet wired into a report cascade)
6. `rubingh_beta_regression` (raw binary composition series)
7. `rodenas_x1_series` (raw binary composition series, model-independent)
8. `select_mixture_model` (NEW, this session -- raw binary composition series, no model named)
9. `eommm_global_fit` (raw binary composition series)
10. `fit_k_and_cc_from_salinity_scan` (HLD-NAC, raw EACN-vs-S* salinity scan)
11. `vant_hoff_multi_point_fit` (raw CMC-vs-T, 5+ points)
12. `szyszkowski_fit_K` (raw C-vs-gamma below CMC)
13. `classify_surfactant_charge_type`, `classify_surfactant_structural_family` (SMILES alone)
14. `digitize_plot` (NEW, this session -- an IMAGE, feeds any of the above as its raw-data source)

## Chains from each entry point (what a real researcher can derive next)

**From CMC-from-tensiometry** (richest chain, `cascade_cmc_surface_tension.py`):
  base (CMC, gamma@CMC, gamma0, premicellar slope, C20/pC20 if the curve reaches it)
  -> **thermodynamics branch** (needs: system_type [never guessed], T, counterion_factor):
     Gamma_max -> A_min -> X_cmc -> deltaG_mic (pseudo-phase)
  -> **geometry branch, Nagarajan's method** (needs: n_carbons, a real headgroup-specific
     prefactor constant, T, dielectric constant -- valid ONLY ionic/no-added-salt):
     tail volume, critical length (Tanford, universal constants)
     -> kappa^-1 (Debye-Huckel, using the CMC ITSELF as the ionic-strength proxy --
        genuinely different from the "given electrolyte concentration" pattern GZ-06 used)
     -> equilibrium area a_e (chain-length/screening-DEPENDENT, the real alternative
        to a fixed Amin -- THIS is the real tool GZ-10 should have chained into
        instead of declaring "not computable here")
     -> CPP -> predicted morphology
  -> **solubilization branch** (needs a real, separate solubilizate-uptake experiment's
     3 concentrations): MSR, Km

**From CMC-from-conductivity** (`cascade_cmc_conductivity.py`): same thermodynamics +
  Nagarajan-geometry + solubilization branches as tensiometry, PLUS the technique's own
  real counterion-binding-degree beta (from the two slopes) is already in hand for free --
  meaning deltaG_mic's counterion_factor doesn't need to be separately supplied/guessed
  here the way it does for a plain tensiometry-only run. **This is a real, clean
  "conductivity gives you something tensiometry can't" entangled question shape.**

**From SSFQ/SLS aggregation number** (`cascade_aggregation_number.py`):
  base (N_agg) -> micellization degree (free, both definitions) 
  -> mass-action free energy (NONIONIC ONLY, needs T + a CMC from a SEPARATE run)

**From Rubingh regression** (`cascade_rubingh.py`):
  base (beta_mean, beta_std, x1 per point) -> activity coefficients + excess free energy
  at each point (free) -> **Maeda free energy** (needs the run's own cmc1/cmc2 + T --
  a real chain FROM Rubingh's beta INTO a different thermodynamic model, not a
  competing cmc_mix predictor, matching why select_mixture_model correctly excludes
  Maeda from its own BIC race but this IS a legitimate follow-on question)

**From Rodenas series** (`cascade_rodenas.py`): base (x1 per point, model-independent)
  -> activity coefficients (needs the two pure CMCs, new inputs)

**From EOMMM fit** (`cascade_eommm.py`): base (W12, W21, x1/cmc_var per point)
  -> activity coefficients at ANY queried x1 (genuinely PREDICTIVE -- x1_query need not
  be one of the originally measured compositions, a real "use the fitted model to
  extrapolate/interpolate" question shape not yet used in GZ-01..10)

**From HLD-NAC salinity scan** (`cascade_hld.py`): base (k, Cc)
  -> predict HLD for a NEW oil/salinity/T this scan never measured
  -> predict optimal (Winsor III) salinity for a NEW oil (EACN) never measured
  -> HLB estimate (CATIONIC QUATERNARY AMMONIUM ONLY -- a real, narrow applicability
     condition, a genuine "does the model know this branch doesn't apply to its own
     surfactant's charge type" trap if paired with a non-quat cationic or any other charge type)

**From van't Hoff multi-point fit** (`cascade_vant_hoff.py`): base (deltaH, deltaS, deltaCp)
  -> thermodynamic character classification (enthalpy- vs entropy-driven, free --
     a real qualitative judgment call from already-computed numbers)

**From Szyszkowski K-fit** (`cascade_szyszkowski.py`): base (K, Gamma_max)
  -> predicted surface tension at NEW concentrations (genuinely predictive)
  -> Frumkin extension (needs a real, independently-known lateral-interaction
     parameter -- NOT fit from the same data, a genuine "this needs something the
     researcher must supply separately" gap, different in kind from select_isotherm_model's
     own internal Langmuir-vs-Frumkin BIC choice)

## Structural/charge classification as a MULTIPLIER, not a standalone question type

`classify_surfactant_charge_type` and `classify_surfactant_structural_family` were used
in GZ-01..10 mostly as one-off identification questions. Their real value is as a
PREREQUISITE GATE on every other branch above:
  - charge_type determines whether `system_type` for Gamma_max/deltaG_mic is even
    askable, and whether Nagarajan's ionic-only geometry branch applies at all.
  - structural_family (gemini/dimeric) should, in principle, change which mixture
    theories are even physically sensible to try in select_mixture_model (none of
    the 9 theories currently distinguish "is this component itself a gemini" --
    a real, currently-unimplemented interaction between the two NEW tools built
    this session, worth flagging as a possible re-module target, see below).
  - glycolipid/sugar-based structural family currently only gates the HLB-Davies
    refusal (GZ-04) -- it does NOT yet gate anything in the mass-action/Nagarajan/
    solubilization branches, even though a real researcher would ask "does the
    Nagarajan ionic-no-added-salt method even make physical sense for a
    zwitterionic-at-low-pH glycolipid?" -- a real, currently-open modeling question,
    not a tool bug, but worth a question that tests whether an LLM raises it.

## Confirmed gaps in the combination graph itself (not just missing MCP exposure)

1. **No tool asks "does this mixture-model choice interact with either component's
   own structural family?"** -- select_mixture_model and classify_surfactant_
   structural_family are both NEW this session and have never been chained together.
   A gemini-monomeric mixture (like a real system in this project's own literature
   validation, e.g. G6 + TX-114) is a real, different physical case from a
   monomeric-monomeric mixture, and none of the 9 mixture theories currently
   encode that distinction explicitly.
2. **Nagarajan's geometry branch and the plain Tanford-fixed-Amin branch have never
   been asked about on the SAME real dataset in the same question** -- GZ-10 only
   used the fixed-Amin path. A real "compute both, do they agree, which is more
   defensible for THIS real system" question is fully buildable right now with
   real real data (e.g. the AOT/SDS tensiometry curves already in this project),
   no new tool work needed.
3. **EOMMM's genuinely predictive x1-at-a-new-composition branch has never been
   asked** -- every mixture question so far only asked about the SAME compositions
   that were measured. A real "extrapolate to an unmeasured composition" question
   is a materially different, harder skill than fitting the measured points.
4. **HLD's cationic-quat-only HLB branch has never been paired with a real
   non-quat cationic surfactant to test whether the applicability limit is
   recognized** -- a clean, buildable refusal trap using real data already in
   this project (e.g. any of the real HLD-NAC systems already validated).
