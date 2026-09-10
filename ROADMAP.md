# SurfactantKit Roadmap

Resumable task list, written so a fresh session (no memory of this conversation) can pick
up exactly where this one left off. Read this file FIRST before starting new work.

Last updated: 2026-09-10.

---

## Reporting pipeline pilot -- DONE 2026-09-10

**New companion package built and working end to end**: `H:\CodeProjects\surfactantkit-reports`
(separate repo, git-initialized, not yet committed -- `pip install -e .` done, all deps
resolved: `surfactantkit`, `openpyxl`, `matplotlib`, `pygrace` 1.7). Followed the architecture
decided in the Pre-Execution Study below exactly (new companion package, pygrace for `.agr`,
matplotlib for PNG as two independent rendering paths).

Pilot technique (`cmc_from_surface_tension_curve`) taken to 100%:
- `excel_io.py` -- reads the 2-column template, validates by COLUMN NAME (not position, so
  reordering is fine but renaming is caught), raises a specific error naming the exact
  problem for every identified failure mode (missing column, non-numeric cell, empty sheet,
  blank rows skipped cleanly).
- `templates.py` -- generates the blank Excel template (units spelled out in the header
  itself: "Concentration (mM)", "Surface Tension (mN/m)"); a real pre-built copy committed at
  `templates/cmc_from_surface_tension_template.xlsx`.
- `report_data.py` + `plotting.py` -- the real fix for the identified "PNG and .agr silently
  disagree" failure mode: both `render_png` and `render_agr` consume ONLY a single shared
  `ReportData` object built once from the fit result, never recomputing anything
  independently.
- `agr_verify.py` -- since pygrace has no load/parse API (confirmed by inspection -- it is
  write-only) and this project cannot visually render Grace/xmgrace/QtGrace, built an
  independent (non-pygrace) regex-based re-parser and a `verify_agr_matches_report_data()`
  check that is a REAL regression gate in the test suite, not just documentation -- the
  disclosed partial substitute for a visual check, per the Pre-Execution Study.
- `pipeline.py` -- orchestrates the whole thing; `run_cmc_surface_tension_report(excel_path,
  output_dir)` is the one public entry point for this technique.

Validated against the exact same real pilot dataset already used in SurfactantKit's own
`curve_analysis.py` `__main__` block (AOT in water, Shah/Das/Bhattarai 2025, Heliyon,
PMC11835642) -- pipeline recovers CMC = 2.512 mM against the paper's own reported 2.51 mM.
PNG visually inspected directly (log-x axis, premicellar fit line exactly overlaid on the
declining data, CMC vertical marker exactly at the break point, postmicellar points visible,
annotation box with the fit's own R^2/point-count) -- matches expectations exactly, not just
"a plot got produced." Real bug found and fixed during the pilot:
`graph.add_drawing_object("DrawLine", ...)` fails (`add_drawing_object` needs the actual
class object, not a string) -- fixed to `add_drawing_object(DrawLine, ...)` after importing
`from pygrace.drawing_objects import DrawLine`. 13 tests in
`tests/test_pipeline_cmc_surface_tension.py`, all passing (Excel read/write/validation, full
pipeline round trip, the PNG/agr-agreement regression gate, axis-label-units check,
missing-file handling).

**Second technique (`cmc_from_conductivity_curve`) also DONE 2026-09-10.** Rolled out the
same pattern: `read_cmc_conductivity_template`/`write_cmc_conductivity_template` (headers
"Concentration (mM)" / "Conductivity" -- unit left generic since the method is scale-
invariant in the y-axis unit), `build_cmc_conductivity_report_data`, `run_cmc_conductivity_report`.
Real design point: `cmc_from_conductivity_curve` doesn't return the two fitted lines'
intercepts (only slopes + point counts), so the two lines are reconstructed for the PLOT ONLY
via an independent local OLS re-fit, split at the EXACT index the result's own
`n_points_below_cmc`/`n_points_above_cmc` specify (not a re-search) -- and the re-fit slopes
are asserted in tests to match the result's own authoritative `slope_below_cmc`/`slope_above_cmc`
almost exactly, a real consistency check. Validated on the same synthetic construction already
used in SurfactantKit's own `test_curve_analysis.py` (no external raw-data table exists for
this technique, disclosed there already) -- recovers the true CMC and beta exactly, PNG
visually confirmed (two clean linear segments meeting at the break point, linear not log
x-axis, correctly different from the tensiometry technique). 7 new tests in
`tests/test_pipeline_cmc_conductivity.py`. Package total: 20/20 tests passing across both
techniques.

**Third and fourth techniques (SSFQ + SLS) also DONE 2026-09-10 -- all 4 originally-planned
techniques now complete.**

`aggregation_number_from_quenching_curve` (SSFQ) and `aggregation_number_from_sls_debye_plot`
(SLS) both needed a template layout the first two techniques didn't: several SCALAR values
(I0/Ct/CMC for SSFQ; dn/dc/wavelength/monomer-molar-mass/refractive-index for SLS) alongside
the (x, y) raw-data series. Solved with a new, reusable template shape -- a labeled metadata
block (label in column A, value in column B, matched by exact label TEXT not fixed row number,
so reordering the metadata rows doesn't break the reader) directly above the data table in the
SAME sheet -- implemented once as `_read_metadata_and_data_template`/`_write_metadata_and_data_template`
in `excel_io.py`/`templates.py`, both specific readers/writers now thin wrappers around it.

Both techniques' plots show the LINEARIZED data the underlying fit actually uses (ln(I0/I) vs
[Q] for SSFQ; Kc/DeltaR vs concentration for SLS Debye plot), not the raw instrument readings
directly -- this is the standard, correct way both are plotted (the whole point of
linearizing is a straight-line fit), computed via the exact same transform formula the library
function itself uses internally (re-derived in the reporting package, not re-fit
independently) so the plotted points are guaranteed consistent with the reported result. Both
validated via the same mathematically-guaranteed synthetic construction already used in
SurfactantKit's own `test_curve_analysis.py` (no external raw-data table exists for either,
already disclosed there) -- exact recovery of the true aggregation number in both cases, PNGs
visually confirmed (clean zero-intercept line for SSFQ, correct Debye-plot line for SLS). 16
new tests across `tests/test_pipeline_ssfq.py` and `tests/test_pipeline_sls.py`.

**Package total: 36/36 tests passing across all 4 techniques** (`cmc_from_surface_tension_curve`,
`cmc_from_conductivity_curve`, `aggregation_number_from_quenching_curve`,
`aggregation_number_from_sls_debye_plot`). Every technique has: a reader with specific,
non-guessing error messages; a blank template committed under `templates/`; a `ReportData`
builder that reconstructs (never re-fits) the fitted line(s) from the library result's own
returned fields; PNG + `.agr` output from that single shared object; and a real
`verify_agr_matches_report_data` regression gate against the "PNG and .agr silently disagree"
failure mode in every test file. `surfactantkit-reports` is git-initialized at
`H:\CodeProjects\surfactantkit-reports` but not yet committed/pushed -- ask before doing so.

**Not started / explicitly out of scope for this package**: any technique not already
implemented in `curve_analysis.py` (i.e. everything in SurfactantKit that still takes a
pre-computed value rather than a raw dataset -- see the "Also still open" list below and
Phase 8 above) has nothing to build a report template FOR yet; a report front-end is only
possible once the underlying raw-data extraction function exists.

---

## Immediate next task, ready to start Thursday: raw-data reporting pipeline

**User's actual ask** (verbatim intent): a researcher without SANS/TRFQ/ultracentrifuge-tier
lab access should be able to run a SIMPLE experiment, fill in an Excel template with the raw
readings, hand it to our tool, and get back: a fitted result, a PNG plot with correct axis
labels for that technique, and a `.agr` file (Grace/xmgrace/QtGrace native format) the
researcher can open and further edit by hand. "Automatic iterations wherever needed" --
i.e. wherever the underlying method already does an iterative/regression fit (most of them
do), the pipeline just needs to call that, not add new iteration logic.

### Architecture decision (made, confirmed by user 2026-09-08)

**New companion package**, not inside SurfactantKit itself. SurfactantKit's own
`pyproject.toml` currently has `dependencies = []` (pure stdlib, zero deps) by design --
keep it that way. The new package (name TBD, e.g. `surfactantkit-reports`) depends on
SurfactantKit itself plus `openpyxl`, `matplotlib`, and `pygrace`.

### Key finding from this session's Pre-Execution Study (do NOT re-research this)

**`pygrace`** (PyPI, `uqfoundation/pygrace`, latest release 1.7, Jan 2026, actively
maintained, `pip install pygrace`, requires Python >=3.9) builds `.agr` files
programmatically via a clean object API -- **do not hand-write Grace's raw file grammar**,
this eliminates that entire risk class. Real verified usage pattern (fetched directly from
`github.com/uqfoundation/pygrace/blob/master/examples/01_singleplot.py`):

```python
from pygrace.project import Project
grace = Project()
graph = grace.add_graph()
graph.xaxis.label.text = 'X label'
graph.yaxis.label.text = 'Y label'
dataset = graph.add_dataset(data, legend='Set 1')  # data = list of (x, y) tuples
graph.set_world_to_limits()
grace.saveall('output.agr')
```

A different, related tool (`goerz/xmgrace_parser`, GitHub) was checked and explicitly
rejected for this purpose -- its own README states it's for editing an EXISTING template
agr file, not creating one from scratch, and it recommends pygrace for that. Don't
reconsider it unless the pygrace approach hits a real wall.

**PNG rendering**: use `matplotlib` directly (already installed in this environment, full
control), independently of the `.agr` generation -- NOT by shelling out to
`gracebat`/xmgrace in batch/hardcopy mode. Two independent rendering paths deliberately, so
PNG generation doesn't require Grace to be installed on the machine running the pipeline
(only needed on whichever machine the researcher later opens the `.agr` on).

**`local_repo_toolkit_inventory.md`** (referenced in the user's global CLAUDE.md as a thing
to check before building anything new) **does not exist on this machine** -- confirmed by
direct search this session. `D:\ClaudeAssets\toolkit-repos\` (the real folder per the
CLAUDE.md moved-folders ledger) has only 6 repos, none relevant here. This is a stale
reference in the user's own instructions -- flag it again if it comes up, don't silently
work around it or assume the file exists.

### Proposed pilot (do this first, before any other tool)

Build the full pipeline (Excel template -> reader -> existing fitting function -> PNG ->
`.agr`) end-to-end for **`cmc_from_surface_tension_curve`** (tensiometric CMC) specifically
-- it's already implemented and tested in `curve_analysis.py`, so the pilot only has to
prove the NEW plumbing works, not debug new science at the same time. Once solid (including
an actual visual check of the PNG, and ideally getting real Grace/QtGrace installed
somewhere to visually confirm the `.agr` renders correctly, not just that it parses), roll
the same pattern out to the other tools one at a time, per the user's explicit "one at a
time, 100% done before moving on" instruction.

### Failure modes already identified (from the Pre-Execution Study, still valid)

- Excel template drift (columns renamed/reordered over time) -- validate on read, reject
  with a clear error naming the missing/misnamed column, never guess.
- PNG and `.agr` silently disagreeing (different axis ranges/fit shown) -- since they're two
  independent rendering paths, they must be built from the exact same intermediate
  data/fit-result object, not recomputed separately.
- Can't visually verify `.agr` rendering without Grace installed locally -- gate on the file
  being re-parseable by a second script and on its embedded data values exactly matching the
  PNG's, as a partial substitute for an actual visual check.
- Units ambiguity in the Excel template (mM vs M, mN/m vs dyn/cm) -- label units explicitly
  in the template itself, don't rely on the researcher inferring them.

---

## Full audit of "raw dataset -> real answer" gaps across the toolkit

Done this session (2026-09-08), verified against actual code (not recalled from memory).
Only `curve_analysis.py` currently goes from a raw multi-point dataset to a real answer
(`cmc_from_surface_tension_curve`, `aggregation_number_from_quenching_curve`,
`aggregation_number_from_sls_debye_plot`, all three built and tested this session). Every
other module's functions take pre-computed/pre-known single values as input. Real,
confirmed gaps, in priority order (feed into the "one at a time" rollout of the reporting
pipeline above -- each of these also needs a fitting function built/verified BEFORE it can
get an Excel template, if it doesn't already exist):

1. **Conductometric CMC extraction** -- DONE 2026-09-10. Built `cmc_from_conductivity_curve()`
   + `CmcFromConductivityResult` in `curve_analysis.py`: two-segment linear regression on
   LINEAR concentration (not log, unlike the surface-tension version), break point at the
   fitted lines' intersection = CMC, plus `counterion_binding_degree` returned directly using
   the exact formula `thermodynamics.counterion_binding_degree()` implements. No clean
   fully-reproducible raw (concentration, conductivity) literature table was found during this
   pass (PMC8007100 gives only a summary CMC of 8.40+/-1.14 mM, no raw table) -- same
   disclosed "no clean external numeric example" situation as `aggregation_number_spherical`
   and the SSFQ/SLS tools; validated instead via a mathematically-guaranteed round-trip plus a
   real literature plausibility range (SDS CMC ~8.1-8.4 mM, beta ~0.6-0.8) and a cross-check
   against `thermodynamics.counterion_binding_degree()`. 5 new tests in
   `tests/test_curve_analysis.py` (round-trip, literature-range plausibility, cross-check
   against thermodynamics module, bad-input rejection, non-ionic-flat-data rejection) --
   14/14 passing. Wired into `__init__.py` (import + `__all__`) and `mcp_server.py` (new
   `cmc_from_conductivity` tool), registered + end-to-end tested in `test_mcp_server.py`.
   Full suite: 188/188 passing.
2. **Multi-point Rubingh beta regression** -- DONE 2026-09-10. Verified the real standard
   practice via WebSearch before implementing (not guessed): published papers report a single
   system beta as the arithmetic MEAN of the pointwise betas computed independently at each
   composition in the series -- there is no single closed-form cmc_mix(alpha1) curve to
   least-squares-fit against directly, since x1 itself must be solved pointwise by
   construction of RST. Built `rubingh_beta_regression()` + `RubinghBetaRegressionResult` in
   `mixed_micelle.py`: loops `solve_rubingh_x`/`rubingh_beta` over the series, returns
   `beta_mean`, `beta_std` (a real diagnostic -- large spread signals RST is only an
   approximation for that system), and the per-point betas/x1 values; points with no valid
   root are skipped (per `solve_rubingh_x`'s existing None contract), not errored, unless
   fewer than 2 points solve. Validated via a mathematically-guaranteed round trip: derived
   the closed-form RST mass-balance construction `cmc_mix = x1*f1*cmc1 + (1-x1)*f2*cmc2`,
   `alpha1 = x1*f1*cmc1/cmc_mix` to build a composition series that all share one true beta by
   construction, then confirmed the regression recovers it exactly. Found and worked around a
   real solver limitation in the process: `solve_rubingh_x`'s grid+bisection scan returns the
   FIRST sign change, which is not always the intended root when beta is large and positive
   (antagonistic) enough to make the residual non-monotonic -- verified numerically (beta=2.2
   showed root mismatches at 3 of 4 test compositions, beta=1.8 was clean), so the round-trip
   test uses beta magnitudes confirmed single-rooted; this pre-existing pointwise-solver
   behavior is unchanged, just documented here rather than silently worked around. 6 new tests
   in `tests/test_mixed_micelle.py` (round-trip recovery, antagonistic-sign case, skip-not-error
   for an unsolvable point via monkeypatch, bad-input rejection) -- all passing. Wired into
   `__init__.py` and `mcp_server.py` (new `rubingh_beta_regression` tool), registered +
   end-to-end tested in `test_mcp_server.py`. Full suite: 193/193 passing.
3. **HLD salinity-scan fitting** -- DONE 2026-09-10. Built `fit_k_and_cc_from_salinity_scan()`
   + `HldSalinityScanFitResult` in `hld.py`: fits ln(S*) vs EACN by OLS (k=slope,
   Cc=alpha*deltaT-intercept, per the source paper's own Eq. 5 inverted), returns k, Cc,
   r_squared (a real diagnostic -- HLD-NAC assumes k is constant across EACN for a surfactant
   class; poor r_squared signals that assumption breaks for the system under test) and the
   points used. No external multi-EACN raw dataset was sourced this session (same disclosed
   "no clean external numeric example" pattern as several other tools built this session);
   validated instead via mathematically-guaranteed round trips (construct S* from a chosen
   true k/Cc via the already-built `optimal_salinity_ionic`, fit them back exactly, including
   a nonzero-temperature case) plus an independent manual-OLS cross-check. 8 new tests in
   `tests/test_hld.py` -- 16/16 passing. Wired into `__init__.py` and `mcp_server.py` (new
   `hld_fit_k_and_cc_from_salinity_scan` tool), registered + end-to-end tested in
   `test_mcp_server.py`. Full suite: 198/198 passing.
4. **Multi-point van't Hoff / Gibbs-Helmholtz fit** -- DONE 2026-09-10. Verified Kantonen,
   Henriksen & Gilson 2018's exact Eq. 4 via WebFetch before implementing (not guessed):
   deltaG(Tj) = deltaH(Tr) - Tj*deltaS(Tr) + deltaCp(Tr)*[(Tj-Tr) - Tj*ln(Tj/Tr)]. Real finding
   made while implementing: although the source paper describes fitting this by "nonlinear
   optimization," the equation is actually LINEAR in the three unknowns (deltaH, deltaS,
   deltaCp) for a fixed reference temperature Tr -- the Tj-dependence lives entirely in the
   three regressor coefficients. Built `vant_hoff_multi_point_fit()` +
   `VantHoffMultiPointFitResult` in `thermodynamics.py` solving the exact closed-form OLS
   normal-equations system (new local `_solve_3x3` Gaussian-elimination helper, no numpy/scipy
   dependency added) instead of iterative nonlinear optimization -- mathematically equivalent
   for fixed Tr, confirmed via exact round-trip recovery. `reference_temperature_K` defaults to
   the mean of the input temperatures (well-conditioned; overridable), with deltaH/deltaS
   correctly reported as reference-dependent while deltaCp and the fitted deltaG(T) curve are
   reference-independent (verified in a dedicated test with two different reference choices).
   9 new tests in `tests/test_thermodynamics.py` (round-trip recovery, default-reference case,
   ionic counterion_factor case, bad-input rejection including a genuinely degenerate/singular
   all-same-temperature case) -- 17/17 passing. Wired into `__init__.py` and `mcp_server.py`
   (new `vant_hoff_multi_point_fit` tool, taking a raw CMC-molarity series + temperature
   series), registered + end-to-end tested in `test_mcp_server.py`. Full suite: 203/203
   passing.
5. **Rodenas' derivative-from-series helper** -- DONE 2026-09-10. Checked via WebSearch/WebFetch
   before implementing whether the Rodenas source papers (or papers using the method, e.g.
   PMC4087020) specify their exact numerical-differentiation procedure -- they don't (they
   describe reading "the slope from the curve" without a specified numerical method), so this
   uses the real, standard, general numerical-differentiation technique for an unevenly-spaced
   discrete series rather than a guessed surfactant-specific one: at each composition, fit the
   unique quadratic (Lagrange basis) through it and its 2 nearest neighbors along alpha1, then
   differentiate that quadratic analytically at that exact alpha1. Built `rodenas_x1_series()`
   + `RodenasSeriesResult` + local `_quadratic_lagrange_derivative_at()` helper in
   `mixed_micelle.py`, chaining straight through the existing `rodenas_x1()` to return x1 at
   every composition, not just the slope. Verified the method is a real generalization (not
   invented) two ways: (1) proved the general Lagrange-derivative formula reduces EXACTLY to
   the textbook central-difference formula (y2-y0)/(2h) at an evenly-spaced interior point,
   and to the textbook 3-point forward-difference formula (-3y0+4y1-y2)/(2h) at an endpoint --
   both confirmed in tests, not just asserted; (2) confirmed exact (not approximate) recovery
   when the underlying ln(cmc_mix) truly is quadratic in alpha1 (a local quadratic fit is then
   not an approximation at all). 8 new tests in `tests/test_mixed_micelle.py` (exact-quadratic
   recovery, central/forward-difference reduction, chaining into `rodenas_x1`, input-order
   sorting, bad-input rejection including duplicate alpha1) -- 45/45 passing. Wired into
   `__init__.py` and `mcp_server.py` (new `rodenas_x1_from_series` tool), registered +
   end-to-end tested in `test_mcp_server.py`. Full suite: 210/210 passing.
6. **`szyszkowski_fit_K`** -- DONE 2026-09-10. Was a real, confirmed loose thread: direct grep
   found `adsorption.py`'s own `szyszkowski_surface_tension` docstring referencing this
   function by name ("K... fit from real data, e.g. via szyszkowski_fit_K") even though it did
   not exist anywhere in the codebase. Built it as a real nonlinear least-squares fit -- a
   DIFFERENT curve regime/purpose than `cmc_from_surface_tension_curve` (which locates the CMC
   break point): this one fits the pre-CMC region's actual monomer-adsorption curvature to
   extract K, with `gamma_max_mol_per_m2` taken as an already-known input from the same curve's
   pre-CMC log-slope (`gibbs_gamma_max`), matching the real, standard two-step surface-tension-
   isotherm workflow. Solved via a new local `_golden_section_minimize()` helper over ln(K)
   (the least-squares objective is unimodal in K since the model's predicted surface tension is
   monotonic in K) rather than adding a general nonlinear-solver dependency. Returns K,
   r_squared, and the fitted curve for residual inspection. Validated via mathematically-
   guaranteed round trips (construct a curve from a chosen true K via
   `szyszkowski_surface_tension` itself, fit it back exactly, for both a nonionic and an
   ionic-no-added-salt system_type) plus a real diagnostic check (a curve that does not follow
   the Szyszkowski shape at all fits with r_squared << 1, not silently ~1). 5 new tests in the
   new `tests/test_szyszkowski_fit_k.py` -- 5/5 passing. Wired into `__init__.py` and
   `mcp_server.py` (new `szyszkowski_fit_K` tool), registered + end-to-end tested in
   `test_mcp_server.py`. Full suite: 216/216 passing.

**All 6 items of the raw-dataset-gap audit are now complete** (2026-09-10), each reaching the
same "100% done, tested, documented, ready to push" bar as everything else built this session
(Frumkin, Motomura, EOMMM, Rodenas, Maeda, Ohshima/Swan-Furst/Qin regimes, the HLD module, the
quaternary-ammonium HLB group numbers, and the two aggregation-number tools): conductometric
CMC extraction, multi-point Rubingh beta regression, HLD salinity-scan fitting, the multi-point
van't Hoff/Gibbs-Helmholtz (deltaCp) fit, Rodenas' derivative-from-series helper, and
szyszkowski_fit_K. SurfactantKit's test suite grew from 188 (session start) to 216 tests, all
passing, across this push.

---

## Phase 8: full alternative-methods sweep (2026-09-10, ongoing)

User's explicit standing ask: "for each tool present in our toolkit, might have some
alternative models/method... i need all alternate method to be built in our toolkit. so that
research can chose which one suitable for them." Working through every item flagged in
`benchmark/METHOD_ALTERNATIVES_LITERATURE_REVIEW.md` as a real, buildable (formula-based, not
needing external training data or being purely experimental) alternative method, one at a time
to the same 100%-done bar, real literature verification before implementation each time.

1. **EOMMM global fit** -- DONE (see Phase 7 above, done just before this phase started).
2. **Nagarajan (2002) ionic-headgroup equilibrium-area model** -- DONE 2026-09-10. Built
   `nagarajan_debye_huckel_kappa_inverse()` and `nagarajan_equilibrium_area_ionic()` in
   `cpp.py`: a real, sourced alternative to feeding a constant/assumed head_area_A2 into
   `critical_packing_parameter` -- shows a_e (and hence CPP) genuinely varies with tail length
   for IONIC surfactants via the Debye screening length (which depends on cmc), contrary to
   the common assumption that only the headgroup sets a_e. Verified via WebFetch (vision-read
   PDF, ScienceDirect-style extraction blocks worked around the same way as earlier this
   session) directly against the source paper's own Table 2 (real numeric worked example,
   sodium alkyl sulfates n_C=8..16) -- kappa^-1 and a_e both reproduce the table to within
   ~0.2%, and the resulting CPP is confirmed to strictly decrease across the homologous series
   matching the table's own v0/(a_e*l0) column. `headgroup_prefactor_A` (the d/sigma
   electrostatic prefactor) is required, not defaulted -- the paper's own 82.0 Angstrom value
   is disclosed as illustrative for sodium-alkyl-sulfate-like headgroups specifically, not a
   universal constant. 6 new tests in `tests/test_hlb_cpp.py` -- 39/39 passing in that file.
   Wired into `__init__.py` and `mcp_server.py` (new `nagarajan_equilibrium_area_ionic` tool),
   registered + end-to-end tested. Full suite: 234/234 passing.

3. **Perrin friction factors** -- DONE 2026-09-10. Built `perrin_friction_factor()` and
   `hydrodynamic_radius_perrin_corrected()` in `dynamics.py`: a real, sourced alternative to
   plain Stokes-Einstein (`hydrodynamic_radius_stokes_einstein`), which assumes a sphere and
   so returns an inflated apparent radius for rodlike/wormlike or disklike aggregates. Formula
   verified via WebFetch (Perrin 1934, as commonly tabulated e.g. Wikipedia's "Perrin friction
   factors"): F_P = 2*p^(2/3)/S, S=2*atanh(xi)/xi (prolate) or 2*atan(xi)/xi (oblate),
   xi=sqrt(|p^2-1|)/p. No external worked numeric table was found/fetchable this pass (two
   PDF attempts failed to yield a clean table), so validated via the one point provable
   directly from the formula itself (sphere limit p=1 gives F_P=1 exactly) plus real physical-
   consistency checks (F_P>1 away from p=1 in both directions, since a sphere is the minimum-
   drag shape for a given volume; strictly monotonic with elongation) -- disclosed as a
   narrower validation than the Nagarajan/EOMMM items above, not overstated. 10 new tests in
   `tests/test_electrostatics_dynamics.py` -- 33/33 passing in that file. Wired into
   `__init__.py` and `mcp_server.py` (new `hydrodynamic_radius_perrin_corrected` tool),
   registered + end-to-end tested. Full suite: 244/244 passing.
4. **Guo/Rong/Ying 2006 nonionic HLB refinement** -- DONE 2026-09-10. Built
   `guo_effective_eo_chain_length()`, `guo_effective_alkyl_chain_length()`,
   `guo_effective_po_chain_length()`, and `hlb_davies_guo_ecl()` in `hlb.py`: a real,
   sourced alternative to `hlb_davies` specifically for NONIONIC polyethoxylated
   surfactants, using "effective chain length" (sub-linear in actual EO count) instead of
   actual chain length, per the source paper's own real validation (224 surfactants, average
   absolute error < 1.5 vs. Davies). Primary source (J. Colloid Interface Sci. 298 (2006)
   441-450) is paywalled; exact coefficients (NEOeff=13.45*ln(NEO)-0.16*NEO+1.26,
   NCH2eff=0.965*NCH2-0.178, NPOeff=2.057*NPO+9.06) verified via a citing patent (US
   11,344,493 B2) that quotes the formulas directly -- NOT reconstructed or guessed.
   Honestly disclosed limitation: whether Guo/Rong/Ying also re-fit new per-unit
   hydrophilic/lipophilic weights (vs. reusing Davies' 1.3/0.475) was not confirmed from the
   primary source -- `hlb_davies_guo_ecl` assumes Davies' original weights apply to the
   effective chain lengths, documented as a reasonable but unverified reading, not an exact
   reproduction of the paper's own numbers. PO chain support is exposed as a standalone
   utility only (no verified Davies-scale PO group number exists in this project to chain it
   into a full HLB). 9 new tests in `tests/test_hlb_cpp.py` -- 50/50 passing in that file.
   Wired into `__init__.py` and `mcp_server.py` (new `hlb_davies_guo_ecl` tool), registered +
   end-to-end tested. Full suite: 256/256 passing.
5. **Rusanov/classical mass-action model** -- DONE 2026-09-10. Built
   `mass_action_free_energy_of_micellization()` in `thermodynamics.py`: a real, sourced
   alternative to `gibbs_free_energy_micellization`'s pseudo-phase-separation convention,
   `deltaG_mic = RT*[(1+1/n)*ln(Xcmc) - (1/n)*ln(n)]`. The primary source flagged in the
   literature review (Rusanov, Langmuir 30 (2014) 14443-14451) is paywalled; the same
   well-known classical single-equilibrium mass-action result was independently verified via
   an open-access paper instead (Boneva-Aroca et al., RSC Advances 13 (2023) 9387, their eq.
   21, quoted verbatim via WebFetch) -- not guessed or reconstructed from memory. Validated by
   confirming the real, paper-stated property that this formula converges EXACTLY to
   `gibbs_free_energy_micellization` as n->infinity (verified numerically: the gap shrinks
   monotonically and vanishes by n=100000), and that the correction is real and non-negligible
   at small n (>0.5 kJ/mol at n=15) -- ties the new model concretely back to the existing one
   rather than floating unconnected. NONIONIC-ONLY, explicitly disclosed: extending to ionic
   surfactants needs a different mass-action treatment (explicit counterion-binding
   equilibrium) not verified this pass. 5 new tests in `tests/test_thermodynamics.py` --
   21/21 passing in that file. Wired into `__init__.py` and `mcp_server.py` (new
   `mass_action_free_energy_of_micellization` tool), registered + end-to-end tested. Full
   suite: 261/261 passing.
6. **O'Brien-White relaxation-effect correction** -- DONE 2026-09-10. Built
   `electrophoretic_mobility_relaxation_corrected()` (forward) and
   `zeta_potential_relaxation_corrected()` (numerical inverse, grid+bisection like
   `solve_rubingh_x`) in `electrostatics.py`: a real, sourced alternative to the plain Henry
   equation for HIGH zeta (>50 mV, common for ionic surfactant micelles) -- Henry is linear in
   zeta and misses the relaxation effect (mobility grows sub-linearly with zeta, can pass
   through a maximum beyond ~100 mV) entirely. The primary source (O'Brien & White 1978) and
   the closed-form OHW83 semi-empirical approximation both proved hard to pin down exactly
   (several paywalled/inaccessible attempts); found the real, usable, citable formula instead
   in an IUPAC Technical Report (Delgado et al., J. Colloid Interface Sci. 309 (2007) 194-224,
   their Eq. 24 -- O'Brien's own simplification of the Dukhin & Semenikhin relaxation
   equation), fetched and read directly via the PDF's real extractable text layer (not vision-
   read this time -- a genuine text PDF, confirmed via pypdf). Validated via a real, strong
   cross-check tying it back to already-validated existing code: at low zeta and large
   kappa_a, this formula reduces to within 0.1% of the plain Smoluchowski limit of
   `zeta_potential_henry`/`henry_function` (f(kappa*a)=1.5) -- confirmed algebraically and
   numerically. Also confirmed the real physical signature (sub-linear mobility growth at high
   zeta) and a clean round-trip (forward -> inverse solver) across 10-90 mV. Restricted to
   kappa_a > ~20 (the source's own stated valid range) and zeta in [0, 150] mV by default for
   the inverse solve (below the reported non-monotonic-mobility maximum) -- both limits
   enforced with explicit errors, not silently guessed past. 8 new tests in
   `tests/test_electrostatics_dynamics.py` -- 39/39 passing in that file. Wired into
   `__init__.py` and `mcp_server.py` (new `predict_mobility_relaxation_corrected` and
   `zeta_potential_relaxation_corrected` tools), registered + end-to-end tested. Full suite:
   268/268 passing.
7. **OWRK / van Oss-Chaudhury-Good surface energy decomposition** -- DONE 2026-09-10. Built
   `owens_wendt_solid_surface_energy()` (2-component, OLS across 2+ test liquids) and
   `van_oss_chaudhury_good_solid_surface_energy()` (3-component LW/acid/base, closed-form 3x3
   linear solve across exactly 3 test liquids) in `wetting.py`: real, sourced alternatives to
   `work_of_adhesion`'s single scalar number, needed when decomposing solid surface energy
   into dispersive/polar (OWRK: Owens & Wendt 1969; Rabel 1971; Kaelble 1970) or Lifshitz-van
   der Waals/Lewis acid-base (vOCG: van Oss, Chaudhury & Good, Chem. Rev. 88 (1988) 927-941)
   components matters -- not a more-correct answer to the same question `work_of_adhesion`
   answers. Both formulas verified via WebSearch before implementing. Real design decision:
   both mixing rules combine with this project's OWN existing `work_of_adhesion` (Young-Dupre)
   rather than reimplementing it, and both reduce to LINEAR systems in sqrt(gamma_S component)
   terms (OWRK: 2-parameter OLS; vOCG: exact 3x3 linear solve, correctly using the vOCG
   cross-terms where the solid's acid component pairs with each liquid's BASE component, not
   acid-acid). No external worked numeric table was sourced this pass; validated via
   mathematically-guaranteed round trips (construct contact angles from chosen true solid
   surface-energy components via the forward mixing rule, recover them back exactly) plus a
   real physical sanity check (a purely-dispersive solid recovers exactly zero polar
   component). 9 new tests in `tests/test_szyszkowski_wetting_solubilization.py` -- 21/21
   passing in that file. Wired into `__init__.py` and `mcp_server.py` (new
   `owens_wendt_solid_surface_energy` and `van_oss_chaudhury_good_solid_surface_energy`
   tools), registered + end-to-end tested. Full suite: 275/275 passing.
8. **K_M partition coefficient** -- DONE 2026-09-10. This module's own docstring had
   PREVIOUSLY made an explicit, reasoned decision NOT to implement Km at all (genuine
   literature ambiguity: mole-fraction-basis vs. molar-concentration-ratio-basis conventions,
   no single dominant one). Re-checked this pass -- the ambiguity is real and still
   unresolved, but rather than continue a blanket omission, built
   `micelle_water_partition_coefficient()` in `solubilization.py` implementing the MOLE-
   FRACTION-BASIS convention specifically (Km = X_solubilizate,micelle/X_solubilizate,water),
   with the convention ambiguity prominently disclosed in the function's own docstring (the
   molar-concentration-ratio convention is explicitly NOT implemented, still a real, disclosed
   gap) -- letting a researcher use it with eyes open rather than not at all, matching the
   "let a researcher choose" goal. Reuses the exact same 4 raw-data inputs as the existing
   `molar_solubilization_ratio` for consistency, and the same dilute-water-molarity constant
   (55.5 mol/L) already used in `thermodynamics.cmc_to_mole_fraction`. 6 new tests in
   `tests/test_szyszkowski_wetting_solubilization.py` (manual-formula match, physical sanity
   checks for strongly-vs-weakly-solubilized species, bad-input rejection) -- 27/27 passing in
   that file. Wired into `__init__.py` and `mcp_server.py` (new
   `micelle_water_partition_coefficient` tool), registered + end-to-end tested. Full suite:
   282/282 passing.
9. **Gouy-Chapman/Stern rigorous electrostatic treatment** -- PARTIAL, 2026-09-10, with an
   important honest finding. Found the real self-consistent n_factor correction (Mchedlov-
   Petrossyan, Electron. Proc. Mater. 50(2) (2014) 71-80, Eq. 8) requires a term
   (d ln Gamma_max/d ln a2) this session's "chain existing tools" approach could not close
   cleanly: a first attempt built `gibbs_gamma_max_gouy_chapman()` in `adsorption.py`,
   chaining `counterion_binding_degree` + `ionic_strength` + a new Grahame-equation function,
   via a closed-form analytic dPsi/d ln a2 that held the surface charge density fixed at the
   evaluation point. Numerically debugging it (checking the beta=0 and beta=1 physical limits
   against expectation) surfaced a real units bug (Faraday's constant vs. bare elementary
   charge when converting a MOLAR Gamma_max to a charge density -- fixed), and THEN a deeper
   physical gap: holding q fixed while differentiating drops the very term (d ln Gamma/d ln a2)
   the full theory needs, and produced numerically unstable / physically implausible results
   at the beta=0 extreme even after adding damping and clamping. Rather than ship a
   plausible-looking but uncertain correction, **removed `gibbs_gamma_max_gouy_chapman`
   entirely** -- a real, deliberate decision not to guess past a genuine remaining gap, not an
   oversight. KEPT and shipped the piece that IS fully correct and independently verified:
   `grahame_equation_surface_potential()` in `electrostatics.py` -- the standard, textbook
   closed-form Gouy-Chapman relation (Grahame 1947), cross-checked against a real fetched
   secondary source's same functional form AND against a real, derivable low-potential limit
   (reduces exactly to eps_r*eps0/debye_length, tying it back to this project's own already-
   validated `debye_length`). 9 new tests in `tests/test_electrostatics_dynamics.py` (round
   trip, zero-charge limit, sign check, the debye_length cross-check, bad-input rejection) --
   44/44 passing in that file. Wired into `__init__.py` and `mcp_server.py` (new
   `grahame_equation_surface_potential` tool). A real, standalone, useful Gouy-Chapman tool
   for computing adsorbed-monolayer surface potential from charge density -- but the
   originally-scoped "alternative n_factor for gibbs_gamma_max" remains genuinely unresolved
   and is NOT implemented; do not assume it exists.
10. FLM reorientation/multistate adsorption model (adsorption.py) -- alternative to plain
    Langmuir/Frumkin for surfactants that can adopt multiple interfacial orientations. Given
    item 9's finding (a structurally similar multi-state/self-consistent model), deferred
    without attempting -- likely has the same class of closure-term difficulty.
11. **Additional aggregation-number raw-data routes** -- DONE (DLS + ultracentrifugation),
    2026-09-10, with one deliberately dropped. Built `aggregation_number_from_dls()` (Stokes-
    Einstein hydrodynamic radius -> volume -> molar mass via partial specific volume,
    directly reusing the same Stokes-Einstein formula as `dynamics.hydrodynamic_radius_stokes_einstein`,
    with a real, disclosed hydration-shell bias caveat) and
    `aggregation_number_from_svedberg_equation()` (the classic analytical-ultracentrifugation
    Svedberg equation, M=RTs/[D*(1-v_bar*rho)], independent physical principle from the
    light-scattering routes -- a genuine cross-check when both are available) in
    `curve_analysis.py`. Both verified via WebSearch against standard textbook forms before
    implementing (not guessed), and both give physically sensible aggregation numbers (30-90
    range) for realistic test inputs. **Deliberately did NOT build a viscosity-based route**:
    checked the real physics (Einstein's viscosity law, [eta]=2.5*v_bar for hard spheres) and
    found intrinsic viscosity alone is INSENSITIVE to aggregation number for a compact
    spherical micelle -- it depends only on shape/hydration (v_bar), not size, so a "viscosity
    -> Nagg" formula would either be circular or require smuggling in an assumption that isn't
    real physics. A genuine, disclosed finding, not an oversight (same "do not guess/force a
    method that isn't real" discipline as the K_M and PIT findings elsewhere in this project).
    8 new tests in `tests/test_curve_analysis.py` (manual-formula match for both, realistic-
    range sanity checks, non-physical-buoyancy rejection, bad-input rejection) -- 21/21
    passing in that file. Wired into `__init__.py` and `mcp_server.py` (new
    `aggregation_number_from_dls` and `aggregation_number_from_svedberg` tools), registered +
    end-to-end tested. Full suite: 297/297 passing.

**This closes the full-day alternative-methods sweep (Phase 8).** Of the original review's
11-item list: 9 fully built (EOMMM, Nagarajan, Perrin, Guo/Rong/Ying, mass-action model,
O'Brien-White/relaxation-corrected mobility, OWRK/vOCG, K_M, DLS+Svedberg aggregation-number
routes), 1 partially built with an honest, documented reason for the unshipped half
(Gouy-Chapman: the Grahame equation shipped, the self-consistent n_factor correction did not),
and 1 deferred without attempting for the same structural reason (FLM). Genuinely out-of-scope
items (Blankschtein MT theory, QSPR/ML HLB) remain flagged, not silently dropped. SurfactantKit's
test suite grew from 188 to 297 tests across this whole push, all passing.

Explicitly OUT OF SCOPE for this phase (flagged, not silently skipped): Blankschtein
molecular-thermodynamic theory (predictive from molecular structure alone -- a categorically
larger feature, would need its own free-energy-of-transfer/interface/packing/headgroup
sub-models, not a drop-in alternative formula) and QSPR/ML-based HLB prediction (Wang & Duan
2009, Luan et al. 2009 -- needs a training dataset; this is in scope for the separate SurfQSPR
paper/repo, not SurfactantKit).

---

## Also still open (older, lower priority than the above)

- **Rodenas, Maeda alternative-method validation against a second independent literature
  system** (beyond the single Azum et al. 2022 source both were verified against) -- would
  strengthen confidence but isn't blocking anything.
- **Second independent literature sources for categories D (geometry), F (dynamics), G
  (thermodynamics)** -- flagged as open in an earlier session, not touched this session,
  real Paper 3/SurfBench scope. Lower priority than the reporting-pipeline work above per
  the user's own explicit redirect this session, but not abandoned.
- **EOMMM's actual multi-point global-fit solver** -- DONE 2026-09-10, with an important
  disclosed caveat. Built `eommm_global_fit()` in `mixed_micelle.py`: fits W12, W21, and every
  point's x1 simultaneously via a genuine (n+2)-parameter nonlinear least-squares, using the
  same mass-balance equations `solve_rubingh_x` is built from, generalized to
  `asymmetric_margules_activity_coefficients`. **This is a first-principles construction, NOT
  a verified transcription of Schulz & Durand 2016's exact published procedure (their Eq.
  3.2 / SI Point 2.3)** -- the primary source is paywalled on ScienceDirect (same block hit
  for the HLD paper), and a related open-access companion paper (Serafini et al.,
  arXiv:1806.09721, a real TX100-DTAB EOMMM application reporting W12=+4.04 kBT,
  W21=-14.02 kBT) explicitly defers the exact procedure to its own SI, unavailable in the
  fetched copy. Flagged this gap to the user via AskUserQuestion rather than guess; user's
  explicit choice: "Build from first principles, clearly labeled." Real numerical findings
  made while validating (both documented in the function's own docstring, not silently
  patched): (1) a single point's mass-balance residual in x1 is NOT always unimodal for fixed
  W12/W21 -- directly analogous to solve_rubingh_x's own already-documented multi-root
  behavior at large beta -- fixed via a coarse grid scan before the golden-section refine
  (`_minimize_with_grid_refine`); (2) the OUTER joint objective is consequently non-convex too,
  so a single coordinate-descent run can converge to a local optimum -- fixed via multi-start
  (5 starting points, keep the lowest-SSE result). This pre-existing "Also still open" note
  (above, before this fix) already flagged that Schulz & Durand's own real implementation uses
  GAMS/BARON -- a global (not local) optimizer -- which independently corroborates that this
  problem genuinely needs global-search technique, consistent with what the multi-start fix
  above found empirically, even without access to their exact equations. Validated via
  mathematically-guaranteed round trips only (a closed-form construction that satisfies BOTH
  mass-balance equations exactly for ANY chosen W12/W21/x1, not just RST's symmetric case) plus
  a cross-check that the symmetric limit (W12=W21) recovers the same beta
  `solve_rubingh_x`/`rubingh_beta` find on identical data. 3 new tests in
  `tests/test_mixed_micelle.py` (round-trip, RST-reduction cross-check, bad-input rejection) --
  slow by the nature of the problem (a few seconds each, genuine multi-parameter optimization,
  no numpy/scipy) but kept to n=3 points to bound runtime; 48/48 passing in that file. Wired
  into `__init__.py` and `mcp_server.py` (new `eommm_global_fit` tool, which documents the
  same first-principles disclosure and expected slowness), registered + end-to-end tested in
  `test_mcp_server.py`. Full suite: 220/220 passing.
- **Ohshima's electrostatics gap for the ionic-dissociation `r` parameter** -- see
  `benchmark/METHOD_ALTERNATIVES_LITERATURE_REVIEW.md` for full detail, not repeated here.

---

## Standing reminder for whoever resumes this (from the user, 2026-09-08)

Once ALL work across the whole PhD project (not just SurfactantKit) is complete, the user
wants a comprehensive account compiled of every failure faced, every workaround attempted,
and the final success achieved, for their thesis defense presentation in front of their
guide and colleague. Saved in memory as `project_thesis_defense_documentation.md` -- **not
needed yet**, but don't lose track of it. The raw material for this already exists in
`PROJECT_STATUS.md`, each paper's own `ROADMAP.md`, `literature_validation_notes.md`, and
`METHOD_ALTERNATIVES_LITERATURE_REVIEW.md` -- when asked, it's a compilation job from real
existing records, not new research.

## MicelleMD compute status as of 2026-09-08 (for context, not this repo's concern)

Both Ubuntu-hosted runs finished cleanly the same day this file was written:
- CG production (150 gemini-analog surfactants, corrected physics): reached its full 200 ns
  target exactly, no crash. Final: largest aggregate 22 molecules, Rg ~12.7 nm (this Rg
  value looks large/loose for that aggregation number across the whole back half of the
  run -- flagged as worth a second look, not diagnosed).
- Atomistic extension (100 SDS anions): finished cleanly at 40 ns total simulated time
  (original 20 ns + this 20 ns extension). Final aggregation number not yet computed --
  needs `atomistic_analysis.py` run against the finished trajectory, a clean next step.
See MicelleMD's own `ROADMAP.md` for the authoritative record of this work.
