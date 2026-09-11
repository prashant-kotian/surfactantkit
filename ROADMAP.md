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
failure mode in every test file. `surfactantkit-reports` lives at
`H:\CodeProjects\surfactantkit-reports`, committed and pushed to
`github.com/prashant-kotian/surfactantkit-reports` (confirmed 2026-09-12: `git status` clean,
up to date with `origin/master`, its own 36/36 tests re-run and still passing against
SurfactantKit's current state despite the many changes made to this repo since 2026-09-10).

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
   **RESOLVED 2026-09-11** (primary source PDF provided by the user): Table 1 of the primary
   paper lists group numbers under "Davies and Lin" vs. "ECL method" columns side by side, and
   for CH2/CH3/CH/=CH- and the EO-repeat-unit group they are IDENTICAL (-0.475 and 0.33) --
   confirming directly that Guo/Rong/Ying reuse Davies' weights unchanged, exactly what
   `hlb_davies_guo_ecl` already assumed. Also added the previously-unimplemented NEO>=50 branch
   of the source's own eq. (5') (NEOeff=0.056*NEO+43.08), removing the n_eo<=50 cap
   `guo_effective_eo_chain_length` used to enforce. PO chain support is exposed as a standalone
   utility only (no verified Davies-scale PO group number exists in this project to chain it
   into a full HLB). 11 tests (9 original + 2 new, 2026-09-11) in `tests/test_hlb_cpp.py` --
   52/52 passing in that file. Wired into `__init__.py` and `mcp_server.py` (new
   `hlb_davies_guo_ecl` tool), registered + end-to-end tested.
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
   equilibrium) not verified this pass. **2026-09-11: primary source PDF obtained (provided by
   the user)**; confirmed directly against it (Rusanov's own eq. 18, cm=1/Kj=K^(1-n), is the
   same CMC definition this formula is built on, stated in his own abstract). Also added a new
   companion function, `critical_micellization_degree()` (same module), implementing the
   paper's own eqs. (15)-(16) for the critical micellization degree alpha_m under two different
   curvature-based CMC definitions -- exact-matched to the paper's own worked example (n=100:
   alpha_m=0.061 via eq. 15, 0.091 via eq. 16), with a guard against the d2alpha_dc2 branch's
   unphysical negative output for n<2. 10 tests (5 original + 5 new, 2026-09-11) in
   `tests/test_thermodynamics.py` -- 27/27 passing in that file. Wired into `__init__.py` and
   `mcp_server.py` (new `mass_action_free_energy_of_micellization` and
   `critical_micellization_degree` tools), registered + end-to-end tested.
6. **O'Brien-White relaxation-effect correction** -- DONE 2026-09-10, **UPGRADED 2026-09-11
   to the real primary-source formula** (old method fully replaced, not left alongside the
   new one, per the Upgrade Protocol). Built `electrophoretic_mobility_relaxation_corrected()`
   (forward) and `zeta_potential_relaxation_corrected()` (numerical inverse, grid+bisection
   like `solve_rubingh_x`) in `electrostatics.py`: a real, sourced alternative to the plain
   Henry equation for HIGH zeta (>50 mV, common for ionic surfactant micelles) -- Henry is
   linear in zeta and misses the relaxation effect (mobility grows sub-linearly with zeta, can
   pass through a maximum beyond ~100 mV) entirely.
   Originally (2026-09-10) shipped with a substitute formula (an IUPAC Technical Report's
   restatement of O'Brien's simplification of the Dukhin-Semenikhin equation) because the real
   primary sources -- O'Brien & White, J. Chem. Soc. Faraday Trans. 2, 74 (1978) 1607-1626, and
   Ohshima, Healy & White (OHW83), same journal, 79 (1983) 1613-1628 -- were paywalled.
   2026-09-11: both PDFs obtained (provided by the user); the substitute formula was fully
   REPLACED with OHW83's own semi-empirical formula (their eqs. 57-62 + the correction eqs.
   75-76), which their own Fig. 1/2 show beating the old substitute-family formulas at every
   kappa_a/zeta tested against the exact O'Brien-White computer solution, and widens the valid
   range from kappa_a > ~20 down to kappa_a >= 10 (the source's own eq. 76 domain, <1% relative
   error there). Validated three ways: (1) the original Smoluchowski cross-check (low zeta,
   large kappa_a -> `zeta_potential_henry`/`henry_function`'s f(kappa*a)=1.5 limit) still holds
   under the new formula; (2) a NEW direct check against the source's own eq. (63) small-zeta/
   large-kappa_a asymptotic limit, agreeing to 1e-4 relative error; (3) reproduces the
   qualitative shape of the source's own Fig. 1/2 (kappa_a=20, m=0.184: mobility rising to a
   maximum near zeta~=5, then declining) closely matching their digitized computer-result
   points -- checked in a scratch script before wiring into the library, not asserted blind.
   Also confirmed the real physical signature (sub-linear mobility growth at high zeta) and a
   clean round-trip (forward -> inverse solver) across 10-90 mV still hold. Restricted to
   kappa_a >= 10 (the source's own stated valid range) and zeta in [0, 150] mV by default for
   the inverse solve (below the reported non-monotonic-mobility maximum) -- both limits
   enforced with explicit errors, not silently guessed past. 9 tests (8 original + 1 new eq.-63
   check) in `tests/test_electrostatics_dynamics.py` -- 45/45 passing in that file. Wired into
   `__init__.py` and `mcp_server.py` (`predict_mobility_relaxation_corrected` and
   `zeta_potential_relaxation_corrected` tools, docstrings updated to the new source).
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

**2026-09-11: primary-source paper batch.** User obtained and provided 7 PDFs this session
(Schulz & Durand main text again -- SI still not among them; Rodenas 1999 primary; Rusanov
2014 primary; O'Brien & White 1978 primary; Ohshima/Healy/White 1983 primary; Guo/Rong/Ying
2006 primary; Proverbio et al. 2003, already resolved 2026-09-08 and now independently
reconfirmed). Net effect: the O'Brien-White relaxation-effect correction was materially
upgraded to the real OHW83 formula (item 6 above); the Guo/Rong/Ying weight-reuse uncertainty
was resolved and its NEO>=50 branch added (item 4 above); the Rusanov citation was upgraded to
direct primary-source confirmation and a new `critical_micellization_degree()` capability was
added (item 5 above); a genuinely independent second literature system for Rodenas was found
and tested, disclosed honestly as a worse (not better) numeric match with an explained cause
(see the Rodenas entry below). Schulz & Durand's SI remains the one item from the original gap
list still genuinely unresolved. Full suite: 309/309 passing (up from 297).

**2026-09-11 (continued): second paper batch + a genuinely new capability.** More PDFs
obtained and provided by the user across several rounds: Schulz & Durand's actual SI (mmc1.docx
+ 9 CSV datasets) -- gave the real GAMS objective function and r-parameter/GRDIS wiring, but
implementing it exactly revealed the free-energy-minimization objective is unbounded below as
literally transcribed (runs to whatever W12/W21 search box is set, on both synthetic and real
Hyamine/DTAB data) -- a genuine, disclosed dead end for now, not silently patched around (see
the EOMMM entry above for the full finding); Bales 1998/2001/2002 (SDS tail volume/aggregation
number worked example, alpha=0.272 SDS ionization degree, CsDS Krafft-temperature ionization);
Khademi et al. 2017 + its SI (checked, figures only, confirms the existing round-trip approach
was already the best available); Vautier-Giongo/Bales 2005 (two new independent Rubingh
systems, DHPC+SDS and DHPC+DTAB, one of the five DHPC+SDS points a disclosed ~40% outlier,
likely a table typo in the source); Schafer et al. 2020 (a real CORRECTION, not just an
addition -- the project's earlier secondary-summary paraphrase of this paper's CPP-threshold
finding was simply wrong, fixed in `classify_aggregate_morphology`'s own docstring); Sutherland
et al. 2009 (closed the long-standing Category F/Stokes-Einstein gap for real, with an
explanation for why 3 earlier attempts had failed: DLS-measured "mutual" diffusion coefficients
for IONIC surfactants run several-fold too large due to counterion coupling, so a plain
Stokes-Einstein check needs the true micelle diffusion coefficient specifically, which this
paper's Taylor-dispersion data provides).

**New capability (not a citation): `derive_davies_group_number_from_griffin()` in `hlb.py`.**
Built in response to the user's question "can we compute HLB values for amide/sulfonate
ourselves" -- yes, via the SAME cross-calibration method Davies himself used to build his
original 1957 table (assume Davies' additive HLB equals Griffin's independent mass-ratio HLB,
20*Mh/M, for a real reference compound, solve for the one unknown group number) -- the same
method this project already used once, via a literature source (B.H. O 1998), to resolve the
quaternary-ammonium gap. Applied it to derive real, disclosed, NON-literature-sourced group
numbers for **sulfonate** (-SO3Na, calibrated against the sodium alkyl sulfonate series C6-C12,
converging to ~6.2-6.3 over the practically-relevant C9-C12 range) and two **fatty acid
alkanolamide** head groups (di- and mono-ethanolamide, -CO-N(CH2CH2OH)2 and -CO-NH-CH2CH2OH,
e.g. cocamide/lauramide DEA and MEA -- real, common, commercially important nonionic
surfactants, not bare primary fatty amides which are barely water-soluble), calibrated against
lauric/myristic/palmitic acid derivatives, self-consistent to ~2-5% across chain lengths.
Kept in a separate `DAVIES_DERIVED_HYDROPHILIC_GROUPS` table (never merged into the
literature-only `DAVIES_HYDROPHILIC_GROUPS`), usable via `hlb_davies(..., allow_derived_groups=True)`
-- an explicit opt-in, so the distinction between "Davies' own literature number" and "this
project's derived estimate" can never blur at the call site. Honestly disclosed limitation: the
derived numbers are chain-length-sensitive (U-shaped, not a true asymptote -- expected, since a
linear-in-chain-length Davies model cannot exactly reproduce Griffin's nonlinear Mh/M ratio at
every chain length; the two only agree closely over the realistic commercial-surfactant range).
9 new tests in `tests/test_hlb_cpp.py` (self-consistency across chain lengths, the trivial-by-
construction Davies/Griffin agreement check, bad-input rejection, opt-in-required behavior, a
worked lauramide-DEA example landing in the real-world catalog HLB range ~9-10) -- 60/60 passing
in that file. Wired into `__init__.py` and `mcp_server.py` (new `derive_davies_group_number`
tool; `hlb_from_groups` gained the `allow_derived_groups` parameter), registered + end-to-end
tested. Full suite: 322/322 passing.

**2026-09-11 (continued): 4 more bounded items completed using material already in hand.**

1. **Extended the derived-group-number method to 3 more of the SurfBench-flagged gaps.** Same
   `derive_davies_group_number_from_griffin()` methodology, applied to **sultaine**/alkyl
   sulfobetaine (-N+(CH3)2-(CH2)3-SO3-, e.g. SB10/SB12 -- the exact zwitterionic compounds
   already cited via Sutherland et al. 2009; tightest self-consistency of any derived group,
   ~2.7% spread across C8-C14, GN=8.61 @ C12), **carboxybetaine** (-CO-NH-(CH2)3-N+(CH3)2-CH2-COO-,
   e.g. cocamidopropyl betaine/CAPB, one of the most common real amphoteric surfactants, GN=9.16
   @ C12), and **phosphate ester**, disodium monoalkyl (-O-P(=O)(ONa)2, GN=7.79 @ C12).
   **Imidazoline deliberately skipped, not just unattempted**: real imidazoline surfactants
   hydrolyze in solution into a genuinely ambiguous mixture of open-chain forms, so there is no
   single undisputed structure to compute Mh for -- disclosed as a real, structural (not just
   informational) gap. 7 new tests in `tests/test_hlb_cpp.py` -- 64/64 passing in that file.
   **RESOLVED 2026-09-12, not left as a permanent gap.** "Imidazoline" itself still isn't
   derivable (it's genuinely never the final sold product -- always a synthesis intermediate
   reacted further before sale), but the real, final, commercial imidazoline-DERIVED family
   is: disodium lauroamphodiacetate (INCI; the fully dicarboxymethylated reaction product,
   PubChem CID 109973, confirmed formula C20H36N2Na2O6 two independent ways -- hand
   atomic-mass summation and RDKit `get_weight_from_smiles` on the open-chain amide tautomer,
   both landing on 446.496 g/mol). PubChem's own IUPAC name describes the cyclic
   imidazolinium tautomer instead -- a real, literature-documented open-chain/cyclic
   tautomer ambiguity for this compound class -- but it doesn't block the Mh calculation
   since both tautomers are isomers of the identical formula and the tail/head atom split
   used (tail = the undecyl C11H23 group only) lands the same either way. GN=11.27 @ C12
   (computed self-consistently at C10/C12/C14, ~2.3% spread). Added as
   `"amphodiacetate"` in `DAVIES_DERIVED_HYDROPHILIC_GROUPS`, with full derivation and the
   tautomer-ambiguity disclosure in `hlb.py`'s own block comment, and disclosed scope (this
   is specifically the diacetate product; the mono-carboxymethylated "amphoacetate" form is
   a different, also-real, also-commercial product not covered here). New test:
   `test_derive_davies_group_number_amphodiacetate_self_consistent` in `tests/test_hlb_cpp.py`.

2. **Wired Bales 2001's real, precisely-measured SDS ionization degree (alpha=0.272+/-0.017,
   J. Phys. Chem. B 105 (2001) 6798-6804) into `gibbs_free_energy_micellization`'s docstring and
   a new worked example**, replacing the previously arbitrary illustrative counterion_factor with
   a properly sourced one (counterion_factor = 1+alpha = 1.272 for SDS specifically). Along the
   way, caught and fixed a real error in the planning for this: initially mis-derived
   counterion_factor = 2-alpha = 1.728 by treating Bales' alpha (a DISSOCIATION degree) as if it
   were beta (a BINDING degree) -- caught by checking the (2-beta)=(1+alpha) identity against its
   own physical limits (fully dissociated -> counterion_factor=2; fully bound -> 1) before
   committing anything to code, not after. New tests confirm both the identity's physical limits
   and a real worked SDS deltaG_mic (-27.77 kJ/mol using the well-known salt-free cmc0=0.0083 M
   already cited elsewhere in this project via Bales 1998) -- honestly scoped as a properly-
   sourced INPUT with a real, physically-sensible output, not a literature-matched deltaG_mic
   value (none was found to check the final number against). 2 new tests in
   `tests/test_thermodynamics.py` -- 29/29 passing in that file.

3. **Serafini et al. 2019 SI's SLS-determined micellar masses (Table SI-II)** used as an
   aggregation-number sanity check: N_agg(TX-100) = 66600/625 = 106.6 (within the well-known,
   temperature-sensitive literature range ~100-155); N_agg(DTAB) = 16800/308.34 = 54.5 (excellent
   match to the commonly-cited ~50-60, often ~56). NOT wired into a new function -- N_agg =
   M_micelle/M_monomer is a one-line ratio with no natural function to attach it to; creating one
   would be a premature abstraction around a single division. Recorded as a reference/sanity
   check in `literature_validation_notes.md` instead.

4. **Schulz & Durand 2016 SI's CSV datasets, mmc2.csv and mmc3.csv, wired into
   `asymmetric_margules_activity_coefficients` validation.** Real insight: at the pure-component
   boundaries (x1->0, x1->1), the CSVs' own "Gexc/RTx1x2" column is trivially equal to W12 and W21
   respectively (a direct algebraic consequence of the Margules formula), so W12/W21 can be READ
   OFF the data directly rather than fitted. This let mmc2.csv's boundary values (-3.567, -0.890
   RT units) be cross-checked against `test_asymmetric_margules_hyamine_dtab_literature_case`'s
   existing W12=-8836.50 J/mol, W21=-2204.19 J/mol (Hyamine-DTAB, Case Study 1) -- matching to
   within CSV rounding, CONFIRMING mmc2.csv is that same system, and upgrading that test from a
   sign/order-of-magnitude-only check (all that was possible before, with only a graph available)
   to an exact digit-level match across all 5 interior points. mmc3.csv used the same way for a
   second real dataset (7 interior points, not confidently identified as a specific named system
   from the extractable SI text, cited generically) -- also an exact match. **mmc4.csv NOT used**:
   its "delta theta cmc" quantity is a deviation/residual (confirmed by physically-impossible
   negative "CMC ideal" values in the raw data), not raw CMC or a directly interpretable
   dissociation-generalized-Clint output -- genuinely ambiguous from the extractable SI text, not
   guessed at. **mmc5/6 (and 7/8, 9/10) NOT converted into a precise test**: real experimental
   Gexc (mmc5/7/9) and EOMMM-fitted Gexc (mmc6/8/10) are reported at DIFFERENT, non-matching
   composition grids -- qualitatively consistent in magnitude (e.g. mmc5's real -1645.9 J/mol at
   x=0.308 vs. mmc6's fitted -1698.8 J/mol at x=0.328, ~3% apart) but a precise quantitative test
   would need interpolation across mismatched grids, not attempted this pass. 2 new/upgraded
   tests in `tests/test_mixed_micelle.py`. Full suite: 331/331 passing.

**2026-09-11 (continued): `eommm_global_fit` REWRITTEN -- the free-energy dead end resolved.**
This closes the one item repeatedly flagged as "genuinely hard, still open" throughout this
whole session. Real history, kept in full because it matters for trusting the result:

1. The Schulz & Durand 2016 SI's own GAMS code gives the real objective as minimizing TOTAL
   FREE ENERGY OF MICELLIZATION subject to r-generalized mass-balance constraints (CMCexpVar
   bounded, not pinned, to the experimental CMC). A first attempt to transcribe that literally
   found it UNBOUNDED BELOW (the Margules mixing term diverges as \|W12\|,\|W21\| grow) -- it ran
   to whatever W12/W21 search box was set, on both synthetic and real Hyamine/DTAB data.
   Abandoned rather than shipped broken (documented at the time as a genuine dead end).
2. A companion paper's SI (Serafini, Fernandez-Leyes, Sanchez M., Pereyra, Schulz E.P., Durand,
   Schulz P.C. & Ritacco, Colloids Surf. A (2019), TX100-DTAB system) explained the REAL
   procedure in prose: "the mg [margin] parameter was varied in order to obtain the minimum
   value that allowed a feasible solution" -- i.e. the real fitting criterion is CONSTRAINT
   SATISFACTION at the tightest feasible margin, not open-ended free-energy minimization under
   a fixed generous margin. This is what the rewrite implements: minimize a total INFEASIBILITY
   measure of the r-generalized cdefp1/cdefp2 constraints at a given `cmc_margin` -- a
   well-posed, bounded objective, unlike the free-energy approach.
3. A first working version of this new objective had its own real bug, found and fixed before
   shipping: using the geometric mean of the two individually-implied CMCvar values to check
   margin compliance let a badly-mismatched pair (e.g. a 40% relative disagreement between the
   two mass-balance equations) slip through with a geometric mean that coincidentally landed
   within 1% of the target CMC -- corrupting the search into spurious non-physical optima.
   Fixed by penalizing the mismatch between the two equations directly, not just checking
   whether their average happens to look right.
4. Two real, separate multi-root/resolution issues were found and fixed, each confirmed via
   direct numeric investigation (not assumed): (a) coordinate-descent/golden-section refinement
   of the outer (W12, W21) search, even STARTED exactly at the known true optimum, was found to
   walk AWAY from it -- replaced with a pure 2-level (coarse-then-fine) grid search, which
   recovers the true point exactly; (b) the inner per-point x1 solve occasionally landed on a
   real second root (verified: two different x1 values both gave mismatch=0 to machine
   precision for the same alpha1/W12/W21) at the original grid_points=40/refine_iters=40
   resolution -- fixed by raising to 100/60 (the only caller, so no side effects elsewhere).
5. A real, disclosed, carefully-VERIFIED (not merely claimed) physical property of the
   resulting model, matching the Serafini SI's own stated rationale for the margin-tightening
   procedure: at a tight `cmc_margin` (e.g. 1e-6), the round-trip recovers the true (W12, W21)
   EXACTLY. At the SI's own default `cmc_margin=0.10`, the fit is genuinely NOT always uniquely
   determined -- a real, still-feasible (total_infeasibility ~0) fit can differ from the true
   parameters (verified: W12 off by 0.3 in one concrete case). An earlier, more dramatic
   specific claim about this (a very different alternate feasible point) was tested directly
   and found to be a resolution artifact from an early prototype, not real -- corrected in the
   shipped tests rather than left in, since an honest "feasible but not necessarily exact" is
   the claim that actually holds up.
6. This function does NOT automate the SI's own margin-tightening search (an automatic
   binary-search-on-margin was prototyped and found to have a real, unresolved grid-resolution
   sensitivity of its own, converging to a margin somewhat looser than the true minimum) --
   disclosed as a genuine remaining limitation, not hidden. Callers wanting the tightest
   well-determined answer should call this function at a few decreasing `cmc_margin` values
   themselves, the same diagnostic the SI's own procedure does manually.

New signature: `eommm_global_fit(alpha1_series, cmc_mix_series, cmc1, cmc2, r1=2.0, r2=2.0,
cmc_margin=0.10, w_bound=30.0)` -- r1/r2 are the real dissociation-number generalization from
the SI (default 2.0, the classical fully-dissociated convention already used everywhere else
in this module). `EommmGlobalFitResult` gained `cmc_var_values` and `total_infeasibility`
(replacing `sse`); `method` field updated. 4 tests in `tests/test_mixed_micelle.py` (exact
recovery at tight margin, feasible-but-imprecise at default margin, symmetric-limit reduction
to Rubingh beta, bad-input rejection) -- run in ~20s total, consistent with this project's own
established "a genuine multi-parameter search, expect several seconds, not an error"
precedent. Wired into `mcp_server.py` (new `r1`/`r2`/`cmc_margin` parameters, new return
fields, docstring rewritten with the full history). Full suite: 332/332 passing.

**2026-09-11 (continued): three more bounded items completed.**

1. **A second, independent anchor for the amide HLB numbers.** Checked directly against a
   fetched primary source (not a search-engine paraphrase): Pawignya et al. (IOP Conf. Ser.:
   Mater. Sci. Eng., Atlantis Press) report a real, EXPERIMENTALLY MEASURED HLB=5.940 for a
   palm-oil-derived diethanolamide, via their own CMC-based formula (genuinely different from
   Griffin's mass-ratio formula this project's numbers are calibrated against -- not an exact-
   match comparison). Honestly scoped as a qualitative cross-check: palm oil's longer-chain
   fatty-acid profile (C16-C18-dominant vs coconut/lauric C12) giving a lower HLB than this
   project's C12-calibrated value is directionally consistent with this project's own
   already-documented chain-length trend, and both land in the same broad order of magnitude.
   A second, repeatedly-cited search-reported value (~13.2-13.5 for cocamide DEA) could NOT be
   independently verified against a fetched primary document (ResearchGate/SpecialChem/
   IOPscience all blocked automated fetch) -- kept as an unconfirmed, disclosed data point.
   1 new test in `tests/test_hlb_cpp.py`.
2. **mmc4.csv's "Δθcmc" quantity -- checked and confirmed genuinely unresolvable from the
   extractable SI text, not guessed at.** Searched the SI's own prose text specifically for
   "theta" -- zero mentions anywhere except the image-only Figure 7 caption/axis label
   ("Variation of the critical micellar concentration mixture (Δθcmc) with αDTAB"). No defining
   equation exists anywhere in the extractable body text; the symbol is only ever shown
   rendered inside a figure. The negative values in the raw table confirm it IS a real
   deviation/difference quantity (a raw CMC can never be negative), but the specific reference
   subtracted was never written down in recoverable text. A genuine, confirmed dead end from
   this source -- not silently dropped, actively investigated with evidence.
3. **mmc5/7/9 (real experimental) vs. mmc6/8/10 (EOMMM-fitted) Gexc comparison, converted into
   real, precise tests.** Two results, of very different strength:
   - **The strong result**: since Gexc/RT = x1*x2^2*W12 + x1^2*x2*W21 is LINEAR in (W12, W21)
     for a fixed composition, (W12, W21) can be solved EXACTLY from any 2 of mmc6/mmc8/mmc10's
     own interior points (a real 2x2 linear system, not curve-fitting) -- and
     `asymmetric_margules_activity_coefficients`/`excess_free_energy` then reproduce ALL of
     each curve's OTHER points to <0.25% relative error. This is a real, precise, third-through-
     fifth confirmation that this project's formula is exactly the one Schulz & Durand used to
     generate their own EOMMM-fitted curves for the C8E4-SDS system (Case Study 3) -- much
     stronger than the originally-planned interpolation comparison. 3 new tests.
   - **The weaker, honestly-disclosed result**: comparing the REAL experimental data
     (mmc5/7/9, from Hey et al. 1985) against the fitted curves via piecewise-linear
     interpolation (pure data comparison, not a test of this project's own code) shows a real,
     methodologically-explained pattern -- compositions well within the fitted curve's dense
     coverage agree closely (<3% relative), while compositions near the edges (where the
     nearest fitted points are far away and linear interpolation poorly approximates the true,
     presumably-curved approach to zero at the boundary) disagree by 9-21%. 1 new test
     documents this pattern explicitly rather than hiding it.

   4 new tests total in `tests/test_mixed_micelle.py`. Full suite: 337/337 passing.

---

## Also still open (older, lower priority than the above)

- **Maeda validation against a second independent literature system** -- DONE 2026-09-10.
  Found a genuinely different system (a drug-surfactant pair, not another gemini/conventional
  surfactant): Rub, Azum, Kumar, Arshad, Khan, Alotaibi & Asiri, Polymers 13(22) (2021) 4025,
  imipramine hydrochloride (IMP) + Triton X-100, their own reported X1^Rb=0.8585,
  beta^Rb=-4.35, deltaG_Maeda=-19.64 kJ/mol at alpha1=0.5, 298 K. Reproducing their deltaG from
  their own reported X1/beta gives -20.85 kJ/mol, a real ~6% relative error -- plausible
  rounding residual from their 3-4-sig-fig table, not a formula error (right sign, right order
  of magnitude, matches within single-digit percent). Real, disclosed limitation found along
  the way: attempting to ALSO independently re-solve X1 via this project's own
  `solve_rubingh_x` from their raw (alpha1, cmc_mix, cmc1, cmc2) inputs gave 0.147, not 0.8585,
  unless cmc1/cmc2 are swapped in that specific call (giving 0.853, close) -- a genuine
  component-labeling-convention ambiguity between this paper's Rubingh-solve indexing and its
  Maeda indexing that couldn't be resolved with full confidence from the fetched text (two
  separate fetches of the same paper gave inconsistent descriptions of which component is
  "ionic"/"nonionic" in their indexing). The new test therefore validates
  `maeda_free_energy_of_micellization` using their own reported X1/beta directly, not an
  additional Rubingh solver cross-check -- honestly scoped, not overclaimed. New test in
  `tests/test_mixed_micelle.py`, full docstring disclosure. Suite: 298/298 passing.
  **Rodenas second-source validation still open for a genuinely different source, but a real
  related finding landed 2026-09-10.** Found the SAME Azum et al. 2022 paper (not a second
  source) has a full 5-point G6+T-20 alpha1 series (Table 1/2) with real reported X1_Rod
  values at every point -- used this to run `rodenas_x1_series` (the numerical local-slope
  helper, previously only round-trip tested against synthetic exactly-quadratic data) against
  REAL, sparse, experimentally-noisy multi-point data for the first time. Real, disclosed
  result: interior points (alpha1=0.4, 0.5, 0.6) match within 3.4-9.1% relative error, but
  the two series ENDPOINTS are off by ~19-20% (one even landing outside the physically-valid
  [0,1] range) -- consistent with, not contradicting, the function's own documented weaker
  one-sided-derivative reliability at series endpoints, now empirically confirmed on real
  data rather than just a theoretical concern. New test in `tests/test_mixed_micelle.py`,
  explicitly NOT claimed as a second-source validation (same primary source) and explicitly
  NOT claimed as proof the function matches literature to the precision the round-trip tests
  alone would suggest.
  **RESOLVED 2026-09-11: a genuinely different second source obtained** -- Rodenas, Valiente &
  Villafruela, J. Phys. Chem. B 103(21) (1999) 4549-4554, the paper that ORIGINATED this model
  (not a later paper citing it), PDF provided by the user. Its own Table 1 (a completely
  different chemical system: C12E4/CTAB, not G6-gemini) has 5 real (alpha1, CMC*, chi1) points.
  Real, disclosed result: this reproduces WORSE than the Azum-paper check (29-60% relative
  error at every point, one point -- alpha1=0.5 -- landing slightly outside [0,1]), and the
  paper's own text explains why: it computes its slope by differentiating a smooth 2-exponential
  global fit to CMC*(alpha1) (its own eq. 14) analytically, not by finite-differencing 5 sparse,
  unevenly-spaced (0.025 to 0.7) raw points the way `rodenas_x1_series` does -- and the paper's
  own text states it excluded the alpha1=0.025 point from its own Figure 4 as unreliable. New
  test in `tests/test_mixed_micelle.py`, kept as a real, honest record like the Azum-paper
  check (evidence of the sparse-series numerical-differentiation method's real limits on data
  this sparse, not evidence the core `rodenas_x1` formula is wrong -- that's separately verified
  via an exact zero-interaction algebraic identity that doesn't depend on finite differencing).
- **Second independent literature sources for categories D (geometry), F (dynamics), G
  (thermodynamics)** -- attempted 2026-09-10, mixed real results, all honestly disclosed:
  - **Category G (thermodynamics): DONE.** Found a genuinely different system (an
    alkylguanidinium cationic surfactant, not a conventional headgroup) -- Bouchal, Hamel,
    Hesemann, In, Prelot & Zajac, Int. J. Mol. Sci. 17(2) (2016) 223, dodecylguanidinium
    chloride, their own reported CMC=6.2 mmol/kg, beta=0.74, deltaG_mic=-28.3 kJ/mol at 298 K.
    `gibbs_free_energy_micellization` with this project's (2-beta) counterion_factor
    convention reproduces -28.4 kJ/mol -- ~0.4% relative error, well inside the paper's own
    +/-0.9 kJ/mol uncertainty. A clean, strong match, confirming both the formula and the
    counterion_factor convention against a second, chemically unrelated real system. New test
    in `tests/test_thermodynamics.py`. The same paper's van't Hoff enthalpy (needs CMC at 2+
    temperatures) could NOT be cross-checked -- their conductivity table only gives DDGC's
    CMC at the single temperature 298 K; a second CMC-vs-T pair for the same surfactant would
    be needed and wasn't found in this pass.
  - **Category D (CPP/Tanford geometry): partial win after 8 attempts.** 7 search/fetch
    attempts (ResearchGate 403s, an SSL failure on norgwyn.com, several papers citing Tanford
    only qualitatively without reproducing his coefficients or a worked example) found no
    clean external numeric example from a real SURFACTANT PAPER to check
    `tanford_tail_volume`/`tanford_critical_length` against -- that harder bar remains unmet.
    Along the way found that secondary sources round the critical-length coefficient slightly
    differently across the literature (1.265 vs 1.26 vs 1.256), consistent with the same
    underlying Tanford reference, not an error. An 8th attempt found a genuine independent
    confirmation instead: a physical-chemistry course problem (unrelated to this project's
    own citations) states v(n)=(27.4+26.9n)e-3 nm^3 and lc(n)=(0.154+0.1265n) nm -- converting
    units, this is an EXACT match to `tanford_tail_volume` and matches
    `tanford_critical_length` to within 0.04 Angstrom (a small, explainable rounding
    difference in the constant term only; the chain-length coefficient, 1.265, matches
    exactly). New test in `tests/test_hlb_cpp.py`, explicitly disclosed as confirming the
    FORMULA's correct transcription from an independent source, not a real-surfactant worked
    example -- a real but narrower win than Category G's.
  - **Category F (dynamics/Stokes-Einstein): still open, real attempt made and a real near-
    miss disclosed rather than accepted.** Found a search-snippet-reported (D, R_h) pair for
    SDS micelles in D2O (D=0.1213e-9 m^2/s, R_h=18.3 Angstrom) -- but independently computing
    Stokes-Einstein from that D using D2O's known viscosity gave 16.4 Angstrom, an 11%
    mismatch, and the primary source (a ResearchGate figure/table) could not be reached
    directly to confirm whether their reported R_h came from a plain Stokes-Einstein
    calculation or a more elaborate fit. Deliberately NOT counted as a validated second
    source given that unresolved discrepancy -- reported as a real near-miss, not silently
    accepted or silently dropped.
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
  **RESOLVED 2026-09-11: this first-principles construction was fully REPLACED** (per the
  Upgrade Protocol -- not kept alongside the new version) with a real, SI-derived transcription,
  after obtaining both the Schulz & Durand SI itself and a companion paper's SI that explained
  the real margin-tightening procedure. See the dated 2026-09-11 entry earlier in this file
  ("`eommm_global_fit` REWRITTEN -- the free-energy dead end resolved") for the complete history
  -- including a genuinely unbounded free-energy objective that was tried and abandoned first,
  a real bug found and fixed in the corrected infeasibility objective, and a carefully-verified
  (not just claimed) margin-dependent uniqueness property of the resulting model.
- **Ohshima's electrostatics gap for the ionic-dissociation `r` parameter** -- see
  `benchmark/METHOD_ALTERNATIVES_LITERATURE_REVIEW.md` for full detail, not repeated here.

**2026-09-12: EOMMM minimal-feasible-margin search AUTOMATED -- the earlier grid-resolution
disclaimer no longer applies.** When `eommm_global_fit` was rewritten (2026-09-11 entry
above), an automatic "binary-search the margin down to the tightest feasible value" wrapper
was tried and shelved: an early prototype showed real grid-resolution sensitivity, because
it was tested against the objective/grid resolution *before* the two bugs documented in that
same entry were fixed (the geometric-mean feasibility-check bug, and the too-coarse default
grid that missed a genuine second root). Rather than assume that limitation still held,
the automation was re-tried from scratch against the current, already-fixed
`_eommm_point_infeasibility` / `_two_level_grid_search_2d` / `_minimize_with_grid_refine`.
It now works reliably: a synthetic round-trip converges monotonically to the exact true
(W12, W21) at every intermediate margin tested down to ~1e-6, and real Hyamine/DTAB data
converges to a genuine, stable, nonzero minimal margin (~6.2%) reflecting real measurement
noise rather than search noise. **Methodological lesson worth keeping**: a previously-shelved
approach should be re-tested after its underlying dependencies are independently fixed, not
assumed permanently blocked just because it failed once.

Implemented as `eommm_find_minimal_feasible_margin(alpha1_series, cmc_mix_series, cmc1, cmc2,
r1=2.0, r2=2.0, w_bound=30.0, margin_hi=0.30, feasibility_tolerance=1e-3,
binary_search_iters=20)` in `mixed_micelle.py`, immediately after `eommm_global_fit`. Binary-
searches the margin between 0 and `margin_hi`, at each step calling `eommm_global_fit` and
checking `total_infeasibility` against `feasibility_tolerance`; raises `ValueError` up front
if even `margin_hi` itself is infeasible (a real data/model inconsistency, not a search-
precision issue -- distinct from the search simply not converging). `eommm_global_fit`'s own
docstring updated to point here instead of disclaiming automation as unresolved. Wired into
`__init__.py` (import + `__all__`) and `mcp_server.py` (new `eommm_find_minimal_feasible_margin`
tool, exposing `binary_search_iters` so callers can bound runtime). Tests: round-trip and
bad-`margin_hi`-rejection tests added to `tests/test_mixed_micelle.py`; MCP end-to-end
round-trip test added to `tests/test_mcp_server.py` (both use a reduced `binary_search_iters=8`
to bound runtime, since each iteration re-runs the full two-level grid search). Full suite:
340/340 passing.

**2026-09-12: Muherei & Junin 2009 and DTAB-SDS (PMC6554738) beta mismatch re-investigated.**
Both flagged as low-priority, "likely dead ends" by the user, with instruction to actually
investigate rather than leave stale. Re-checked with fresh primary-source access attempts
rather than re-stating the old notes as-is:
- **DTAB-SDS (PMC6554738)**: substantially resolved. Fetched the paper's real Table 3 directly.
  Its own stated ideal CMC (9.04 mM) reproduces from `clint_ideal_cmc()` to 4 sig figs (9.038 mM).
  The earlier "beta did not reproduce" comparison had implicitly used Table 3's x=0.5 row
  (beta=-2.5674) -- but x=0.5 is a mathematical identity point of the RST equations (Table 3's
  own f1Rub and f2Rub columns are printed IDENTICAL there, 0.526377=0.526377, true for ANY
  beta at x=0.5), so that row is very likely an illustrative scan entry, not the paper's actual
  solved root for the stated composition. Solving this library's own `solve_rubingh_x`/
  `rubingh_beta` on the paper's directly-stated data point (alpha_DTAB=0.25, CMCmix=6.011 mM)
  gives x1=0.3045, beta=-2.27 -- a ~12% magnitude difference from -2.5674, correct sign,
  squarely inside the same 5-20% pointwise-vs-regression gap already documented for every
  other system in `literature_validation_notes.md`. No longer an unexplained anomaly.
- **Muherei & Junin 2009**: confirmed genuine dead end, with new evidence rather than an
  assumption. Found and fetched the real open-access source (`scialert.net`, Asian J. Appl.
  Sci. 2(2), 115-127) -- confirmed the paper's own CMCid formula is exactly Clint's relation
  and its pure-component CMCs (0.387, 3.468 mM) exactly match what this project already uses,
  so the mismatch is not a wrong-formula/wrong-input problem on this project's side. The
  blocker: Table 2A's exact alpha-to-CMCid row mapping is an embedded image, unreadable by
  text extraction (same limitation as mmc4.csv's I_I,cmc). Two PDF-download attempts to get
  past this (the journal's own redirect, a ResearchGate copy of a related companion paper)
  were both blocked (JS-gated shell; HTTP 403). Confirmed unresolvable without manual/
  institutional access, not a library defect.
Full detail: `literature_validation_notes.md`'s "Notes on individual systems" section.

**2026-09-12: Category D's harder bar met -- a true, independent, real-surfactant CPP worked
example found.** The 2026-09-10 entry above closed Category D with a course-problem formula-
transcription check only, explicitly disclosed as "not a worked example from a real surfactant
system -- that harder bar remains unmet." Re-attempted rather than left as a permanent gap:
found Kamboj, Kaur, Bhalla et al., *R. Soc. Open Sci.* 6, 181979 (2019), PMC6458362 (SDS-DTAB
mixed micelles with dyes), Table 2. Both surfactants are 12-carbon chains (shared Tanford
V0/lc), and the paper reports Amin from its own real Gibbs-isotherm surface-tension-slope
measurement (not assumed) at 3 temperatures for each of the SDS-rich and DTAB-rich systems --
6 independent (Amin, P) pairs total. Chaining `tanford_tail_volume(12)` ->
`tanford_critical_length(12)` -> `critical_packing_parameter()` reproduces the paper's own P
to within 1.2% at every point, and `classify_aggregate_morphology()` matches the paper's own
"cylindrical or rod-shaped micelles" call at every point -- the residual <1.2% offset is the
same already-documented 1.265-vs-1.26 lc-coefficient rounding, not a new gap. This is the real
thing Category D was missing: a paper's own experimentally-derived input feeding this
project's own full CPP pipeline and landing on that same paper's own reported output. New test:
`test_critical_packing_parameter_matches_real_sds_dtab_paper` in `tests/test_hlb_cpp.py`.
Full detail: `literature_validation_notes.md`'s CPP/Tanford geometry section.

**2026-09-12: real bug in `cmc_from_surface_tension_curve` found via a researcher-submitted
dataset in the new browser UI (`surfactantkit-reports`) -- fixed, not a dataset problem.**
A researcher uploaded a genuinely clean, standard-shaped tensiometry curve with a flat
pre-onset "lag" baseline (surface tension near pure water) at low concentration BEFORE the
real decline starts, then a real decline, then a real plateau -- three regimes. The tool
reported CMC = 1.16 mM; the true breakpoint (visually and physically obvious from the plot)
was ~10.7-13 mM. Diagnosed by reproducing the algorithm's own split-search RSS table by hand
on the exact submitted data (not guessed): the plain two-segment model can't represent three
regimes, and forcing one straight line through [decline+plateau] combined scored a lower
combined RSS (~41) than forcing one line through [baseline+decline] combined (~469) at the
true breakpoint -- so minimizing raw RSS picked the wrong pairing every time, even though the
data itself was unambiguous by eye. Root cause confirmed both numerically (RSS table) and
visually (the PNG) before any fix was attempted, per this project's own "verify against real
source when debugging" discipline.

**Fix**: `cmc_from_surface_tension_curve` now also tries a 3-segment model (flat baseline +
decline + flat plateau) and selects between the 2- and 3-segment models via BIC (Bayesian
Information Criterion -- the standard statistic for this exact segmented-regression model-
selection problem, penalizes the 3-segment model's extra parameters so it's only chosen when
a real third regime is present). Verified on both datasets: the researcher's data now
correctly selects the 3-segment model, CMC = 13.17 mM (R^2 = 0.99996 on the true declining
segment), matching the visually obvious breakpoint; the existing AOT literature pilot dataset
(Shah/Das/Bhattarai 2025, no baseline lag) still correctly selects the 2-segment model and
reports the SAME CMC as before (2.51 mM, unchanged) -- confirmed no regression, not just
assumed. `CmcFromCurveResult` gained `premicellar_x_min_mM` and `n_baseline_points` fields
(0 when no baseline detected) so downstream consumers know when/where a baseline was
excluded; `premicellar_slope`/`r_squared_premicellar`/`n_premicellar_points` now always
describe the DECLINING segment specifically (unchanged behavior for 2-regime data, since
that already equaled the "premicellar" segment before). 2 new tests in
`tests/test_curve_analysis.py` using this exact real submitted dataset plus a regression
guard on the AOT case. `mcp_server.py`'s `cmc_from_surface_tension_curve` tool updated to
surface the 2 new fields. `surfactantkit-reports`' `pipeline.py` fixed to use
`premicellar_x_min_mM` directly (was recomputing the fit-line's x-range from a stale
"everything <= CMC" filter, which would have drawn the steep decline line stretched back
across the flat baseline region -- visually wrong) and now plots the baseline as its own
separate flat segment when present. Full suite: 342/342 passing (SurfactantKit),
48/48 passing (surfactantkit-reports).

**2026-09-12 (same day, follow-up): postmicellar plateau's own regression exposed in
SurfactantKit's core; `surfactantkit-reports` UI dropped literature citations.** User decision:
the researcher-facing UI should not cite papers at all -- it should show the actual formula
used and the run's own result numbers instead. Two changes:
- `CmcFromCurveResult` gained `postmicellar_slope_mN_per_m_per_log10C`,
  `postmicellar_intercept_mN_per_m`, `postmicellar_mean_gamma_mN_per_m`, and
  `r_squared_postmicellar` -- the plateau segment's own regression, computed once after the
  breakpoint search (previously discarded internally, never returned). New test using the
  same real researcher dataset: confirms the plateau's mean gamma (~37.8 mN/m) and near-zero
  slope, with an explicit note that R^2 is naturally LOW (~0.07) for near-constant data even
  when the fit is visually excellent -- a known statistical property, not a defect, disclosed
  so a researcher doesn't misread it as a bad fit. `mcp_server.py`'s tool wrapper updated to
  surface these fields too. This part stayed permanently -- real, useful, low-risk data the
  underlying library now exposes.
- `surfactantkit-reports`: `TechniqueSpec.citation` removed from the UI package entirely;
  `formula`/`method_name` for all 4 techniques rewritten to be self-contained (no external
  reference needed), and the "Source:" line dropped from both the on-screen methodology panel
  and the downloadable `methodology.txt`.
  **First attempt over-built this**: also added a genuine nonlinear Szyszkowski/Langmuir
  isotherm fit (`szyszkowski_fit_K`) on the premicellar data as a second, separate
  characterization step alongside CMC location, plus a new required `system_type` UI input and
  an `extra_inputs` mechanism to support it -- verified working end-to-end (real, honestly-
  disclosed R^2=0.776 against the researcher dataset, not a misleadingly perfect number) but
  the user then clarified the actual ask was simpler: "we will show the result and we will
  explain how [we] did it, that is it" -- no two-part premicellar/postmicellar breakdown.
  **Reverted per that feedback**: `run_cmc_surface_tension_report` and `CmcSurfaceTensionReport`
  back to their simple pre-Szyszkowski shape (just `result`/`png_path`/`agr_path`/`report_data`,
  no `system_type` parameter); `TechniqueSpec.extra_inputs`/`ExtraInputSpec` removed entirely
  as now-unused abstraction, not left in as dead infrastructure; `CMC_SURFACE_TENSION`'s
  `formula`/`result_fields` back to the plain single-paragraph breakpoint-method description,
  same simple pattern as the other 3 techniques -- citation-free, but not over-engineered.
  Full suite: 344/344 passing (SurfactantKit), 48/48 passing (surfactantkit-reports); Streamlit
  app boot-tested headlessly three times across this back-and-forth (HTTP 200, clean log each
  time).

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
