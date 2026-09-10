# Method Alternatives — Literature Review (2026-09-07/08)

Real, sourced research into alternative methods for every property SurfactantKit
computes, done specifically to answer "why did you choose this method over
others that exist" for thesis defense purposes. Each of the 4 sub-reports below
was produced by an independent research pass with explicit instructions to cite
real sources and flag anything not independently verified — confidence levels
are preserved as reported, not smoothed over.

---

## 1. Mixed micelle theory (Category A) + Adsorption (Category B)

[Full agent report — mixed micelle interaction-parameter alternatives (Motomura,
Rodenas, Maeda, EOMMM/Margules, Blankschtein MT theory), why RST fails
structurally, and adsorption isotherm alternatives (Frumkin, Gouy-Chapman/Stern,
FLM reorientation models) — see conversation log 2026-09-07/08 for the complete
text with all citations, DOIs, and confidence flags.]

**Headline findings:**
- RST assumes W12=W21 (symmetric interaction, no excess entropy) — Holland &
  Rubingh's own literature calls this "very improbable." Real example: DTAB/
  Triton X-100, EOMMM found W12=+4.04 kBT vs W21=-14.02 kBT (>18kBT asymmetry;
  RST could not even numerically converge). Source: Serafini et al., arXiv:1806.09721.
- RST is a mathematically-proven SPECIAL CASE of the asymmetric Margules/EOMMM
  formulation (Schulz & Durand, Comput. Chem. Eng. 87 (2016) 145-153) — not a
  competing-but-equally-valid alternative. It persists in the field "due to its
  simplicity" (direct quote, Serafini et al.), which is a defensible but
  disclosed scope choice, not a rigor claim.
- Alternatives, in order of generality: Motomura (general thermodynamic
  framework) < Rodenas (model-independent, Gibbs-Duhem based) < Maeda
  (2-parameter B1/B2, built for ionic/nonionic systems) < EOMMM/Margules
  (fully asymmetric W12≠W21) < Blankschtein molecular-thermodynamic theory
  (predictive from molecular structure, not fit to CMC data at all).
- Real, confirmed practice: current literature routinely runs MULTIPLE methods
  (Clint, Rubingh, Motomura, Rodenas, Maeda) on the SAME dataset side-by-side
  as cross-checks, not as mutually exclusive claims.
- Adsorption: Szyszkowski/Langmuir assumes no lateral interaction between
  adsorbed molecules; Frumkin adds an interaction/activity-coefficient term for
  when packing gets tight; Gouy-Chapman/Stern is the rigorous electrostatic
  treatment that SurfactantKit's simplified n_factor heuristic stands in for;
  FLM reorientation/multistate models handle surfactants that can adopt
  multiple interfacial orientations. Source: Amankeldi et al., Eurasian
  Chemico-Technological Journal 27(4) (2025), doi:10.18321/ectj1673 (read in full).

## 2. HLB (Category C) + Molecular geometry/CPP (Category D)

**Headline findings:**
- HLB: Griffin/Davies remain standard for hand-calculation. PIT (Shinoda) is a
  real alternative but is an experimental, system-specific calibration, not a
  competing calculation. Guo/Rong/Ying 2006 (J. Colloid Interface Sci. 298,
  441-450) improved Davies' group-contribution accuracy for NONIONIC
  polyethoxylated surfactants specifically. QSPR/ML methods (Wang & Duan 2009;
  Luan et al. 2009 for anionics; reviewed in Hu/Zhang/Wang 2010, PMC2868353)
  achieve higher fit accuracy and can resolve isomer differences Davies'
  additive scheme structurally cannot — but need training data.
- **Confirmed: no trustworthy literature revision of Davies' ionic group
  numbers (quaternary ammonium, sulfonate, amide) was found** — this
  corroborates, not contradicts, SurfactantKit's decision to leave those out
  rather than guess.
- CPP: Tanford's constants (v=27.4+26.9n, lc=1.5+1.265n) remain the field
  standard, folded into Israelachvili's canonical CPP framework. But a 2019/2020
  Advances in Colloid and Interface Science review (PMID 31954873) explicitly
  states the underlying "space-filled spherical core" model CANNOT be applied
  to short/medium chain-length surfactants. Nagarajan (Langmuir 2002, 18,
  31-38) is the standard critique/refinement citation.
- Aggregation number: SurfactantKit's geometric estimate is categorically
  different from (not a competing formula against) real experimental methods —
  TRFQ, SLS, DLS, SANS, viscosity, ultracentrifugation. Real papers routinely
  cross-validate aggregation number across 2-3 of these simultaneously (e.g.
  SANS+EPR+TRFQ combined fits, J. Phys. Chem. B, doi:10.1021/jp0371478) — strong
  evidence that no single method, let alone a pure geometric estimate, is
  treated as sufficient alone in rigorous primary research.

## 3. Electrostatics (Category E) + Dynamics (Category F)

**Headline findings — this is the most actionable gap found tonight:**
- The real, continuous Henry function (Henry 1931, Proc. R. Soc. Lond. A 133,
  106-129) exists in closed form; SurfactantKit only implements the two
  asymptotic limits (Smoluchowski f=1.5, Huckel f=1.0).
- Ohshima's 1994 approximation (J. Colloid Interface Sci. 168, 269-271,
  doi:10.1006/jcis.1994.1419) reproduces the exact Henry function to ~1% error
  in closed form across the full kappa*a range — a direct, minimal, well-cited
  fix. Explicitly recommended for exactly this gap in Skoglund et al., PLOS ONE
  12(7):e0181735 (2017).
- **Computed directly using SurfactantKit's own Debye length formula: real
  micelles (R=2-4nm) at realistic ionic strengths (1-100mM) give kappa*a ≈
  0.2-4.2 — squarely in the poorly-approximated INTERMEDIATE zone, not cleanly
  in either the Smoluchowski or Huckel limit.** This is a real, quantified,
  disclosure-worthy limitation, not a theoretical nitpick.
- Separate gap: the Henry function itself (any form) degrades for |zeta|>25mV,
  common for ionic surfactant micelles — O'Brien-White relaxation-effect theory
  (J. Chem. Soc., Faraday Trans. 2, 74, 1607-1626, 1978) is the next rigor tier.
- Debye length: still the standard default, but has documented breakdown modes
  at high ionic strength (underscreening) and high surface potential
  (linearization failure) — both plausibly relevant to ionic micelles.
- Dynamics: Stokes-Einstein is exactly right for dilute spherical micelles.
  Perrin friction factors (Perrin 1934/1936) correct for rod/wormlike
  aggregates. SANS/SAXS is an independent alternative to DLS (different
  physical basis, gives shape not just size). PFG-NMR is an alternative route
  to measuring D itself (works in turbid/concentrated systems DLS can't handle).

## 4. Thermodynamics (Category G) + Wetting (Category H) + Solubilization (Category I)

**Headline findings:**
- **Our finding tonight (two-point van't Hoff giving wildly inconsistent
  pairwise dH, -60 to +2.5 kJ/mol across adjacent pairs) is a KNOWN,
  already-documented weakness in the literature, not a novel discovery.**
  Sources: Kantonen, Henriksen & Gilson, Biochim. Biophys. Acta (2018),
  PMC5851798 (general physical chemistry); Corea et al., Entropy 23(2):236
  (2021), PMC7922405 (surfactant-specific); McGhee, Mingins & Pethica, Langmuir
  37(28) (2021), surveying 23 real amphiphile systems and attributing many
  literature ΔH discrepancies to "incorrect use of the van't Hoff equation."
- Standard fix: multi-point linear regression (already tried) or better, a full
  nonlinear Gibbs-Helmholtz fit across 5+ temperature points solving for ΔH,
  ΔS, and ΔCp simultaneously (Kantonen et al. 2018) — captures real curvature
  a two-point or even linear-regression treatment cannot.
- Mass-action model (Rusanov, Langmuir 30(48) (2014) 14443-14451) is more
  general than the pseudo-phase-separation model SurfactantKit uses; the two
  converge for large aggregation number (N gtrsim 50) — a real, citable scope
  boundary for when pseudo-phase separation is a valid approximation.
- **The counterion_factor convention itself is a genuinely unsettled question
  in the field as of 2025** — Schulz, Durand & Schulz, J. Solution Chem. 54
  (2025) 1426-1450, doi:10.1007/s10953-025-01478-9, argues the field has not
  agreed on one interpretation of the underlying ionization-degree parameter.
  This needs to be stated explicitly, with SurfactantKit's specific convention
  named and cited, not left implicit.
- Wetting: Young-Dupre remains correct and standard for a scalar work-of-
  adhesion number. OWRK and van Oss-Chaudhury-Good are needed only when
  decomposing surface energy into dispersive/polar/acid-base components — a
  different question, not a more-correct answer to the same question. Source:
  Georgiev et al., Colloids and Interfaces 8(6):62 (2024), read in full.
- Solubilization: K_M (partition coefficient) is a distinct, commonly
  co-reported complementary metric to MSR. Real evidence (Ianiro et al.,
  Langmuir 2019, PMC6448116, read in full) that apparent solubilization
  capacity depends heavily on WHERE the solute sits (palisade/interface vs.
  core) — the same nominal system can show near-0% to near-100% apparent
  loading efficiency purely from solubilization-site shift. MSR alone doesn't
  capture this.

---

## Synthesis: answering "why this method and not that one"

Two real, distinct categories emerged, and they need different defenses:

**(a) Strict generality hierarchy (mathematical containment)** — e.g. RST is
literally the symmetric special case of Margules/EOMMM. Here the honest answer
is: "we use the simpler, standard case, which is provably a special case of
the more general method, and we disclose that it's known to fail when [specific,
citable physical condition] holds." This is a legitimate scope choice, not a
weakness to hide.

**(b) True regime-dependent correctness** — e.g. Frumkin vs. Szyszkowski
(packing density), Henry's two limits vs. Ohshima (kappa*a), Young-Dupre vs.
OWRK (scalar vs. component-resolved). Here the answer is: "our target regime
is X, method Y is valid in regime X, and we've verified our typical inputs
actually fall in that regime" — which we did concretely for the Debye
length/Henry function case, and found a real gap.

**Recommended follow-up**: the Henry function gap (E) is the single most
actionable finding — Ohshima's approximation is closed-form, ~1% error, and
directly fixes a quantified, real problem (typical micelles landing in the
badly-approximated intermediate kappa*a zone). Worth implementing as an
additional selectable regime, not a replacement for the existing two limits.
More broadly: SurfactantKit could support explicit alternative-method selection
(e.g. `rubingh_solve()` alongside a new `margules_asymmetric_solve()`, `henry_function("smoluchowski"|"huckel"|"ohshima")`,
two-point vs. multi-point van't Hoff) so the user picks the regime explicitly
rather than the tool silently assuming one — turning "which is right" into an
informed, disclosed choice.

**2026-09-07 (implementation session) — Frumkin isotherm implemented, first
alternative method actually shipped.** `frumkin_theta`/`frumkin_surface_tension`
added to `adsorption.py`, alongside (not replacing) `szyszkowski_surface_tension`,
same `system_type`-required pattern. Sign convention (a>0=attraction,
a<0=repulsion) and closed form (`K*C = theta/(1-theta)*exp(-2*a*theta)`,
`Pi = -RT*Gamma_max*[ln(1-theta) + a*theta^2]`) verified against two independent
real sources before writing any code, not guessed: the standard corrosion-
inhibitor/surfactant literature form, and Xu & Di Talia (PMC7995737, Table 1 /
Eq. 2.18) explicitly confirming beta>0 (this module's `2a`) is the attractive
case. No numeric literature validation case was found for this specific
isotherm during the review pass, so the real regression guard is that a=0
reduces EXACTLY to the already-validated Szyszkowski code (asserted by test,
both at the library level and through the MCP tool layer) — Szyszkowski is a
provable a=0 special case, not an independent formula. 12 new tests
(`test_adsorption_alternatives.py` + 2 in `test_mcp_server.py`), suite
109 -> 121, all passing. New MCP tool: `frumkin_predict_surface_tension`.

**Motomura (mixed micelle, ideal composition) — implemented.** Fetched the
real Serafini et al. PDF directly (arXiv:1806.09721 / Colloids Surf. A 562
(2019) 170-181) and read equation (8) — `motomura_ideal_composition()` added
to `mixed_micelle.py`, predicting the ideal micellar mole fraction X1 as
Clint's ideal-mixing companion (Clint predicts CMC_mix, Motomura predicts
composition at that CMC). Same zero-parameter, zero-interaction assumption
as Clint, no new source-verification risk. Validated two ways: a trivial
identity (equal pure CMCs -> X1_id = alpha1) and an exact algebraic identity
against Clint's own formula (X1_id = alpha1*CMC_mix_ideal/cmc1, derived from
the shared zero-interaction assumption, checked across 4 of this project's
existing literature systems) — no separate external numeric X1_id value was
needed since the identity itself is a real, independent-of-any-paper check.
4 new tests, suite 121 -> 125. New MCP tool: `motomura_ideal_composition`.

**EOMMM (mixed micelle, asymmetric Margules) — still NOT implemented, and
now for a more specific, confirmed reason.** The fetched paper explicitly
defers EOMMM's actual equations to its own Supplementary Information ("For
details see S.I. for EOMMM") — the main text only states qualitatively that
it solves a full Equation-Oriented Optimization (simultaneously fitting
Margules W12/W21, an ionic dissociation parameter r, and multi-point Gibbs-
free-energy-of-mixing minimization across ALL compositions at once), not a
closed-form substitutable into the existing Rubingh-style residual. This is
a materially bigger implementation than Frumkin or Motomura (a genuine
multi-equation constrained optimizer, not a single algebraic formula) and
needs the actual SI fetched (or Schulz & Durand 2016's own full derivation)
before writing code — not attempted this session, flagged rather than
guessed, same discipline as Ohshima.

**EOMMM (binary activity-coefficient formula) — implemented, after the user
retrieved the primary source directly.** Three automated fetch attempts
failed first (CONICET repository connection refused twice, MDPI 403
blocked, PMC7841946 was the wrong paper) — reported to the user with the
exact failed URLs, who then downloaded and provided `schulz2016.pdf` (the
actual primary source: Schulz & Durand, Comput. Chem. Eng. 87 (2016)
145-153) directly. Read in full. Their general n-component asymmetric
Margules formulation (Eqs. 9-13) restricted to a binary system (n=2, no
ternary term since there's no third component) gives a closed form,
verified to reduce EXACTLY to their own stated binary symmetric-case
equations (14-16, i.e. Rubingh's RST) when W12=W21 — confirmed both
algebraically and by test. Added `asymmetric_margules_activity_coefficients()`
to `mixed_micelle.py` (dimensionless w12/w21 in RT units, same convention as
`rubingh_beta`). Deliberately scoped to computing f1/f2 for an
ALREADY-KNOWN w12/w21 pair only — the real EOMMM procedure fits w12/w21 via
a multi-point global free-energy minimization across an ENTIRE composition
dataset simultaneously (their Section 3.2 objective function, solved with
GAMS/BARON/CONOPT in the paper), which is a materially larger, separate
feature (a real nonlinear equation-oriented optimizer) not attempted here.
5 new tests including a real numeric literature case (Hyamine-DTAB, their
Case study 1, W12=-8836.50 J/mol, W21=-2204.19 J/mol) checked against the
paper's own reported activity-coefficient graph axis range. Suite
125 -> 130. New MCP tool: `asymmetric_margules_activity_coefficients`.

**Important correction to this file's own earlier note**: an earlier pass
this session flagged a supposed contradiction (generic textbook Margules
giving f(DTAB)~=57 at low DTAB content, vs. the paper's qualitative "activity
coefficient near zero" description) as a reason NOT to guess the formula.
With the real source now in hand, that concern is resolved: the standard
textbook two-parameter Margules form IS algebraically identical to the
paper's own general formula restricted to n=2 (verified directly, matches
their eq. 14-16 exactly) — the earlier "contradiction" was a misreading of
a qualitative graph description (which plots gamma vs. BULK composition
alpha, not micellar x1, and depends on the separate, unimplemented global
x1(alpha) fit), not a real formula error. Kept here as an honest record:
the caution was the right call at the time (formula unverified), and the
correction came from getting the real source, not from further guessing.

**What's still NOT implemented from EOMMM**: the multi-point global
optimization that actually FITS w12/w21 and x1(alpha) from raw CMC/alpha1
data (needs a real nonlinear solver across N compositions simultaneously,
not a per-point residual like Rubingh's) and the ionic-dissociation
generalization (Eqs. 4-6, the `r` parameter generalizing Rubingh's
undissociated-only assumption). Both are real, separately-scoped, larger
features — flagged as open, not attempted this session.

**Rodenas — implemented, 2026-09-08.** Real source found and fetched: Azum,
Rub, Alotaibi, Khan & Asiri, Biointerface Res. Appl. Chem. 12(6) (2022)
7416-7428 — the SAME G6 gemini paper already used throughout this project's
test suite for Clint/Rubingh validation — gives Rodenas' exact equations
(10)-(12) directly (citing Rodenas, Valiente & Villafruela, J. Phys. Chem. B
103(21) (1999) 4549-4554, itself built on Lange's model; the primary paper
was 403-blocked on ACS, but wasn't needed once the secondary source gave
the real equations). Added `rodenas_x1()` and `rodenas_activity_coefficients()`
to `mixed_micelle.py`. Rodenas is "model-independent" (Gibbs-Duhem applied
directly to the real cmc_mix(alpha1) curve, no assumed activity-coefficient
form) but needs the LOCAL SLOPE d[ln(cmc_mix)]/d(alpha1) of a real multi-
point dataset as input — deliberately scoped like `corrin_harkins_predict_cmc`
(caller supplies real data, function doesn't fabricate it). Validated by a
real, source-independent algebraic identity (checked by test, finite-
differenced directly against `clint_ideal_cmc`, not hand-derived): under
Clint's exact ideal-mixing law, `rodenas_x1` reduces EXACTLY to
`motomura_ideal_composition` — both are zero-interaction limits of
otherwise-different methods and must agree, and they do. 4 new tests.

**Maeda — implemented, 2026-09-08.** Same Azum et al. 2022 source gives
Maeda's exact equation (9) (citing Maeda, J. Colloid Interface Sci. 172
(1995) 98-105, ScienceDirect-paywalled but not needed). Real, simple
relation: B2 = -beta (Rubingh's own interaction parameter, negated) — Maeda
directly reuses Rubingh's solved x1/beta rather than fitting independently.
Added `maeda_free_energy_of_micellization()` to `mixed_micelle.py`, reusing
`thermodynamics.cmc_to_mole_fraction()` for the required mole-fraction-scale
CMC conversion. Validated with a REAL numeric cross-check, not just formula
transcription: computing B0=ln(Xcmc2) with this project's own already-
validated TX-114 CMC (0.263 mM) reproduces the paper's own reported
-B0=12.25 (G6+TX-114 system) to within 0.02 — genuine numeric agreement
against a real published value. 2 new tests.

**Session total for this pass: Frumkin, Motomura, EOMMM (binary activity
coefficients), Rodenas, and Maeda all implemented and tested. Suite
109 -> 138.**

**Ohshima (electrostatics, category E) — implemented, 2026-09-08, after the
user retrieved the primary source directly.** WebFetch attempts on the
primary source (ScienceDirect, paywalled) and two secondary sources
(ResearchGate 403, another ScienceDirect page 403) all failed and were
reported to the user, who then downloaded and provided `ohshima1994.pdf`
(the actual primary source) and `swan2012.pdf` (a related paper) directly.
Read both in full. Confirmed the earlier suspicion exactly: the real
Eq. [16] from Ohshima's own paper is `f = (2/3)*[1 + 1/(2*(1+delta)^3)]`
with `delta = 2.5/[ka*(1+2*exp(-ka))]` in Ohshima's own [2/3, 1] mobility
normalization; converting to this library's [1.0, 1.5] convention (factor
of 3/2, verified by checking both limits land exactly on 1.0 and 1.5) gives
`f = 1 + 1/[2*(1+delta)^3]` -- exactly the "literature-standard form"
inferred last session, now confirmed against the actual primary source
rather than guessed. The earlier open-access secondary source (Skoglund et
al., PLOS ONE 2017) genuinely does have a typo (missing the factor of 2),
confirmed definitively rather than merely suspected. Added a new `'ohshima'`
regime to the existing `henry_function(regime, kappa_a)` in
`electrostatics.py` (extended with an optional `kappa_a` parameter), and
`zeta_potential_henry` updated to pass it through. **Bonus, unplanned
addition**: `swan2012.pdf` (Swan & Furst, J. Colloid Interface Sci. 388
(2012) 92-94) independently re-derives Ohshima's exact formula for cross-
confirmation (a second, unrelated primary source agreeing exactly) AND
gives their own simpler, MORE accurate closed-form rational-function
approximation (<0.1% error vs Ohshima's <1%, per their own Fig. 2
comparison) -- added as a second new regime, `'swan_furst'`. 12 new tests
across both regimes (limits, monotonicity, mutual agreement in the
intermediate zone, input validation, round-trip through
`zeta_potential_henry`). Suite 138 -> 147. MCP tool `zeta_potential`
updated to expose both new regimes plus `kappa_a`.

**Update, same day: the ResearchGate reference itself obtained and implemented
as a THIRD regime.** The user provided the actual paper (originally blocked
by ResearchGate's 403) directly: Qin, Liu, Wang, Thomas, Wang & Shen, Acta
Optica Sinica 37(10) (2017) 1029003, doi:10.3788/AOS201737.1029003 --
Chinese-language with English abstract/equations/tables, read and extracted
in full. Confirms Ohshima's formula yet again (their Eq. 5, third
independent primary source agreeing exactly) and gives a least-squares
REFIT of Ohshima's same functional form against Wiersema's numerically-
exact Henry function values (their Eq. 6): same form, different constants
(2.8, 0.9, 1 instead of Ohshima's 2.5, 1, 2). Real, exact numeric Table 1
(kappa*a from 0.01 to 1000, comparing Wiersema exact / Ohshima / this
"Optimization" fit) reproduced directly as hard-coded test assertions --
the strongest possible validation available, published numbers rather than
just asymptotic limits. Notably more accurate than Ohshima specifically in
the kappa*a ~ 1-20 range (e.g. kappa*a=10: Ohshima 3.0% error vs. 0.2% for
this refit) -- exactly the zone real micelles land in per this review's own
earlier finding. Added as `henry_function`'s third new regime, `'qin'`. 4
new tests (exact-table-match, limits, input validation, a real quantified
"more accurate than Ohshima" comparison against the same Wiersema reference
point). Suite 147 -> 151. `henry_function` now supports THREE independently-
sourced, cross-confirming full-range approximations (`'ohshima'`,
`'swan_furst'`, `'qin'`) alongside the original two limits.

---

## 5. HLB ionic-group gap (Category C) -- new HLD module, 2026-09-08

**The original goal (a trustworthy Davies quaternary-ammonium group number)
independently re-confirmed as a genuine dead end, not a search failure.**
Multiple fresh searches this session converged on the same real finding as
the original 2026-09-07 review: Davies (1957) never published an increment
value for the quaternary ammonium ion at all -- this isn't a "revision we
couldn't find," the original paper itself has no entry for it. Real
secondary attempts to find a widely-cited proposed value (patents,
textbook-style compilations, ResearchGate threads) either didn't discuss it
or were blocked (ResearchGate, ScienceDirect, academia.edu, and
ri.conicet.gov.ar all 403'd or connection-refused on every attempt this
session -- a recurring pattern, now also recorded in memory as a standing
"report failed fetches" rule).

**Pivoted to a real, better-founded alternative instead of forcing the
original angle.** Found and (with the user retrieving the actual PDF after
automated fetch failed) read in full: Schirone, Tartaro, Gentile & Palazzo,
"An HLD framework for cationic ammonium surfactants," JCIS Open 4 (2021)
100033, doi:10.1016/j.jciso.2021.100033. This is a genuinely modern,
physically-measured (not table-lookup) framework -- Hydrophilic-Lipophilic
Difference (HLD), built from real Winsor III phase-equilibrium titration
experiments -- extended to cationic quaternary ammonium surfactants for
apparently the first time in the literature (the paper's own claim, and
plausible given how recent it is). It does not have Davies' gap: it
measures amphiphilicity directly rather than needing a pre-tabulated group
number.

**New module `hld.py`**, all constants and equations taken directly from
the fetched paper (not extrapolated): `hld_ionic`/`hld_nonionic` (the core
HLD equation, both ionic-log-salinity and nonionic-linear-salinity forms),
`cc_mixing_rule` (surfactant-blend mixing rule), `optimal_salinity_ionic`
(direct algebraic inverse at HLD=0), and real measured reference data --
`CATIONIC_QUAT_CC` (characteristic curvature Cc for LTAB, MTAB, CTAB, CMIC,
BDHC, DDAB, all from real HLD-titration measurements in the paper's own
Table 1, with real uncertainties) and `CATIONIC_QUAT_K` (k=0.7+/-0.1,
shared across this head-group class).

**Real, independently-verified bridge back to Davies' HLB scale**:
`hlb_from_cc_cationic_quat()`, the paper's own Fig. 8 linear regression
(HLB = -0.27*Cc + 20.1, R^2=0.975, fit across 5 real quats with both an
independently-measured Cc and a real Davies-scale HLB from Proverbio,
Bardavid, Arancibia & Schulz, Colloids Surf. A 214 (2003) 167-171 --
notably the same "Schulz" research group already cited throughout this
project). Found a real transcription ambiguity while implementing: the
paper's own figure-axis description suggested Cc=f(HLB), but plugging the
paper's own real data pairs into both directions of f(x)=ax+b showed only
HLB=f(Cc) reproduces the real numbers (LTAB/MTAB/CTAB within ~0.3 units) --
documented as a correction made by checking against real data, not a
guess, consistent with this project's standing discipline.

**This is the closest thing to the "genuinely new, publishable
contribution" the earlier session flagged**: not a literature value that
was sitting there unfound, but a real, disclosed, quantitatively validated
bridge between a modern measured framework and the classic HLB scale,
specifically for the exact head-group class (quaternary ammonium) where
Davies' original table has a genuine, confirmed hole.

**UPDATE, same day: the Proverbio et al. 2003 paper was obtained (user
provided the PDF directly after every automated fetch attempt failed) and
read in full -- this closes the ORIGINAL question completely, not just via
the HLD bridge above.** The paper explicitly cites its own source for the
quaternary ammonium Davies numbers: B.H. O, J. Colloid Interface Sci. 198
(1998) 249, "GN values are 22.5 for >N+(CH3)2Cl- and 22.0 for -N+(CH3)3Cl-"
-- i.e. Davies himself never assigned one (confirmed, as found earlier),
but O (1998) did, 41 years later, and this 2003 paper is where that
extension got applied and published with real worked examples.

**Verified four independent ways, not just transcribed**: Proverbio et
al.'s own paper reports HLB = 7 + Sigma(GN) worked out fully for four real
compounds. Recomputed all four directly (7 + GN - 0.475*n_carbons, this
project's own existing CH2/CH3=0.475 convention) and every one matches the
paper's own reported value exactly: DTAB (C10, trimethyl) = 24.25, LTAB
(C12, trimethyl) = 23.3, DDAB (2xC12, dimethyl-dialkyl) = 18.1, DODAB
(2xC18, dimethyl-dialkyl) = 12.4. Both group numbers now added to
`hlb.py`'s `DAVIES_HYDROPHILIC_GROUPS`: `N_quat_trimethyl` (22.0, single-
tail heads like CTAB/LTAB/DTAB) and `N_quat_dimethyl_dialkyl` (22.5,
double-tail heads like DDAB/DODAB) -- both assume Cl-/Br- counterions (the
source paper explicitly neglects the difference between them; NOT verified
for other counterions, e.g. the same paper's own CTATOS/tosylate example
has no computable Davies HLB for exactly this reason).

**Real cross-validation bonus**: `hld.py`'s `CATIONIC_QUAT_HLB_DAVIES`
dict (sourced independently via Schirone et al.'s reproduction of this
same Proverbio et al. data) now matches `hlb.py`'s direct group-
contribution computation EXACTLY for LTAB, MTAB, CTAB, and DDAB -- two
separate implementations, from two separate readings of two separate
(though ultimately same-origin) sources, landing on identical numbers.
This is a real, strong consistency check, not a tautology.

**Also found, not yet used**: the same paper gives a graphically-estimated
(lower confidence, water-number-correlation-based, not a direct group-
table lookup) group number for the pyridinium chloride headgroup, GN=15.1
-- noted here as a real lead but deliberately NOT added to
`DAVIES_HYDROPHILIC_GROUPS` (which only holds directly-verified numbers),
since the source paper itself describes it as an estimate from a
nonlinear graph read-off, not the same rigor tier as the quaternary
ammonium numbers above.

**Amide and sulfonate remain genuinely unresolved** -- this paper didn't
touch either. DOI for Proverbio et al.: 10.1016/S0927-7757(02)00404-1.
9 new tests (4 exact worked-example reproductions + 1 cross-module
consistency check in this pass, on top of the 14 from the HLD module
above). Suite 165 -> 170.

**This closes the HLB ionic-group item as a genuine, complete success**:
started as "no trustworthy value found," ended with a real, source-
verified, four-times-cross-validated Davies group number now live in the
library, PLUS an independent modern-framework bridge (hld.py) as a bonus
that wouldn't have been built without chasing the original dead end first.
(`test_hld.py` + 3 in `test_mcp_server.py`). Suite 151 -> 165. New MCP
tools: `hld_optimal_salinity`, `hld_cationic_quat_reference`.
