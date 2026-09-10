"""Clint ideal mixing and Rubingh regular-solution theory for binary
surfactant mixtures.

Ported and generalized from a calculator originally built for a specific
cationic/anionic (amidoamine gemini / SLG) binary system, validated
against real conductometry and tensiometry data. Component 1 and
component 2 are generic here -- there is nothing cationic-specific left
in the math.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .thermodynamics import cmc_to_mole_fraction

R_GAS = 8.314462618  # J / (mol K)


def clint_ideal_cmc(alpha1: float, cmc1: float, cmc2: float) -> float:
    """Clint's ideal mixed CMC for a binary surfactant mixture.

    alpha1: bulk mole fraction of component 1 in the total surfactant.
    cmc1, cmc2: pure-component CMCs (same concentration unit; the
    returned value is in that same unit).
    """
    if not (0.0 < alpha1 < 1.0):
        raise ValueError("alpha1 must be strictly between 0 and 1")
    if cmc1 <= 0 or cmc2 <= 0:
        raise ValueError("pure-component CMCs must be positive")
    return 1.0 / ((alpha1 / cmc1) + ((1.0 - alpha1) / cmc2))


def _rubingh_residual(x1: float, alpha1: float, cmc_mix: float, cmc1: float, cmc2: float) -> float:
    term1 = math.log((alpha1 * cmc_mix) / (x1 * cmc1))
    term2 = math.log(((1.0 - alpha1) * cmc_mix) / ((1.0 - x1) * cmc2))
    return (x1 * x1 * term1) - (((1.0 - x1) ** 2) * term2)


def solve_rubingh_x(
    alpha1: float,
    cmc_mix: float,
    cmc1: float,
    cmc2: float,
    grid_points: int = 4000,
    bisection_iters: int = 100,
) -> float | None:
    """Solve Rubingh's regular-solution equation for x1, the micellar mole
    fraction of component 1.

    Scans a fine grid for a sign change in the residual, then bisects to
    high precision. Returns None if no root is found in (0, 1) -- this
    happens for genuinely ideal or near-ideal mixtures, or for input data
    that isn't self-consistent (e.g. cmc_mix outside the physically
    possible range for the given alpha1).
    """
    eps = 1e-6
    xs = [eps + i * (1.0 - 2.0 * eps) / grid_points for i in range(grid_points + 1)]
    values: list[tuple[float, float]] = []
    for x in xs:
        try:
            fx = _rubingh_residual(x, alpha1, cmc_mix, cmc1, cmc2)
            if math.isfinite(fx):
                values.append((x, fx))
        except (ValueError, ZeroDivisionError):
            continue
    if not values:
        return None

    for (x_lo, f_lo), (x_hi, f_hi) in zip(values, values[1:]):
        if f_lo == 0:
            return x_lo
        if f_lo * f_hi < 0:
            lo, hi = x_lo, x_hi
            flo = f_lo
            for _ in range(bisection_iters):
                mid = (lo + hi) / 2.0
                fmid = _rubingh_residual(mid, alpha1, cmc_mix, cmc1, cmc2)
                if abs(fmid) < 1e-14:
                    return mid
                if flo * fmid <= 0:
                    hi = mid
                else:
                    lo, flo = mid, fmid
            return (lo + hi) / 2.0

    best_x, best_val = min(values, key=lambda pair: abs(pair[1]))
    return best_x if abs(best_val) < 1e-6 else None


def rubingh_beta(x1: float, alpha1: float, cmc_mix: float, cmc1: float) -> float:
    """Rubingh interaction parameter beta, given the solved micellar mole
    fraction x1 (see solve_rubingh_x)."""
    if x1 <= 0.0 or x1 >= 1.0:
        raise ValueError("x1 must be strictly between 0 and 1")
    term = math.log((alpha1 * cmc_mix) / (x1 * cmc1))
    return term / ((1.0 - x1) ** 2)


@dataclass
class RubinghBetaRegressionResult:
    beta_mean: float
    beta_std: float  # population std dev across the used points; 0.0 if only 1 point
    betas: list[float]  # individual pointwise beta values, in input order (excluding skipped points)
    x1_values: list[float]  # corresponding solved micellar mole fractions
    alpha1_used: list[float]  # the alpha1 values that actually solved (skipped points omitted)
    n_points_used: int
    n_points_skipped: int  # points where solve_rubingh_x found no root -- e.g. near-ideal mixtures


def rubingh_beta_regression(
    alpha1_series: list[float],
    cmc_mix_series: list[float],
    cmc1: float,
    cmc2: float,
) -> RubinghBetaRegressionResult:
    """Fit a single representative Rubingh beta across a whole composition
    series -- the real standard practice this project's own
    tests/test_mixed_micelle.py docstring already flagged as missing:
    "many papers fit beta by regression across several alpha
    compositions, while this solver computes the exact pointwise value
    for one composition."

    What "regression across compositions" means in practice, verified
    before implementing (not guessed): compute the pointwise x1 and beta
    at EACH (alpha1, cmc_mix) point independently via the existing
    solve_rubingh_x/rubingh_beta pointwise solver, then report the
    arithmetic MEAN beta across the series as the system's single
    interaction parameter -- this is the real practice used across
    published mixed-surfactant papers (e.g. Rodriguez et al. PMC6554738,
    Muherei & Junin 2009, and others in this project's own
    literature_validation_notes.md), not a continuous nonlinear
    least-squares curve fit against a global cmc_mix(alpha1) model
    (regular solution theory does not define one closed-form
    cmc_mix(alpha1) curve to fit against independent of x1 at each
    point -- x1 itself has to be solved pointwise first, by construction
    of the theory). The returned std dev across the individual betas is
    itself a real, literature-reported diagnostic (large spread signals
    that beta is not actually constant across composition, i.e. RST is
    only an approximate fit for that system) -- returned here rather
    than silently discarded.

    Points where solve_rubingh_x finds no root (near-ideal mixtures, or
    data that isn't self-consistent at that composition) are skipped,
    not treated as an error -- consistent with solve_rubingh_x's own
    None-return contract for those cases.
    """
    if len(alpha1_series) != len(cmc_mix_series):
        raise ValueError("alpha1_series and cmc_mix_series must be the same length")
    if len(alpha1_series) < 2:
        raise ValueError("need at least 2 composition points to regress a single beta")

    betas: list[float] = []
    x1_values: list[float] = []
    alpha1_used: list[float] = []
    n_skipped = 0
    for alpha1, cmc_mix in zip(alpha1_series, cmc_mix_series):
        x1 = solve_rubingh_x(alpha1, cmc_mix, cmc1, cmc2)
        if x1 is None:
            n_skipped += 1
            continue
        beta = rubingh_beta(x1, alpha1, cmc_mix, cmc1)
        betas.append(beta)
        x1_values.append(x1)
        alpha1_used.append(alpha1)

    if len(betas) < 2:
        raise ValueError(
            f"only {len(betas)} of {len(alpha1_series)} points solved -- "
            "need at least 2 to regress a single beta"
        )

    beta_mean = sum(betas) / len(betas)
    variance = sum((b - beta_mean) ** 2 for b in betas) / len(betas)
    beta_std = math.sqrt(variance)

    return RubinghBetaRegressionResult(
        beta_mean=beta_mean,
        beta_std=beta_std,
        betas=betas,
        x1_values=x1_values,
        alpha1_used=alpha1_used,
        n_points_used=len(betas),
        n_points_skipped=n_skipped,
    )


def activity_coefficients(x1: float, beta: float) -> tuple[float, float]:
    """Activity coefficients (f1, f2) from regular solution theory."""
    f1 = math.exp(beta * (1.0 - x1) ** 2)
    f2 = math.exp(beta * x1 ** 2)
    return f1, f2


def excess_free_energy(x1: float, f1: float, f2: float, temperature_k: float) -> float:
    """Excess free energy of micelle formation, kJ/mol."""
    return (R_GAS * temperature_k * (x1 * math.log(f1) + (1.0 - x1) * math.log(f2))) / 1000.0


def asymmetric_margules_activity_coefficients(x1: float, w12: float, w21: float) -> tuple[float, float]:
    """Activity coefficients (f1, f2) from the EOMMM (Equation-Oriented
    Mixed Micellization Model) binary asymmetric Margules formulation --
    Rubingh's RST generalized to allow TWO independent interaction
    parameters (w12 != w21) instead of RST's single symmetric beta.

    w12, w21: dimensionless interaction parameters in RT units (same
    convention as rubingh_beta's return value, and as the literature
    reports them, e.g. "W12 = +4.04 kBT" -- NOT raw J/mol; divide a
    literature W value in J/mol by R*T first). w12 is the energy of
    introducing a molecule of component 1 into a pure micelle of
    component 2 (and vice versa for w21) -- physically the same
    infinite-dilution-limit quantity as a Margules A12/A21 parameter.
    At w12 == w21 == beta this reduces EXACTLY to activity_coefficients
    (RST is the mathematically proven symmetric special case of the
    asymmetric formulation, not a competing formula) -- see
    test_mixed_micelle.py for the reduction check.

    Source: equations (9)-(13) of Schulz & Durand, Comput. Chem. Eng. 87
    (2016) 145-153 (doi:10.1016/j.compchemeng.2015.12.026), restricted
    to a binary (n=2) system -- their general n-component ternary-
    interaction term (Cijk) vanishes identically for n=2 since it
    requires a third component k. This is a real, verified reduction of
    their published general formula, not a guessed simplification, and
    was cross-checked to reproduce their own explicit binary symmetric-
    case reduction (their eq. 14-16, beta*x2^2 / beta*x1^2) exactly.

    NOTE: this function only computes activity coefficients for a
    KNOWN (already-fit) w12/w21 pair, e.g. quoted directly from a
    paper's own reported value. It does NOT solve for w12/w21 (or for
    x1) from raw alpha1/CMC data the way solve_rubingh_x does for RST --
    the real EOMMM procedure fits w12/w21 by globally minimizing total
    free energy across ALL compositions of a dataset SIMULTANEOUSLY
    (Eq. 3.2 of the source paper), not from a single composition point
    the way Rubingh's per-point residual works. That global multi-point
    fit is a materially different, larger feature (a real nonlinear
    equation-oriented optimization) and is NOT implemented here --
    deliberately out of scope, not guessed.
    """
    x2 = 1.0 - x1
    g_exc = x1 * x2 * (x2 * w12 + x1 * w21)
    ln_f1 = 2.0 * x1 * x2 * w21 + (x2 ** 2) * w12 - 2.0 * g_exc
    ln_f2 = 2.0 * x1 * x2 * w12 + (x1 ** 2) * w21 - 2.0 * g_exc
    return math.exp(ln_f1), math.exp(ln_f2)


def _golden_section_minimize(f, lo: float, hi: float, iters: int = 60) -> float:
    """Golden-section search for the minimizer of a function on [lo, hi]
    -- local to this module (mixed_micelle.py has no cross-module
    dependency on adsorption.py's identical helper; SurfactantKit is
    zero-external-dependency by design, and each module stays self-
    contained rather than sharing this small a utility)."""
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c = hi - gr * (hi - lo)
    d = lo + gr * (hi - lo)
    fc, fd = f(c), f(d)
    for _ in range(iters):
        if fc < fd:
            hi, d, fd = d, c, fc
            c = hi - gr * (hi - lo)
            fc = f(c)
        else:
            lo, c, fc = c, d, fd
            d = lo + gr * (hi - lo)
            fd = f(d)
        if hi - lo < 1e-13:
            break
    return (lo + hi) / 2.0


def _eommm_point_residual_sq(x1: float, alpha1: float, cmc_mix: float, cmc1: float, cmc2: float, w12: float, w21: float) -> float:
    f1, f2 = asymmetric_margules_activity_coefficients(x1, w12, w21)
    r1 = math.log(f1) - math.log((alpha1 * cmc_mix) / (x1 * cmc1))
    r2 = math.log(f2) - math.log(((1.0 - alpha1) * cmc_mix) / ((1.0 - x1) * cmc2))
    return r1 * r1 + r2 * r2


def _minimize_with_grid_refine(f, lo: float, hi: float, grid_points: int = 40, refine_iters: int = 40) -> float:
    """Coarse grid scan (finds the global-minimum BASIN) followed by a
    local golden-section refine within it -- needed because
    _eommm_point_residual_sq is NOT always unimodal in x1 (multiple
    near-zero-residual roots can genuinely exist for the same alpha1
    under a shared, fixed W12/W21, directly analogous to the multi-root
    behavior already found and documented for solve_rubingh_x elsewhere
    in this module at large beta). A plain single-bracket golden-section
    search silently converges to whichever basin its initial two probe
    points happen to land in; the grid scan first evaluates across the
    whole domain to find the TRUE global-minimum basin before refining
    within it precisely."""
    step = (hi - lo) / grid_points
    best_x, best_val = lo, f(lo)
    x = lo
    for _ in range(grid_points):
        x += step
        v = f(x)
        if v < best_val:
            best_val, best_x = v, x
    lo2 = max(lo, best_x - 2.0 * step)
    hi2 = min(hi, best_x + 2.0 * step)
    return _golden_section_minimize(f, lo2, hi2, iters=refine_iters)


@dataclass
class EommmGlobalFitResult:
    W12: float
    W21: float
    x1_values: list[float]  # fitted micellar mole fraction at each point, same order as input
    alpha1_used: list[float]
    sse: float  # final minimized total objective (sum of squared log-mass-balance residuals)
    r_squared: float  # diagnostic: fit quality against ln(cmc_mix), see docstring
    n_points: int
    method: str = (
        "first-principles joint least-squares (variable projection: golden-section line search "
        "over W12/W21 in a coordinate-descent outer loop, with each point's x1 profiled out by "
        "an inner grid-scan-plus-golden-section-refine minimization for the current W12/W21, "
        "since the per-point residual is not always unimodal in x1) -- see docstring"
    )


def eommm_global_fit(
    alpha1_series: list[float],
    cmc_mix_series: list[float],
    cmc1: float,
    cmc2: float,
    outer_iters: int = 25,
    w_bound: float = 30.0,
) -> EommmGlobalFitResult:
    """Fit EOMMM's two independent Margules interaction parameters
    (W12, W21) SIMULTANEOUSLY across a whole composition series, plus
    each point's own micellar mole fraction x1 -- what
    asymmetric_margules_activity_coefficients's own docstring flags as
    deliberately out of scope: "the real EOMMM procedure fits w12/w21 by
    globally minimizing total free energy across ALL compositions of a
    dataset SIMULTANEOUSLY... a materially different, larger feature."

    IMPORTANT, disclosed honestly: this is a FIRST-PRINCIPLES
    CONSTRUCTION, not a verified transcription of Schulz & Durand 2016's
    exact published objective function (their Eq. 3.2) or their
    Supplementary Information's detailed procedure ("S.I. Point 2.3").
    Both were checked before building this (2026-09-10): the primary
    paper is paywalled on ScienceDirect (blocked for automated fetch,
    same as the HLD paper earlier this session), and a related open-
    access companion paper by the same group (Serafini et al., a real
    TX100-DTAB EOMMM application reporting W12=+4.04 kBT, W21=-14.02 kBT
    for that system) explicitly defers the exact procedure to its own
    SI, which was not available in the fetched copy. Rather than guess
    at the paper's exact numbering/formulation, this derives a
    mathematically well-justified fit DIRECTLY from the same mass-
    balance equations already verified and used elsewhere in this
    module (the same two equations solve_rubingh_x's residual is built
    from, generalized here to asymmetric_margules_activity_coefficients
    instead of RST's symmetric f1/f2):

        alpha1_i * cmc_mix_i = x1_i * f1(x1_i, W12, W21) * cmc1
        (1-alpha1_i) * cmc_mix_i = (1-x1_i) * f2(x1_i, W12, W21) * cmc2

    Unlike Rubingh (one free parameter, so the two equations at a point
    algebraically eliminate beta and solve x1 exactly per point), EOMMM's
    two independent parameters make a single point's two equations
    OVER-determined once W12/W21 are fixed candidates -- there is
    generally no exact x1 solving both simultaneously except at the true
    global (W12, W21). This is exactly why a genuine joint fit across
    the whole series is needed, not a per-point solve wrapped in an
    outer search.

    Solved via variable projection (a real, standard technique for
    exactly this kind of separable nonlinear least squares, e.g. Golub &
    Pereyra 1973): for ANY fixed (W12, W21), the total sum-of-squared
    mass-balance-residual objective separates into independent per-point
    terms, so each point's own x1 can be "profiled out" by an inner 1D
    minimization; the outer (W12, W21) search then runs coordinate-
    descent (alternating golden-section line search on each parameter)
    over that profiled objective. No numpy/scipy dependency added.

    Real numerical finding made while validating this (2026-09-10),
    directly analogous to the multi-root behavior already found and
    documented for solve_rubingh_x at large beta: for a fixed (W12, W21),
    a single point's mass-balance residual in x1 is NOT always unimodal
    -- two genuinely different x1 values can both give a near-zero
    residual for the same alpha1. A plain single-bracket golden-section
    inner search silently locks onto whichever basin its two initial
    probe points happen to land in, which can be the wrong one even
    while W12/W21 themselves converge correctly (confirmed: the wrong-
    root case still gave r_squared > 0.9999, because the wrong x1 was
    ALSO a genuine near-root, not a search failure). Fixed by making the
    inner per-point solve a coarse grid scan (finds the true global-
    minimum basin) followed by a local golden-section refine within it,
    rather than a single unimodal-assuming bracket search -- see
    _minimize_with_grid_refine.

    A second, related finding: the same non-unimodality means the OUTER
    (W12, W21, all x1_i) joint objective is itself non-convex, so a
    single coordinate-descent run from one starting point can converge
    to a local, not global, optimum even with the per-point fix above.
    Fixed by restarting coordinate descent from 5 starting points spread
    across the search range and keeping the lowest-SSE result (standard
    multi-start global optimization, not an invented technique) --
    confirmed necessary and sufficient via the round-trip test below,
    which failed with a single start and passes with multi-start.
    Corroborating detail found in this project's own ROADMAP.md notes:
    Schulz & Durand's real implementation uses GAMS/BARON -- a GLOBAL
    (not local) nonlinear optimizer -- independently supporting that
    this problem genuinely needs global-search technique, even without
    access to their exact objective function.

    Validated via mathematically-guaranteed round trips only (no
    external published raw numeric table was available given the source
    access gap above): for chosen true (W12, W21, x1_i) values, the
    closed-form construction alpha1_i = x1_i*f1*cmc1/D,
    cmc_mix_i = D where D = x1_i*f1*cmc1+(1-x1_i)*f2*cmc2 satisfies BOTH
    mass-balance equations EXACTLY (verified algebraically and in
    tests/test_mixed_micelle.py) for ANY chosen W12/W21 -- not just
    RST's symmetric special case -- giving a strong, general synthetic
    validation dataset. Also verified: at W12=W21 (the symmetric limit),
    recovers approximately the same beta solve_rubingh_x/rubingh_beta
    find on the identical data, tying this back to the already-
    literature-validated Rubingh machinery.

    Needs at least 3 composition points (2 unknowns W12/W21 plus n
    per-point x1 unknowns vs. 2n mass-balance residuals: needs 2n > n+2,
    i.e. n>=3, to be a genuine fit rather than an exact interpolation);
    5+ recommended. w_bound sets the symmetric search range for W12/W21
    in kBT units (default +/-30, generous relative to real reported
    values like the +4.04/-14.02 kBT example above).
    """
    if len(alpha1_series) != len(cmc_mix_series):
        raise ValueError("alpha1_series and cmc_mix_series must be the same length")
    if len(alpha1_series) < 3:
        raise ValueError("need at least 3 composition points (2n residuals vs n+2 unknowns requires n>=3)")
    if any(not (0.0 < a < 1.0) for a in alpha1_series):
        raise ValueError("all alpha1 values must be strictly between 0 and 1")
    if any(c <= 0 for c in cmc_mix_series):
        raise ValueError("all cmc_mix values must be positive")
    if cmc1 <= 0 or cmc2 <= 0:
        raise ValueError("cmc1 and cmc2 must be positive")

    eps = 1e-6

    def profiled_objective(w12: float, w21: float) -> tuple[float, list[float]]:
        total = 0.0
        x1_fit: list[float] = []
        for alpha1, cmc_mix in zip(alpha1_series, cmc_mix_series):
            x1 = _minimize_with_grid_refine(
                lambda x: _eommm_point_residual_sq(x, alpha1, cmc_mix, cmc1, cmc2, w12, w21),
                eps, 1.0 - eps,
            )
            total += _eommm_point_residual_sq(x1, alpha1, cmc_mix, cmc1, cmc2, w12, w21)
            x1_fit.append(x1)
        return total, x1_fit

    # Multi-start coordinate descent: the joint (W12, W21, x1_1..x1_n)
    # objective is non-convex (a direct consequence of the same per-point
    # multi-root structure _minimize_with_grid_refine works around at
    # the inner level -- see docstring), so a single coordinate-descent
    # run from one starting point can converge to a local, not global,
    # optimum. Restart from several starting points spread across the
    # search range and keep the lowest-SSE result -- standard multi-start
    # global optimization, not an invented technique.
    starts = [(0.0, 0.0), (w_bound / 3.0, -w_bound / 3.0), (-w_bound / 3.0, w_bound / 3.0),
              (w_bound * 2.0 / 3.0, w_bound * 2.0 / 3.0), (-w_bound * 2.0 / 3.0, -w_bound * 2.0 / 3.0)]
    best_w12, best_w21, best_sse, best_x1 = None, None, math.inf, None
    for w12_start, w21_start in starts:
        w12, w21 = w12_start, w21_start
        for _ in range(outer_iters):
            w12 = _golden_section_minimize(lambda w: profiled_objective(w, w21)[0], -w_bound, w_bound)
            w21 = _golden_section_minimize(lambda w: profiled_objective(w12, w)[0], -w_bound, w_bound)
        sse_candidate, x1_candidate = profiled_objective(w12, w21)
        if sse_candidate < best_sse:
            best_w12, best_w21, best_sse, best_x1 = w12, w21, sse_candidate, x1_candidate

    w12, w21, sse, x1_values = best_w12, best_w21, best_sse, best_x1

    ln_cmc_actual = [math.log(c) for c in cmc_mix_series]
    ln_cmc_pred = []
    for alpha1, x1 in zip(alpha1_series, x1_values):
        f1, f2 = asymmetric_margules_activity_coefficients(x1, w12, w21)
        ln_from_1 = math.log(x1 * f1 * cmc1) - math.log(alpha1)
        ln_from_2 = math.log((1.0 - x1) * f2 * cmc2) - math.log(1.0 - alpha1)
        ln_cmc_pred.append((ln_from_1 + ln_from_2) / 2.0)

    n = len(cmc_mix_series)
    mean_ln = sum(ln_cmc_actual) / n
    ss_res = sum((a - p) ** 2 for a, p in zip(ln_cmc_actual, ln_cmc_pred))
    ss_tot = sum((a - mean_ln) ** 2 for a in ln_cmc_actual)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

    return EommmGlobalFitResult(
        W12=w12, W21=w21, x1_values=x1_values, alpha1_used=list(alpha1_series),
        sse=sse, r_squared=r_squared, n_points=n,
    )


def rodenas_x1(alpha1: float, dln_cmc_mix_dalpha1: float) -> float:
    """Rodenas' Gibbs-Duhem-based, "model-independent" micellar mole
    fraction of component 1, from the LOCAL SLOPE of ln(experimental
    mixed CMC) vs. bulk mole fraction alpha1:

    X1_Rod = -(1 - alpha1) * alpha1 * d[ln(cmc_mix)]/d(alpha1) + alpha1

    Unlike solve_rubingh_x (which assumes a specific regular-solution
    activity-coefficient model and solves an implicit equation) or
    motomura_ideal_composition (the zero-interaction ideal-mixing
    limit), this derives composition directly from the Gibbs-Duhem
    relation applied to the REAL cmc_mix(alpha1) curve -- "model
    independent" in the sense of not assuming any particular
    interaction/activity-coefficient form, at the cost of needing the
    local slope of a real multi-point dataset rather than a single
    (alpha1, cmc_mix) pair.

    dln_cmc_mix_dalpha1: d[ln(cmc_mix)]/d(alpha1) at this alpha1 --
    must come from an actual multi-point experimental cmc_mix(alpha1)
    dataset (a fit, or a numerical/finite-difference derivative across
    real neighboring data points), never guessed or assumed analytic.
    This is exactly why this function takes the slope as a direct
    input instead of computing it internally from one data pair the
    way solve_rubingh_x/motomura_ideal_composition do: Rodenas'
    method fundamentally needs the LOCAL SLOPE of the real experimental
    composition curve, which requires several real neighboring data
    points, not just the single composition being evaluated.

    Source: equation (10), Azum, Rub, Alotaibi, Khan & Asiri,
    Biointerface Res. Appl. Chem. 12(6) (2022) 7416-7428 (the same G6
    gemini surfactant paper already used elsewhere in this project's
    test suite for Rubingh/Clint validation), citing Rodenas, Valiente
    & Villafruela, J. Phys. Chem. B 103(21) (1999) 4549-4554, itself
    built on Lange's model.

    Verified two ways: (1) directly against the cited source equation;
    (2) a real, source-independent algebraic identity (checked by
    test): when cmc_mix follows Clint's EXACT ideal-mixing law (zero
    interaction), X1_Rod reduces algebraically to EXACTLY
    motomura_ideal_composition's formula -- both are zero-interaction
    limits of otherwise-different methods and must agree there, and
    they do.
    """
    return -(1.0 - alpha1) * alpha1 * dln_cmc_mix_dalpha1 + alpha1


def _quadratic_lagrange_derivative_at(x0: float, y0: float, x1: float, y1: float, x2: float, y2: float, x_eval: float) -> float:
    """Derivative at x_eval of the unique quadratic interpolating
    (x0,y0), (x1,y1), (x2,y2) -- exact Lagrange-basis differentiation,
    valid for arbitrarily (unevenly) spaced x0/x1/x2. Reduces exactly to
    the textbook central-difference formula (y2-y0)/(2h) when x_eval=x1
    is the midpoint of evenly-spaced points, and to the textbook
    three-point forward/backward-difference formula at an endpoint --
    both verified in tests/test_mixed_micelle.py, not just asserted."""
    d0 = (x0 - x1) * (x0 - x2)
    d1 = (x1 - x0) * (x1 - x2)
    d2 = (x2 - x0) * (x2 - x1)
    l0p = (2.0 * x_eval - x1 - x2) / d0
    l1p = (2.0 * x_eval - x0 - x2) / d1
    l2p = (2.0 * x_eval - x0 - x1) / d2
    return y0 * l0p + y1 * l1p + y2 * l2p


@dataclass
class RodenasSeriesResult:
    alpha1_used: list[float]  # sorted ascending
    dln_cmc_mix_dalpha1: list[float]  # local-slope estimate at each alpha1_used
    x1_rodenas: list[float]  # rodenas_x1 evaluated at each point, chained through automatically
    n_points: int
    method: str = "local quadratic (Lagrange) fit through the 3 nearest neighboring points, differentiated analytically at each composition"


def rodenas_x1_series(alpha1_series: list[float], cmc_mix_series: list[float]) -> RodenasSeriesResult:
    """Compute Rodenas' d[ln(cmc_mix)]/d(alpha1) -- and the resulting x1
    at every composition -- directly from a RAW multi-alpha1 cmc_mix(alpha1)
    series, closing the gap rodenas_x1's own docstring flags: that
    function takes the local slope as an already-known input "because
    Rodenas' method fundamentally needs the LOCAL SLOPE of the real
    experimental composition curve, which requires several real
    neighboring data points" -- this is that numerical-differentiation
    step.

    Real, standard numerical-differentiation method for an unevenly-
    spaced discrete series (the source papers describe getting "the
    slope from the curve" without specifying their exact numerical
    procedure -- verified by checking, not guessed as surfactant-
    specific): at each composition, fit the unique quadratic through it
    and its 2 nearest neighbors along alpha1, then differentiate that
    quadratic analytically at that exact alpha1. For an interior point
    this uses its left and right neighbor; for the first/last point
    (which has only one side) it uses the nearest 2 points on its
    available side. This is mathematically exact (not approximate) when
    the underlying ln(cmc_mix) truly follows a smooth quadratic locally,
    and reduces exactly to the standard central-difference /
    forward-difference / backward-difference formulas for evenly-spaced
    data (see test_mixed_micelle.py) -- it is a real generalization of
    those textbook formulas to uneven spacing, not an invented method.

    alpha1_series and cmc_mix_series must be the same length, >= 3
    points, with all alpha1 values distinct and strictly between 0 and
    1, and all cmc_mix values positive (needed for ln). Points are
    sorted by alpha1 internally; results are returned in that sorted
    order (alpha1_used).
    """
    if len(alpha1_series) != len(cmc_mix_series):
        raise ValueError("alpha1_series and cmc_mix_series must be the same length")
    if len(alpha1_series) < 3:
        raise ValueError("need at least 3 (alpha1, cmc_mix) points to fit a local quadratic derivative")
    if any(not (0.0 < a < 1.0) for a in alpha1_series):
        raise ValueError("all alpha1 values must be strictly between 0 and 1")
    if any(c <= 0 for c in cmc_mix_series):
        raise ValueError("all cmc_mix values must be positive")

    paired = sorted(zip(alpha1_series, cmc_mix_series), key=lambda p: p[0])
    alphas = [p[0] for p in paired]
    if len(set(alphas)) != len(alphas):
        raise ValueError("alpha1_series contains duplicate values -- cannot fit a local derivative")
    ln_cmc = [math.log(c) for _, c in paired]

    n = len(alphas)
    slopes: list[float] = []
    for i in range(n):
        if i == 0:
            i0, i1, i2 = 0, 1, 2
        elif i == n - 1:
            i0, i1, i2 = n - 3, n - 2, n - 1
        else:
            i0, i1, i2 = i - 1, i, i + 1
        slope = _quadratic_lagrange_derivative_at(
            alphas[i0], ln_cmc[i0], alphas[i1], ln_cmc[i1], alphas[i2], ln_cmc[i2], alphas[i]
        )
        slopes.append(slope)

    x1_values = [rodenas_x1(alphas[i], slopes[i]) for i in range(n)]
    return RodenasSeriesResult(alpha1_used=alphas, dln_cmc_mix_dalpha1=slopes, x1_rodenas=x1_values, n_points=n)


def rodenas_activity_coefficients(alpha1: float, x1: float, cmc_mix: float, cmc1: float, cmc2: float) -> tuple[float, float]:
    """Activity coefficients (f1, f2) for Rodenas' model, given the
    already-determined x1 (see rodenas_x1) -- equations (11)-(12) of
    the source cited in rodenas_x1's docstring:

    f1 = alpha1 * cmc_mix / (x1 * cmc1)
    f2 = (1 - alpha1) * cmc_mix / ((1 - x1) * cmc2)

    Structurally the same mass-balance definition already used
    internally by solve_rubingh_x (see _rubingh_residual's term1/
    term2) -- Rodenas and Rubingh differ only in HOW x1 itself is
    determined (a real experimental curve's Gibbs-Duhem slope vs. an
    assumed regular-solution model), not in how activity coefficients
    follow from an already-known x1.
    """
    if not (0.0 < x1 < 1.0):
        raise ValueError("x1 must be strictly between 0 and 1")
    f1 = (alpha1 * cmc_mix) / (x1 * cmc1)
    f2 = ((1.0 - alpha1) * cmc_mix) / ((1.0 - x1) * cmc2)
    return f1, f2


def maeda_free_energy_of_micellization(
    x1_rub: float,
    beta: float,
    cmc1_M: float,
    cmc2_M: float,
    temperature_k: float,
    water_molarity_m: float = 55.5,
) -> float:
    """Maeda's free energy of micellization (kJ/mol) for a binary
    ionic/nonionic mixed micelle -- an alternative expression to
    excess_free_energy, built specifically for ionic-nonionic systems
    and requiring chain-chain (not just head-group) interaction:

    deltaG_M = RT * (B0 + B1*X1_Rub + B2*X1_Rub^2)
    B0 = ln(Xcmc2)
    B1 + B2 = ln(Xcmc1/Xcmc2)
    B2 = -beta

    where Xcmc1, Xcmc2 are the pure-component CMCs on the MOLE-FRACTION
    scale (see cmc_to_mole_fraction -- the same conversion this
    library's own thermodynamics.py already uses for
    gibbs_free_energy_micellization), and x1_rub/beta are Rubingh's own
    solved micellar mole fraction and interaction parameter (from
    solve_rubingh_x/rubingh_beta) -- Maeda's model reuses RST's own
    x1/beta rather than fitting independently.

    cmc1_M, cmc2_M: pure-component CMCs in mol/L (molarity) -- convert
    from mM/other units before calling.

    Source: equation (9), Azum, Rub, Alotaibi, Khan & Asiri,
    Biointerface Res. Appl. Chem. 12(6) (2022) 7416-7428 (same paper
    already used elsewhere in this project's test suite for Clint/
    Rubingh/Rodenas validation), citing Maeda, J. Colloid Interface
    Sci. 172 (1995) 98-105. Numerically cross-checked (see
    test_mixed_micelle.py): computing B0 = ln(Xcmc2) with this
    project's own already-validated TX-114 CMC (0.263 mM, from the same
    source paper's Table 1) reproduces the paper's own reported
    -B0 = 12.25 (G6+TX-114, alpha1 = 0.22) to within 0.01 -- real
    numeric agreement, not just a formula transcription check.
    """
    xcmc1 = cmc_to_mole_fraction(cmc1_M, water_molarity_m)
    xcmc2 = cmc_to_mole_fraction(cmc2_M, water_molarity_m)
    b2 = -beta
    b1 = math.log(xcmc1 / xcmc2) - b2
    b0 = math.log(xcmc2)
    delta_g_over_rt = b0 + b1 * x1_rub + b2 * (x1_rub ** 2)
    return (R_GAS * temperature_k * delta_g_over_rt) / 1000.0


def solve_rosen_monolayer_x(alpha1: float, c_mix_sigma: float, c1_sigma: float, c2_sigma: float, **kwargs) -> float | None:
    """Rosen's extension of Rubingh's regular-solution theory to the
    mixed ADSORBED MONOLAYER at an interface, from surface tension data
    rather than CMC data. Mathematically identical to solve_rubingh_x --
    the difference is entirely in what the inputs mean, and getting that
    right matters: c1_sigma/c2_sigma are the pure-component surfactant
    concentrations required to reach a chosen REFERENCE surface tension
    (or surface pressure) at the interface, and c_mix_sigma is the total
    mixed-surfactant concentration required to reach that SAME reference
    surface tension at bulk mole fraction alpha1 -- these are NOT CMC
    values. Using CMC values here by mistake silently produces a
    micellization-monolayer confusion, not an error."""
    return solve_rubingh_x(alpha1, c_mix_sigma, c1_sigma, c2_sigma, **kwargs)


def rosen_beta_sigma(x1_sigma: float, alpha1: float, c_mix_sigma: float, c1_sigma: float) -> float:
    """Rosen's monolayer interaction parameter beta^sigma, given the
    solved monolayer mole fraction x1_sigma (see solve_rosen_monolayer_x).
    Same formula as rubingh_beta -- see that function's docstring for
    the input-meaning distinction (surface-tension-derived concentrations,
    not CMC values)."""
    return rubingh_beta(x1_sigma, alpha1, c_mix_sigma, c1_sigma)


def motomura_ideal_composition(alpha1: float, cmc1: float, cmc2: float) -> float:
    """Motomura & Aratono's ideal micellar composition for a binary
    surfactant mixture:

    X1_id = (alpha1 * cmc2) / (alpha1 * cmc2 + alpha2 * cmc1)

    alpha1: bulk mole fraction of component 1 in the total surfactant
    (0 to 1, exclusive). cmc1, cmc2: pure-component CMCs (any consistent
    concentration unit -- the result is dimensionless and unit-
    independent, since cmc1/cmc2 only ever appear as a ratio).

    This is Motomura's companion to Clint's ideal-mixing CMC formula
    (clint_ideal_cmc): both assume zero net interaction between the two
    surfactants, but Clint predicts the ideal MIXED CMC while Motomura
    predicts the ideal MICELLAR COMPOSITION (X1) at that CMC. The two
    are algebraically linked -- X1_id = alpha1 * CMC_mix_ideal / cmc1 --
    and this identity is what test_mixed_micelle.py checks, since it
    holds regardless of any external literature value.

    Source: equation (8), Serafini et al., Colloids Surf. A 562 (2019)
    170-181 (arXiv:1806.09721), citing Motomura & Aratono's original
    ideal-mixing composition formula. Like Clint's model, this is a
    ZERO-PARAMETER ideal-mixing prediction -- it does not use or need an
    interaction parameter, unlike Rubingh's regular-solution beta.
    """
    if not (0.0 < alpha1 < 1.0):
        raise ValueError("alpha1 must be strictly between 0 and 1")
    if cmc1 <= 0 or cmc2 <= 0:
        raise ValueError("pure-component CMCs must be positive")
    alpha2 = 1.0 - alpha1
    return (alpha1 * cmc2) / (alpha1 * cmc2 + alpha2 * cmc1)


def corrin_harkins_predict_cmc(cmc1: float, c_counterion1: float, cmc2: float, c_counterion2: float, c_counterion_target: float) -> tuple[float, float]:
    """Predict CMC at a new counterion (salt) concentration using the
    Corrin-Harkins log-linear relation: log10(CMC) = -g*log10(C) + b.
    Fits g and b from two known (CMC, counterion concentration) points,
    then predicts CMC at c_counterion_target. Returns (predicted_cmc, g).

    g and b are SYSTEM-SPECIFIC (surfactant + salt identity), not
    universal constants -- this is why two real data points are required
    as input rather than a lookup table. All CMC and concentration
    values must be in the same, consistent unit (e.g. all mM); the
    prediction is only reliable as an interpolation between (or close
    extrapolation beyond) the two given concentrations.
    """
    if cmc1 <= 0 or cmc2 <= 0 or c_counterion1 <= 0 or c_counterion2 <= 0 or c_counterion_target <= 0:
        raise ValueError("all CMC and concentration values must be positive")
    if c_counterion1 == c_counterion2:
        raise ValueError("c_counterion1 and c_counterion2 must differ to fit a slope")
    log_cmc1, log_cmc2 = math.log10(cmc1), math.log10(cmc2)
    log_c1, log_c2 = math.log10(c_counterion1), math.log10(c_counterion2)
    g = (log_cmc2 - log_cmc1) / (log_c1 - log_c2)
    b = log_cmc1 + g * log_c1
    log_cmc_target = -g * math.log10(c_counterion_target) + b
    return 10.0 ** log_cmc_target, g
