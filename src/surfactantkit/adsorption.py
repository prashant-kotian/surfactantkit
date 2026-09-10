"""Gibbs adsorption isotherm: surface excess concentration and minimum
area per molecule from a surface-tension-vs-log(concentration) slope."""

from __future__ import annotations

import math
from dataclasses import dataclass

AVOGADRO = 6.02214076e23

# Gibbs prefactor n by system type. Real, SurfBench-benchmark-confirmed
# regression 2026-09-04: this used to be a raw n_factor float with a
# docstring-only "do not guess" warning. Unaugmented Claude Opus 4.8 got it
# wrong 5/6 times (silently defaulting to n=2 for "a bare ionic surfactant"),
# and -- more importantly -- so did tool-augmented Claude (0/6 under the
# surfmcp condition), because a float parameter with no default still lets
# a caller supply ANY plausible-looking number with no way for the function
# to detect it was guessed. Forced to an explicit categorical string
# instead (same pattern as electrostatics.py's henry_function regime),
# raising on anything unrecognized, so at minimum a caller cannot pass a
# bare, unlabelled number through by accident -- see gibbs_gamma_max's own
# docstring for why this alone is not a complete fix.
N_FACTOR_BY_SYSTEM_TYPE = {
    "nonionic": 1.0,
    "ionic_excess_electrolyte": 1.0,
    "ionic_no_added_salt": 2.0,
}


def gibbs_gamma_max(slope_mn_per_ln_c: float, system_type: str, temperature_k: float, r_gas: float = 8.314462618) -> float:
    """Maximum surface excess concentration (mol/m^2) from the pre-CMC
    slope of surface tension (mN/m) vs ln(concentration).

    slope_mn_per_ln_c: dGamma/d(ln C), in mN/m (negative for a surfactant
    lowering surface tension with increasing concentration).
    system_type: one of 'nonionic', 'ionic_excess_electrolyte', or
    'ionic_no_added_salt' -- determines the Gibbs prefactor n (1 for
    nonionic or ionic-with-excess-electrolyte, 2 for a 1:1 ionic
    surfactant with no added salt). No default and no bare numeric
    n_factor accepted: whether excess inert electrolyte is present is
    real, system-specific information that must come from the question/
    experiment, not be assumed. IMPORTANT for callers (including LLMs
    using this via a tool): if the surfactant's ionic character or its
    electrolyte condition is not stated, do not guess a system_type --
    report that Gamma_max cannot be determined instead of calling this
    with an assumed value.
    """
    if system_type not in N_FACTOR_BY_SYSTEM_TYPE:
        raise ValueError(
            f"system_type must be one of {sorted(N_FACTOR_BY_SYSTEM_TYPE)}, got {system_type!r} -- "
            "do not guess; if the surfactant's ionic character or electrolyte condition isn't "
            "stated, Gamma_max cannot be computed"
        )
    n_factor = N_FACTOR_BY_SYSTEM_TYPE[system_type]
    slope_n_per_m = slope_mn_per_ln_c / 1000.0
    return (-1.0 * slope_n_per_m) / (n_factor * r_gas * temperature_k)


def gibbs_a_min(gamma_max_mol_per_m2: float) -> float:
    """Minimum area per molecule (nm^2) from Gamma_max (mol/m^2)."""
    if gamma_max_mol_per_m2 <= 0:
        raise ValueError("Gamma_max must be positive to compute A_min")
    return 1.0e18 / (AVOGADRO * gamma_max_mol_per_m2)


def szyszkowski_surface_tension(concentration: float, gamma0_mN_m: float, gamma_max_mol_per_m2: float, K: float, system_type: str, temperature_K: float = 298.15, r_gas: float = 8.314462618) -> float:
    """Predict surface tension (mN/m) at a given concentration via the
    Szyszkowski/Langmuir equation:
    gamma(C) = gamma0 - n_factor*R*T*Gamma_max*ln(1 + K*C).

    concentration: surfactant concentration, in whatever unit K's inverse
    is defined in (e.g. mM if K is in 1/mM) -- K*concentration must be
    dimensionless, so K and concentration must use consistent units.
    gamma_max_mol_per_m2: saturation surface excess (mol/m^2, same
    quantity as gibbs_gamma_max's output).
    K: Szyszkowski/Langmuir adsorption constant, fit from real
    surface-tension-vs-concentration data (e.g. via szyszkowski_fit_K) --
    not a universal constant, do not guess a value.
    system_type: same Gibbs prefactor selector as gibbs_gamma_max
    ('nonionic', 'ionic_excess_electrolyte', or 'ionic_no_added_salt') --
    required, no default. Real regression found 2026-09-04: this used to
    default silently to n_factor=1.0, which is simply wrong for any 1:1
    ionic surfactant with no added salt (n=2) and was never flagged to a
    caller who didn't override it. See gibbs_gamma_max's docstring for
    the full rationale; if the surfactant's ionic character or
    electrolyte condition isn't stated, do not guess -- report that
    surface tension cannot be predicted instead.

    At C=0 this returns exactly gamma0. As C increases toward and beyond
    the CMC, gamma decreases monotonically toward a plateau -- this
    equation is only meant to be applied below the CMC, where surfactant
    exists as free monomer at the interface.
    """
    if concentration < 0:
        raise ValueError("concentration must be non-negative")
    if gamma_max_mol_per_m2 <= 0:
        raise ValueError("gamma_max_mol_per_m2 must be positive")
    if K < 0:
        raise ValueError("K must be non-negative")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if system_type not in N_FACTOR_BY_SYSTEM_TYPE:
        raise ValueError(
            f"system_type must be one of {sorted(N_FACTOR_BY_SYSTEM_TYPE)}, got {system_type!r} -- "
            "do not guess; if the surfactant's ionic character or electrolyte condition isn't "
            "stated, surface tension cannot be predicted"
        )
    n_factor = N_FACTOR_BY_SYSTEM_TYPE[system_type]
    term_N_per_m = n_factor * r_gas * temperature_K * gamma_max_mol_per_m2 * math.log(1.0 + K * concentration)
    term_mN_per_m = term_N_per_m * 1000.0
    return gamma0_mN_m - term_mN_per_m


@dataclass
class SzyszkowskiFitKResult:
    K: float
    r_squared: float
    n_points: int
    fitted_surface_tensions_mN_per_m: list[float]  # model prediction at each input concentration, same order


def _golden_section_minimize(f, lo: float, hi: float, iters: int = 200) -> float:
    """Golden-section search for the minimizer of a UNIMODAL function f
    on [lo, hi]. Local, pure-Python (no scipy) -- SurfactantKit is
    zero-external-dependency by design."""
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
        if hi - lo < 1e-14:
            break
    return (lo + hi) / 2.0


def szyszkowski_fit_K(
    concentrations: list[float],
    surface_tensions_mN_per_m: list[float],
    gamma0_mN_m: float,
    gamma_max_mol_per_m2: float,
    system_type: str,
    temperature_K: float = 298.15,
    r_gas: float = 8.314462618,
) -> SzyszkowskiFitKResult:
    """Fit the Szyszkowski/Langmuir adsorption constant K from a RAW
    PRE-CMC surface-tension-vs-concentration curve -- closes the gap
    szyszkowski_surface_tension's own docstring references by name
    ("K... fit from real data, e.g. via szyszkowski_fit_K") but which
    did not exist until now (confirmed via direct grep before building
    this, not assumed).

    Real, standard two-step surface-tension-isotherm workflow (not
    invented): Gamma_max is determined independently from the SAME
    curve's pre-CMC ln(concentration) slope (gibbs_gamma_max, via
    cmc_from_surface_tension_curve's premicellar_slope output or a
    direct log-linear fit) and passed in here as an already-known,
    fixed value -- this function then does a genuine nonlinear
    least-squares fit for the ONE remaining free parameter, K, that
    best reproduces the observed gamma(C) curve via the Szyszkowski
    equation itself (szyszkowski_surface_tension). This is a
    DIFFERENT curve regime/purpose than cmc_from_surface_tension_curve,
    which locates the CMC break point -- this one fits the pre-CMC
    region's actual monomer-adsorption curvature.

    concentrations and surface_tensions_mN_per_m must be pre-CMC data
    only (points at or above the CMC violate the Szyszkowski model's own
    assumption and will bias the fit -- do not include them). Solved via
    golden-section search over ln(K) (the least-squares objective is
    unimodal in K here, since the model's predicted surface tension is
    monotonic in K for any fixed positive concentration) rather than a
    general-purpose nonlinear solver, to avoid adding a dependency.
    Returns the fitted K, the fit's r_squared, and the model-predicted
    surface tension at each input concentration for direct visual/
    residual inspection.

    gamma0_mN_m, gamma_max_mol_per_m2, system_type, temperature_K, and
    r_gas all have the same meaning as in szyszkowski_surface_tension.
    """
    if len(concentrations) != len(surface_tensions_mN_per_m):
        raise ValueError("concentrations and surface_tensions_mN_per_m must be the same length")
    if len(concentrations) < 3:
        raise ValueError("need at least 3 pre-CMC (concentration, surface tension) points to fit K meaningfully")
    if any(c <= 0 for c in concentrations):
        raise ValueError("concentrations must be positive (K*0 contributes no information to the fit)")
    if gamma_max_mol_per_m2 <= 0:
        raise ValueError("gamma_max_mol_per_m2 must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if system_type not in N_FACTOR_BY_SYSTEM_TYPE:
        raise ValueError(
            f"system_type must be one of {sorted(N_FACTOR_BY_SYSTEM_TYPE)}, got {system_type!r} -- "
            "do not guess; if the surfactant's ionic character or electrolyte condition isn't "
            "stated, K cannot be fit"
        )

    def predict(k: float) -> list[float]:
        return [szyszkowski_surface_tension(c, gamma0_mN_m, gamma_max_mol_per_m2, k, system_type, temperature_K, r_gas)
                for c in concentrations]

    def sse(ln_k: float) -> float:
        k = math.exp(ln_k)
        pred = predict(k)
        return sum((p - g) ** 2 for p, g in zip(pred, surface_tensions_mN_per_m))

    ln_k_best = _golden_section_minimize(sse, math.log(1e-8), math.log(1e8))
    k_fit = math.exp(ln_k_best)

    fitted = predict(k_fit)
    n = len(surface_tensions_mN_per_m)
    mean_gamma = sum(surface_tensions_mN_per_m) / n
    ss_res = sum((g - p) ** 2 for g, p in zip(surface_tensions_mN_per_m, fitted))
    ss_tot = sum((g - mean_gamma) ** 2 for g in surface_tensions_mN_per_m)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

    return SzyszkowskiFitKResult(K=k_fit, r_squared=r_squared, n_points=n, fitted_surface_tensions_mN_per_m=fitted)


def frumkin_theta(concentration: float, K: float, a: float, grid_points: int = 2000, bisection_iters: int = 100) -> float:
    """Solve the Frumkin adsorption isotherm for surface coverage theta
    at a given bulk concentration:

    K*concentration = [theta / (1 - theta)] * exp(-2*a*theta)

    K: Langmuir-type adsorption equilibrium constant (same role and same
    unit convention as szyszkowski_surface_tension's K), fit from real
    data -- not a universal constant, do not guess a value.
    a: dimensionless Frumkin lateral-interaction parameter. a > 0 means
    NET ATTRACTIVE interaction between adsorbed surfactant molecules
    (favors closer packing, steeper isotherm); a < 0 means net
    repulsive interaction; a == 0 recovers the plain Langmuir isotherm
    exactly (theta = K*C / (1 + K*C), the same coverage-concentration
    relation implicit inside szyszkowski_surface_tension). Sign
    convention and closed form confirmed against two independent real
    sources: the standard K*C = theta/(1-theta)*exp(-2*a*theta) form
    (Frumkin 1925, widely cited in corrosion-inhibitor and surfactant
    adsorption literature), and Xu & Di Talia's review (PMC7995737,
    Table 1 / Eq. 2.18), which states explicitly that beta > 0 (this
    module's 2*a = beta*Gamma_max/(kB*T)) reduces the free energy to
    assemble a monolayer of mutually ATTRACTIVE surfactants.

    For |a| < 2 the isotherm is monotonic in theta and has a single
    root in (0, 1). For a >= 2 the Frumkin isotherm is genuinely
    multi-valued -- a real physical 2D-condensation/phase-transition
    effect (a first-order-like jump in coverage), not a solver bug.
    This function scans for and returns the smallest-theta root in that
    regime (same grid-scan-then-bisect pattern as this module's sibling
    solve_rubingh_x); a caller landing in a >= 2 should treat the result
    as "near/at a condensation transition", not a single trustworthy
    coverage value.
    """
    if concentration < 0:
        raise ValueError("concentration must be non-negative")
    if K < 0:
        raise ValueError("K must be non-negative")
    if concentration == 0.0 or K == 0.0:
        return 0.0

    def residual(theta: float) -> float:
        return (theta / (1.0 - theta)) * math.exp(-2.0 * a * theta) - K * concentration

    eps = 1e-9
    thetas = [eps + i * (1.0 - 2.0 * eps) / grid_points for i in range(grid_points + 1)]
    values = [(t, residual(t)) for t in thetas]

    for (t_lo, f_lo), (t_hi, f_hi) in zip(values, values[1:]):
        if f_lo == 0:
            return t_lo
        if f_lo * f_hi < 0:
            lo, hi, flo = t_lo, t_hi, f_lo
            for _ in range(bisection_iters):
                mid = (lo + hi) / 2.0
                fmid = residual(mid)
                if abs(fmid) < 1e-14:
                    return mid
                if flo * fmid <= 0:
                    hi = mid
                else:
                    lo, flo = mid, fmid
            return (lo + hi) / 2.0

    raise ValueError(
        "no root found for theta in (0, 1) -- check concentration, K, and a for physical consistency"
    )


def frumkin_surface_tension(
    concentration: float,
    gamma0_mN_m: float,
    gamma_max_mol_per_m2: float,
    K: float,
    a: float,
    system_type: str,
    temperature_K: float = 298.15,
    r_gas: float = 8.314462618,
) -> float:
    """Predict surface tension (mN/m) via the Frumkin isotherm, which
    generalizes szyszkowski_surface_tension with a lateral-interaction
    term a (see frumkin_theta's docstring for its sign convention and
    sourcing):

    gamma(C) = gamma0 + n_factor*R*T*Gamma_max*[ln(1 - theta) + a*theta^2]

    where theta solves frumkin_theta(concentration, K, a). At a = 0
    this reduces EXACTLY to szyszkowski_surface_tension -- Szyszkowski/
    Langmuir is the a=0 special case of Frumkin, not two independent
    formulas (verified by test_adsorption_alternatives.py).

    gamma_max_mol_per_m2, system_type, temperature_K, r_gas: same
    meaning as szyszkowski_surface_tension. system_type is required, no
    default, no guessing -- see gibbs_gamma_max's docstring for the
    full rationale.
    """
    if concentration < 0:
        raise ValueError("concentration must be non-negative")
    if gamma_max_mol_per_m2 <= 0:
        raise ValueError("gamma_max_mol_per_m2 must be positive")
    if K < 0:
        raise ValueError("K must be non-negative")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if system_type not in N_FACTOR_BY_SYSTEM_TYPE:
        raise ValueError(
            f"system_type must be one of {sorted(N_FACTOR_BY_SYSTEM_TYPE)}, got {system_type!r} -- "
            "do not guess; if the surfactant's ionic character or electrolyte condition isn't "
            "stated, surface tension cannot be predicted"
        )
    n_factor = N_FACTOR_BY_SYSTEM_TYPE[system_type]
    theta = frumkin_theta(concentration, K, a)
    term_N_per_m = n_factor * r_gas * temperature_K * gamma_max_mol_per_m2 * (math.log(1.0 - theta) + a * theta * theta)
    term_mN_per_m = term_N_per_m * 1000.0
    return gamma0_mN_m + term_mN_per_m
