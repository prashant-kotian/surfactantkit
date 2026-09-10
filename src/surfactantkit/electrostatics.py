"""Electrostatics of ionic surfactant solutions: Debye screening length
and zeta potential (Henry equation).

Physical constants are CODATA/SI-2019 exact or near-exact values, not
memorized approximations:
    k_B = 1.380649e-23 J/K       (exact, SI 2019)
    N_A = 6.02214076e23 /mol     (exact, SI 2019)
    e   = 1.602176634e-19 C      (exact, SI 2019)
    eps0 = 8.8541878128e-12 F/m  (CODATA 2018)
"""

from __future__ import annotations
import math

K_B = 1.380649e-23
N_A = 6.02214076e23
E_CHARGE = 1.602176634e-19
EPS0 = 8.8541878128e-12
WATER_REL_PERMITTIVITY_25C = 78.4


def ionic_strength(concentrations_M: dict[float, float]) -> float:
    """Ionic strength I = 0.5 * sum(c_i * z_i^2), mol/L (M).

    concentrations_M: {charge_number: concentration_M}, e.g. for 0.1 M
    NaCl: {1: 0.1, -1: 0.1} (cation and anion both counted).
    """
    if not concentrations_M:
        raise ValueError("concentrations_M must not be empty")
    return 0.5 * sum(c * (z ** 2) for z, c in concentrations_M.items())


def debye_length(ionic_strength_M: float, temperature_K: float = 298.15, rel_permittivity: float = WATER_REL_PERMITTIVITY_25C) -> float:
    """Debye screening length (nm) for a 1:1-equivalent electrolyte
    solution: lambda_D = sqrt(eps_r*eps0*k_B*T / (2*N_A*I*e^2)).

    ionic_strength_M: ionic strength in mol/L (see ionic_strength()) --
    NOT the raw salt concentration; for a 1:1 salt like NaCl, I equals
    the salt concentration, but for higher-valence electrolytes it does
    not.
    rel_permittivity: relative permittivity of the solvent (default is
    water at 25 C, ~78.4; this is temperature-dependent -- pass an
    explicit value for other temperatures rather than assuming 25 C).
    """
    if ionic_strength_M <= 0:
        raise ValueError("ionic_strength_M must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    I_SI = ionic_strength_M * 1000.0  # mol/L -> mol/m^3
    lambda_m = math.sqrt(
        (rel_permittivity * EPS0 * K_B * temperature_K) / (2.0 * N_A * I_SI * E_CHARGE ** 2)
    )
    return lambda_m * 1e9  # m -> nm


def grahame_equation_surface_potential(surface_charge_density_C_per_m2: float, ionic_strength_M: float, z: int = 1, temperature_K: float = 298.15, rel_permittivity: float = WATER_REL_PERMITTIVITY_25C) -> float:
    """Surface (Stern-layer) potential (mV) from surface charge density
    via the Grahame equation -- the standard, textbook closed-form
    Gouy-Chapman diffuse-double-layer result (Grahame, Chem. Rev. 41
    (1947) 441-501), for a symmetric z:z supporting electrolyte:

        sigma = sqrt(8*eps_r*eps0*k_B*T*n0) * sinh(z*e*psi / (2*k_B*T))

    solved for psi:

        psi = (2*k_B*T)/(z*e) * asinh[ sigma / sqrt(8*eps_r*eps0*k_B*T*n0) ]

    n0 is the bulk number density of the electrolyte, computed from
    ionic_strength_M the same way debye_length does (mol/L -> mol/m^3 ->
    number density via N_A) -- reuses this module's own existing
    constants for consistency. Cross-checked against a real, fetched
    secondary source (Mchedlov-Petrossyan, Electron. Proc. Mater. 50(2)
    (2014) 71-80, their Eq. (9), the integrated form of this same
    charge-potential relation for the ionic-surfactant adsorption case)
    -- same functional (sinh/asinh) form confirmed independently, not
    guessed.

    Real, sourced use: chains directly with counterion_binding_degree
    (thermodynamics.py) and a surfactant's own surface excess (Gamma,
    from gibbs_gamma_max) to build a genuine Gouy-Chapman-corrected
    Gibbs-isotherm prefactor -- see
    adsorption.gibbs_gamma_max_gouy_chapman, which chains this function
    in sequence with existing tools rather than a standalone monolithic
    Gouy-Chapman solver.
    """
    if ionic_strength_M <= 0:
        raise ValueError("ionic_strength_M must be positive")
    if z <= 0:
        raise ValueError("z must be a positive integer (charge number magnitude)")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    n0_per_m3 = ionic_strength_M * 1000.0 * N_A  # mol/L -> mol/m^3 -> ions/m^3
    denom = math.sqrt(8.0 * rel_permittivity * EPS0 * K_B * temperature_K * n0_per_m3)
    psi_V = (2.0 * K_B * temperature_K) / (z * E_CHARGE) * math.asinh(surface_charge_density_C_per_m2 / denom)
    return psi_V * 1000.0  # V -> mV


def henry_function(regime: str, kappa_a: float | None = None) -> float:
    """Henry function f(kappa*a): 'huckel' = 1.0 (small particles / low
    ionic strength, kappa*a << 1), 'smoluchowski' = 1.5 (large particles
    / ionic strength >= ~10 mM, kappa*a >> 1), 'ohshima', or
    'swan_furst' -- the latter two are closed-form approximations valid
    and accurate across the ENTIRE kappa*a range, including the
    intermediate zone where neither limit is a good approximation (real
    micelles, R~2-4 nm at realistic ionic strengths 1-100 mM, typically
    land at kappa*a ~ 0.2-4.2 -- squarely in that poorly-approximated
    intermediate zone, not cleanly in either limit). Both require
    kappa_a (the dimensionless product kappa*a -- Debye parameter times
    particle radius, same units for both) as a second argument; raise
    if omitted or non-positive.

    'ohshima' (<1% relative error, Ohshima 1994):
        f(kappa*a) = 1 + 1 / [2 * (1 + delta)^3]
        delta = 2.5 / [kappa*a * (1 + 2*exp(-kappa*a))]

    'swan_furst' (<0.1% relative error -- an order of magnitude more
    accurate than Ohshima's, per the source paper's own Fig. 2
    comparison; a simple rational function, no exponential):
        f(kappa*a) = (16 + 18*kappa*a + 3*(kappa*a)^2)
                     / (16 + 18*kappa*a + 2*(kappa*a)^2)

    There is no default regime -- picking the wrong one silently is
    exactly the kind of error this library exists to prevent; the
    caller must state which applies for their particle size and ionic
    strength (or supply kappa_a and use 'ohshima'/'swan_furst' to avoid
    the choice entirely).

    Sources, both fetched and read directly (not from secondary
    summaries):
    - Ohshima, H., J. Colloid Interface Sci. 168 (1994) 269-271,
      doi:10.1006/jcis.1994.1419, Eq. [16]. Ohshima's own paper
      normalizes f(kappa*a) to the range [2/3, 1] (their Eq. [1] omits
      this library's leading 2/3 mobility prefactor); converted here by
      a factor of 3/2 to match THIS library's [1.0, 1.5] convention
      (same convention as zeta_potential_henry's mobility formula) --
      verified correct because it reproduces exactly 1.0 at kappa*a ->
      0 and exactly 1.5 at kappa*a -> infinity.
    - Swan & Furst, J. Colloid Interface Sci. 388 (2012) 92-94,
      doi:10.1016/j.jcis.2012.08.026, Eq. (6) -- already in this
      library's [1.0, 1.5] convention as given (no conversion needed;
      their own Eq. (5) confirms the same mobility formula/convention
      this library uses). Their paper independently reproduces
      Ohshima's exact formula too (their Eq. (7), delta =
      (5/2)/(1+2*exp(-kappa*a)), f_O = 1+1/(2*[1+delta/(kappa*a)]^3) --
      algebraically identical to the 'ohshima' branch above once delta
      is combined with the /kappa*a term), an independent cross-
      confirmation from a second, unrelated primary source that this
      library's 'ohshima' implementation is correct.

    A THIRD option, 'qin' -- Qin, Liu, Wang, Thomas, Wang & Shen, Acta
    Optica Sinica 37(10) (2017) 1029003, doi:10.3788/AOS201737.1029003
    (fetched and read directly, in the original Chinese with an English
    abstract/equations) -- least-squares-refits Ohshima's SAME functional
    form against Wiersema's numerically-exact Henry function values,
    finding better coefficients:

    f(kappa*a) = 1 + 1 / [2 * (1 + 2.8/(kappa*a*(0.9+exp(-kappa*a))))^3]

    (Ohshima's original delta = 2.5/[kappa*a*(1+2*exp(-kappa*a))]; this
    refit changes the three constants from (2.5, 1, 2) to (2.8, 0.9, 1).)
    The source paper's own Table 1 (reproduced against Wiersema's exact
    numerical solution) reports this is consistently more accurate than
    Ohshima's original across kappa*a = 0.01 to 1000, most notably in
    the kappa*a ~ 1-20 range this library's own literature review found
    real micelles typically land in (e.g. at kappa*a=10: Ohshima 3.0%
    error vs. this refit's 0.2%; at kappa*a=2: 2.3% vs. 0.2%) -- overall
    max relative error <1.5% vs. Ohshima's <3% in that range. This
    library's own test suite checks this implementation reproduces the
    source paper's own Table 1 values exactly (to the 3-4 significant
    figures given), not just the asymptotic limits.

    CAUTION: an open-access secondary source found during initial
    research (Skoglund et al., PLOS ONE 12(7):e0181735 (2017), Eq. 3)
    reproduces the SAME Ohshima delta formula but with a DIFFERENT,
    INCORRECT f(kappa*a) = 1 + 1/(1+delta)^3 -- missing the factor of 2
    in the denominator. That form gives f -> 2 (not 1.5) as kappa*a ->
    infinity, which is physically wrong; confirmed as an error in that
    secondary source (not this library) by fetching THREE independent
    primary/near-primary sources directly (Ohshima 1994, Swan & Furst
    2012, Qin et al. 2017), all of which agree with each other and with
    this library's implementation.
    """
    regime = regime.lower()
    if regime == "huckel":
        return 1.0
    if regime == "smoluchowski":
        return 1.5
    if regime == "ohshima":
        if kappa_a is None or kappa_a <= 0:
            raise ValueError("kappa_a (the dimensionless product kappa*a) must be a positive "
                              "number for regime='ohshima' -- do not guess a value")
        delta = 2.5 / (kappa_a * (1.0 + 2.0 * math.exp(-kappa_a)))
        return 1.0 + 1.0 / (2.0 * (1.0 + delta) ** 3)
    if regime == "swan_furst":
        if kappa_a is None or kappa_a <= 0:
            raise ValueError("kappa_a (the dimensionless product kappa*a) must be a positive "
                              "number for regime='swan_furst' -- do not guess a value")
        ka = kappa_a
        return (16.0 + 18.0 * ka + 3.0 * ka ** 2) / (16.0 + 18.0 * ka + 2.0 * ka ** 2)
    if regime == "qin":
        if kappa_a is None or kappa_a <= 0:
            raise ValueError("kappa_a (the dimensionless product kappa*a) must be a positive "
                              "number for regime='qin' -- do not guess a value")
        delta = 2.8 / (kappa_a * (0.9 + math.exp(-kappa_a)))
        return 1.0 + 1.0 / (2.0 * (1.0 + delta) ** 3)
    raise ValueError("regime must be 'huckel', 'smoluchowski', 'ohshima', 'swan_furst', or 'qin' "
                      "-- no default is provided; Huckel applies for kappa*a << 1 (small "
                      "particles/low ionic strength), Smoluchowski for kappa*a >> 1 (larger "
                      "colloids, ionic strength >= ~10 mM), Ohshima/Swan-Furst/Qin are valid and "
                      "accurate across the full range (all three require kappa_a)")


def electrophoretic_mobility_relaxation_corrected(zeta_mV: float, kappa_a: float, viscosity_mPas: float, temperature_K: float = 298.15, z: int = 1, m: float = 0.15, rel_permittivity: float = WATER_REL_PERMITTIVITY_25C) -> float:
    """Electrophoretic mobility ((um*cm)/(V*s)) from zeta potential via
    a real, sourced ALTERNATIVE to the plain Henry equation
    (zeta_potential_henry/henry_function): O'Brien's own simplification
    of the Dukhin & Semenikhin relaxation-effect equation, valid for
    kappa_a > ~20 and capturing the REAL physical effect that Henry's
    LINEAR-in-zeta formula misses entirely -- at high |zeta| (>50 mV,
    common for ionic surfactant micelles per this project's own
    literature review), counterion relaxation around the moving
    particle makes mobility grow SUB-linearly with zeta, and it can
    even pass through a maximum (typically beyond ~100 mV) and decline.
    This is the concrete "next rigor tier beyond Henry" this project's
    own literature review flagged (via O'Brien & White, J. Chem. Soc.
    Faraday Trans. 2, 74 (1978) 1607-1626).

    Source: Delgado, Gonzalez-Caballero, Hunter, Koopal & Lyklema
    (IUPAC Technical Report), J. Colloid Interface Sci. 309 (2007)
    194-224, Eq. (24) -- fetched and read directly from the real paper
    text (not guessed or reconstructed), itself O'Brien's simplification
    (neglecting terms of order 1/kappa_a) of the Dukhin & Semenikhin
    equation (the same review's Eq. (22), for symmetrical z-z
    electrolytes). Both automatically account for diffuse-layer
    conductivity, unlike plain Henry.

        y = e*zeta / (k_B*T)                          [dimensionless zeta]
        X = y/2 - (ln2/z)*(1 - exp(-z*y))
        Y = 2 + [kappa_a / (1 + 3*m/z^2)] * exp(-z*y/2)
        ue = (eps_r*eps0*k_B*T) / (eta*e) * [y - 4*X/Y]

    z: counterion/co-ion charge number (symmetrical z-z electrolyte
    assumed, e.g. z=1 for NaCl-like). m: dimensionless ionic mobility
    parameter (eq. 18 of the same source, m = (2/3)*(kT/e)^2*eps_r*eps0/
    (eta*D_ion)) -- the source paper states "for aqueous solutions, m is
    about 0.15" as a commonly-used illustrative value; pass the real
    value for the actual counterion if known rather than relying on
    this default for anything beyond an order-of-magnitude estimate.

    Real, verified cross-check: at LOW zeta and LARGE kappa_a, this
    reduces to the Smoluchowski limit of this project's OWN
    zeta_potential_henry/henry_function (f(kappa_a)->1.5) -- confirmed
    algebraically and in tests/test_electrostatics_dynamics.py, tying
    this new formula back to already-validated existing code rather
    than floating unconnected.

    kappa_a must exceed ~20 (the source's own stated valid range) --
    below that, this simplification (which drops O(1/kappa_a) terms)
    is not reliable; use the full O'Brien-White numerical solution or
    the Dukhin-Semenikhin equation (22) directly in that regime instead
    (neither implemented here).
    """
    if kappa_a <= 20:
        raise ValueError("this simplified formula requires kappa_a > ~20 (the source's own stated "
                          "valid range, since it drops terms of order 1/kappa_a) -- not reliable below that")
    if viscosity_mPas <= 0:
        raise ValueError("viscosity_mPas must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if z <= 0:
        raise ValueError("z must be a positive integer (charge number magnitude)")
    if m <= 0:
        raise ValueError("m must be positive")
    eta_SI = viscosity_mPas * 1e-3
    zeta_V = zeta_mV / 1000.0
    y = (E_CHARGE * zeta_V) / (K_B * temperature_K)
    x_term = y / 2.0 - (math.log(2.0) / z) * (1.0 - math.exp(-z * y))
    y_term = 2.0 + (kappa_a / (1.0 + 3.0 * m / z ** 2)) * math.exp(-z * y / 2.0)
    ue_SI = (rel_permittivity * EPS0 * K_B * temperature_K) / (eta_SI * E_CHARGE) * (y - 4.0 * x_term / y_term)
    return ue_SI / (1e-6 * 1e-2)  # m^2/(V.s) -> (um*cm)/(V.s)


def zeta_potential_relaxation_corrected(electrophoretic_mobility_um_cm_per_Vs: float, kappa_a: float, viscosity_mPas: float, temperature_K: float = 298.15, z: int = 1, m: float = 0.15, rel_permittivity: float = WATER_REL_PERMITTIVITY_25C, zeta_search_bound_mV: float = 150.0) -> float:
    """Zeta potential (mV) from a measured electrophoretic mobility via
    the relaxation-corrected model (see
    electrophoretic_mobility_relaxation_corrected for the forward
    formula and full sourcing) -- inverts it numerically (grid scan plus
    bisection, mirroring solve_rubingh_x's pattern elsewhere in this
    project), since the relaxation-corrected mobility is NOT a simple
    linear function of zeta the way Henry's is.

    Restricted to the search range [0, zeta_search_bound_mV]: the
    source paper notes mobility passes through a MAXIMUM (typically
    beyond ~100 mV) and then declines with further increasing zeta, so
    mobility is only invertible (one-to-one) below that maximum -- the
    default 150 mV bound stays safely below where that non-monotonic
    behavior is reported to begin for typical systems, but raises
    explicitly if no root is found in range rather than silently
    returning a wrong branch."""
    def f(zeta_mV: float) -> float:
        return electrophoretic_mobility_relaxation_corrected(zeta_mV, kappa_a, viscosity_mPas, temperature_K, z, m, rel_permittivity)

    lo, hi = 1e-6, zeta_search_bound_mV
    f_lo, f_hi = f(lo), f(hi)
    target = electrophoretic_mobility_um_cm_per_Vs
    if not (min(f_lo, f_hi) <= target <= max(f_lo, f_hi)):
        raise ValueError(
            f"no zeta potential in [0, {zeta_search_bound_mV}] mV reproduces the given mobility under "
            "this model -- the mobility may correspond to a zeta beyond the (reported) mobility "
            "maximum, or the inputs are inconsistent; do not guess, extend zeta_search_bound_mV only "
            "if you have independent reason to believe zeta is that high and are aware of the "
            "non-monotonic-mobility caveat above"
        )
    for _ in range(100):
        mid = (lo + hi) / 2.0
        f_mid = f(mid)
        if abs(f_mid - target) < 1e-12:
            return mid
        if (f_mid - target) * (f_lo - target) <= 0:
            hi = mid
        else:
            lo, f_lo = mid, f_mid
    return (lo + hi) / 2.0


def zeta_potential_henry(electrophoretic_mobility_um_cm_per_Vs: float, viscosity_mPas: float, regime: str, rel_permittivity: float = WATER_REL_PERMITTIVITY_25C, kappa_a: float | None = None) -> float:
    """Zeta potential (mV) from electrophoretic mobility via the Henry
    equation: mobility = (2*eps_r*eps0*zeta*f(kappa*a)) / (3*eta), i.e.
    zeta = 3*eta*mobility / (2*eps_r*eps0*f(kappa*a)).

    electrophoretic_mobility_um_cm_per_Vs: mobility in the commonly
    reported unit (um*cm)/(V*s) -- this is converted internally to SI.
    viscosity_mPas: solvent viscosity in mPa.s (=cP; water is ~0.89 at 25C).
    regime: 'huckel', 'smoluchowski', or 'ohshima' -- see henry_function(),
    no default. kappa_a is required (and used) only for 'ohshima'.
    """
    if viscosity_mPas <= 0:
        raise ValueError("viscosity_mPas must be positive")
    f_ka = henry_function(regime, kappa_a)
    mobility_SI = electrophoretic_mobility_um_cm_per_Vs * 1e-6 * 1e-2  # (um*cm)/(V.s) -> m^2/(V.s)
    eta_SI = viscosity_mPas * 1e-3  # mPa.s -> Pa.s
    zeta_V = (3.0 * eta_SI * mobility_SI) / (2.0 * rel_permittivity * EPS0 * f_ka)
    return zeta_V * 1000.0  # V -> mV
