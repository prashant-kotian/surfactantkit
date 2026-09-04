"""Gibbs adsorption isotherm: surface excess concentration and minimum
area per molecule from a surface-tension-vs-log(concentration) slope."""

from __future__ import annotations

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
    import math

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
