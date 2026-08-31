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


def henry_function(regime: str) -> float:
    """Henry function f(kappa*a) limiting values: 'huckel' = 1.0 (small
    particles / low ionic strength, kappa*a << 1), 'smoluchowski' = 1.5
    (large particles / ionic strength >= ~10 mM, kappa*a >> 1). There is
    no default -- picking the wrong regime silently is exactly the kind
    of error this library exists to prevent; the caller must state which
    applies for their particle size and ionic strength."""
    regime = regime.lower()
    if regime == "huckel":
        return 1.0
    if regime == "smoluchowski":
        return 1.5
    raise ValueError("regime must be 'huckel' or 'smoluchowski' -- no default is provided; "
                      "Huckel applies for kappa*a << 1 (small particles/low ionic strength), "
                      "Smoluchowski for kappa*a >> 1 (larger colloids, ionic strength >= ~10 mM)")


def zeta_potential_henry(electrophoretic_mobility_um_cm_per_Vs: float, viscosity_mPas: float, regime: str, rel_permittivity: float = WATER_REL_PERMITTIVITY_25C) -> float:
    """Zeta potential (mV) from electrophoretic mobility via the Henry
    equation: mobility = (2*eps_r*eps0*zeta*f(kappa*a)) / (3*eta), i.e.
    zeta = 3*eta*mobility / (2*eps_r*eps0*f(kappa*a)).

    electrophoretic_mobility_um_cm_per_Vs: mobility in the commonly
    reported unit (um*cm)/(V*s) -- this is converted internally to SI.
    viscosity_mPas: solvent viscosity in mPa.s (=cP; water is ~0.89 at 25C).
    regime: 'huckel' or 'smoluchowski' -- see henry_function(), no default.
    """
    if viscosity_mPas <= 0:
        raise ValueError("viscosity_mPas must be positive")
    f_ka = henry_function(regime)
    mobility_SI = electrophoretic_mobility_um_cm_per_Vs * 1e-6 * 1e-2  # (um*cm)/(V.s) -> m^2/(V.s)
    eta_SI = viscosity_mPas * 1e-3  # mPa.s -> Pa.s
    zeta_V = (3.0 * eta_SI * mobility_SI) / (2.0 * rel_permittivity * EPS0 * f_ka)
    return zeta_V * 1000.0  # V -> mV
