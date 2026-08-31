"""Gibbs adsorption isotherm: surface excess concentration and minimum
area per molecule from a surface-tension-vs-log(concentration) slope."""

from __future__ import annotations

AVOGADRO = 6.02214076e23


def gibbs_gamma_max(slope_mn_per_ln_c: float, n_factor: float, temperature_k: float, r_gas: float = 8.314462618) -> float:
    """Maximum surface excess concentration (mol/m^2) from the pre-CMC
    slope of surface tension (mN/m) vs ln(concentration).

    slope_mn_per_ln_c: dGamma/d(ln C), in mN/m (negative for a surfactant
    lowering surface tension with increasing concentration).
    n_factor: Gibbs prefactor (1 for a nonionic surfactant or in the
    presence of excess electrolyte; 2 for a 1:1 ionic surfactant with no
    added salt, etc.) -- kept explicit because mixed ionic systems are
    model-sensitive; do not guess this value silently.
    """
    slope_n_per_m = slope_mn_per_ln_c / 1000.0
    return (-1.0 * slope_n_per_m) / (n_factor * r_gas * temperature_k)


def gibbs_a_min(gamma_max_mol_per_m2: float) -> float:
    """Minimum area per molecule (nm^2) from Gamma_max (mol/m^2)."""
    if gamma_max_mol_per_m2 <= 0:
        raise ValueError("Gamma_max must be positive to compute A_min")
    return 1.0e18 / (AVOGADRO * gamma_max_mol_per_m2)


def szyszkowski_surface_tension(concentration: float, gamma0_mN_m: float, gamma_max_mol_per_m2: float, K: float, temperature_K: float = 298.15, n_factor: float = 1.0, r_gas: float = 8.314462618) -> float:
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
    n_factor: same Gibbs prefactor as gibbs_gamma_max (1 nonionic/excess
    electrolyte, 2 for a 1:1 ionic surfactant with no added salt) --
    explicit, not auto-detected.

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
    term_N_per_m = n_factor * r_gas * temperature_K * gamma_max_mol_per_m2 * math.log(1.0 + K * concentration)
    term_mN_per_m = term_N_per_m * 1000.0
    return gamma0_mN_m - term_mN_per_m
