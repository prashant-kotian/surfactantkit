"""Standard thermodynamics of micellization: Gibbs free energy, van't
Hoff enthalpy, entropy, and counterion binding degree.

Oddly absent from the rest of this library until now, despite being
central to the field -- these are the ΔG°mic/ΔH°mic/ΔS°mic triad that
essentially every experimental surfactant paper reports.
"""

from __future__ import annotations
import math

R_GAS = 8.314462618  # J/(mol.K)
WATER_MOLARITY_M = 55.5  # mol/L, standard dilute-aqueous-solution approximation


def cmc_to_mole_fraction(cmc_M: float, water_molarity_M: float = WATER_MOLARITY_M) -> float:
    """Convert a CMC in molarity (mol/L) to the mole-fraction scale used
    in the standard Gibbs free energy of micellization formula:
    X_cmc = CMC / (CMC + water_molarity) ~= CMC / water_molarity for
    dilute solutions (the usual approximation, water_molarity ~= 55.5
    mol/L for water at room temperature)."""
    if cmc_M <= 0:
        raise ValueError("cmc_M must be positive")
    return cmc_M / (cmc_M + water_molarity_M)


def counterion_binding_degree(slope_below_cmc: float, slope_above_cmc: float) -> float:
    """Degree of counterion binding to the micelle, beta, from the
    slope-ratio method applied to a conductivity-vs-concentration plot:
    beta = 1 - (slope_above_cmc / slope_below_cmc).

    Both slopes should be positive (conductivity increases with
    surfactant concentration in both regimes, just at different rates).
    Note: the physical interpretation of this slope-ratio method has
    been questioned for some systems in the literature (the arithmetic
    is standard and widely used, but treat beta from this method as an
    approximate, commonly-reported value, not an unimpeachable one --
    see literature_validation_notes.md).
    """
    if slope_below_cmc <= 0 or slope_above_cmc <= 0:
        raise ValueError("both slopes must be positive")
    if slope_above_cmc >= slope_below_cmc:
        raise ValueError(
            "slope_above_cmc should be less than slope_below_cmc (conductivity "
            "rises more slowly above the CMC as counterions associate with "
            "micelles) -- check which slope is which"
        )
    return 1.0 - (slope_above_cmc / slope_below_cmc)


def gibbs_free_energy_micellization(cmc_mole_fraction: float, temperature_K: float, counterion_factor: float = 1.0) -> float:
    """Standard Gibbs free energy of micellization, kJ/mol:
    deltaG_mic = counterion_factor * R * T * ln(X_cmc).

    counterion_factor: explicit, not guessed. Use 1.0 for a nonionic
    surfactant. For an ionic surfactant with counterion binding degree
    beta (see counterion_binding_degree), use (2 - beta) -- this is
    NOT auto-detected from the inputs, because silently assuming
    nonionic behavior for an ionic surfactant (or vice versa) is exactly
    the kind of hidden-assumption error this library exists to prevent.
    """
    if not (0.0 < cmc_mole_fraction < 1.0):
        raise ValueError("cmc_mole_fraction must be strictly between 0 and 1")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    return (counterion_factor * R_GAS * temperature_K * math.log(cmc_mole_fraction)) / 1000.0


def vant_hoff_enthalpy(cmc1_mole_fraction: float, temperature1_K: float, cmc2_mole_fraction: float, temperature2_K: float) -> float:
    """Van't Hoff enthalpy of micellization, kJ/mol, from CMC (mole
    fraction scale) measured at two temperatures:
    deltaH_mic = R * (ln(X_cmc,2) - ln(X_cmc,1)) / (1/T2 - 1/T1).

    Real sign bug found and fixed 2026-09-04: this function previously
    divided by (1/T1 - 1/T2) instead of (1/T2 - 1/T1), giving the exact
    opposite sign of the true van't Hoff enthalpy. Caught via two
    independent, convergent checks: (1) an unaugmented Claude Opus 4.8
    pilot answer for a real SurfBench question matched textbook Gibbs-
    Helmholtz-derivation and Le Chatelier physical reasoning exactly
    (CMC decreasing with rising T => endothermic, positive dH) while
    this function's old output had the opposite sign; (2) re-examining
    this project's own existing literature_validation_notes.md entry for
    this exact function against Fu et al. 2019 (RSC Advances) showed the
    old two-point estimates (+6.446, +0.462 kJ/mol) had the opposite sign
    from the paper's own local-derivative value (-4.472 kJ/mol) -- at the
    time misread as pure two-point-vs-polynomial magnitude noise, but the
    sign mismatch was real and is this same bug. Full re-derivation from
    the Gibbs-Helmholtz equation (d(dG/T)/dT = -dH/T^2, with
    dG = R*T*ln(X_cmc)) gives d(ln X_cmc) = (dH/R) d(1/T), i.e. the
    (1/T2 - 1/T1) denominator used now, not the old (1/T1 - 1/T2).

    Assumes deltaH is constant over the T1-T2 interval -- this
    assumption weakens for wide temperature ranges or when the
    aggregation number itself varies significantly with temperature;
    treat results from a wide T range with appropriate caution.
    """
    if temperature1_K <= 0 or temperature2_K <= 0:
        raise ValueError("temperatures must be positive")
    if temperature1_K == temperature2_K:
        raise ValueError("temperature1_K and temperature2_K must differ")
    if not (0.0 < cmc1_mole_fraction < 1.0) or not (0.0 < cmc2_mole_fraction < 1.0):
        raise ValueError("CMC mole fractions must be strictly between 0 and 1")
    delta_ln_cmc = math.log(cmc2_mole_fraction) - math.log(cmc1_mole_fraction)
    delta_inv_T = (1.0 / temperature2_K) - (1.0 / temperature1_K)
    return (R_GAS * delta_ln_cmc / delta_inv_T) / 1000.0


def entropy_micellization(delta_g_mic_kJ_per_mol: float, delta_h_mic_kJ_per_mol: float, temperature_K: float) -> float:
    """Entropy of micellization, J/(mol.K): deltaS = (deltaH - deltaG) / T.
    Completes the deltaG/deltaH/deltaS triad given the other two."""
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    delta_s_kJ = (delta_h_mic_kJ_per_mol - delta_g_mic_kJ_per_mol) / temperature_K
    return delta_s_kJ * 1000.0  # kJ/(mol.K) -> J/(mol.K)
