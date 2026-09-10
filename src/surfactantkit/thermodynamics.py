"""Standard thermodynamics of micellization: Gibbs free energy, van't
Hoff enthalpy, entropy, and counterion binding degree.

Oddly absent from the rest of this library until now, despite being
central to the field -- these are the ΔG°mic/ΔH°mic/ΔS°mic triad that
essentially every experimental surfactant paper reports.
"""

from __future__ import annotations
import math
from dataclasses import dataclass

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


def mass_action_free_energy_of_micellization(cmc_mole_fraction: float, aggregation_number: float, temperature_K: float) -> float:
    """Standard Gibbs free energy of micellization, kJ/mol, via the
    MASS-ACTION-LAW model (single equilibrium n*S <-> M_n, K =
    [M_n]/[S]^n) -- a real, sourced ALTERNATIVE to
    gibbs_free_energy_micellization's PSEUDO-PHASE-SEPARATION
    convention, for NONIONIC surfactants only (see below):

    deltaG_mic = R*T * [(1 + 1/n)*ln(X_cmc) - (1/n)*ln(n)]

    Verified via an open-access source (Boneva-Aroca et al., RSC
    Advances 13 (2023) 9387, "Refined definition of the critical
    micelle concentration and application to alkyl maltosides used in
    membrane protein research," their eq. 21, quoted directly, not
    guessed or reconstructed) -- this project's own literature review
    had also flagged Rusanov, Langmuir 30 (2014) 14443-14451 as the
    definitive modern mass-action-law treatment, but that primary
    source is paywalled; the RSC Advances paper's eq. 21 is the SAME
    well-known classical single-equilibrium mass-action result (also
    widely cited elsewhere, e.g. Moroi's "Micelles" monograph and
    Rosen's "Surfactants and Interfacial Phenomena"), independently
    confirmed open-access rather than taken from a paywalled source.

    The pseudo-phase-separation model this project otherwise uses
    (gibbs_free_energy_micellization) is the n -> infinity LIMIT of
    this formula: as n grows, (1+1/n) -> 1 and (1/n)*ln(n) -> 0, so
    deltaG_mic -> R*T*ln(X_cmc) exactly (verified in
    tests/test_thermodynamics.py -- confirms this is a real
    generalization, not a competing/inconsistent formula). For small-
    to-moderate aggregation numbers (n below roughly 50-100, per this
    project's own literature review), the correction terms are NOT
    negligible and the pseudo-phase result can be a real, quantifiable
    over- or under-estimate -- this function makes that correction
    available directly, matching the "let a researcher choose the
    method suited to their system" goal.

    NONIONIC SURFACTANTS ONLY: this is the simple single-equilibrium
    mass-action result with no counterion association step. Extending
    to ionic surfactants requires a materially different mass-action
    treatment (an explicit counterion-binding equilibrium alongside
    monomer aggregation) that was not verified this pass -- do not
    apply this function to an ionic surfactant's counterion_factor
    convention; use gibbs_free_energy_micellization with the
    appropriate counterion_factor for that case instead.
    """
    if not (0.0 < cmc_mole_fraction < 1.0):
        raise ValueError("cmc_mole_fraction must be strictly between 0 and 1")
    if aggregation_number <= 1:
        raise ValueError("aggregation_number must be greater than 1")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    n = aggregation_number
    g_per_RT = (1.0 + 1.0 / n) * math.log(cmc_mole_fraction) - (1.0 / n) * math.log(n)
    return (R_GAS * temperature_K * g_per_RT) / 1000.0


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


def _solve_3x3(matrix: list[list[float]], vector: list[float]) -> list[float]:
    """Gaussian elimination with partial pivoting for a 3x3 linear
    system. Small, local, pure-Python -- avoids adding numpy as a
    dependency for a single fixed-size solve (SurfactantKit is
    zero-external-dependency by design)."""
    a = [row[:] + [vector[i]] for i, row in enumerate(matrix)]
    n = 3
    for col in range(n):
        pivot_row = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot_row][col]) < 1e-15:
            raise ValueError("temperature series is degenerate (singular fit matrix) -- "
                              "need at least 3 points with genuinely different temperatures")
        a[col], a[pivot_row] = a[pivot_row], a[col]
        for r in range(col + 1, n):
            factor = a[r][col] / a[col][col]
            for c in range(col, n + 1):
                a[r][c] -= factor * a[col][c]
    x = [0.0, 0.0, 0.0]
    for i in reversed(range(n)):
        x[i] = (a[i][n] - sum(a[i][j] * x[j] for j in range(i + 1, n))) / a[i][i]
    return x


@dataclass
class VantHoffMultiPointFitResult:
    delta_H_kJ_per_mol: float  # at reference_temperature_K
    delta_S_J_per_mol_K: float  # at reference_temperature_K
    delta_Cp_J_per_mol_K: float  # assumed constant over the fitted range
    reference_temperature_K: float
    r_squared: float
    n_points: int
    delta_G_fitted_kJ_per_mol: list[float]  # model-predicted deltaG at each input T, same order as input


def vant_hoff_multi_point_fit(
    cmc_mole_fractions: list[float],
    temperatures_K: list[float],
    counterion_factor: float = 1.0,
    reference_temperature_K: float | None = None,
) -> VantHoffMultiPointFitResult:
    """Full multi-point (3+ temperatures, 5+ recommended) Gibbs-Helmholtz
    fit for deltaH, deltaS, AND deltaCp simultaneously -- the real fix
    this project's own benchmark/METHOD_ALTERNATIVES_LITERATURE_REVIEW.md
    already recommended for vant_hoff_enthalpy's strictly-2-point
    limitation, per Kantonen, Henriksen & Gilson, BBA Gen. Subj. 1862
    (2018) 692-704 (doi:10.1016/j.bbagen.2017.11.020, their Eq. 4):

        deltaG(Tj) = deltaH(Tr) - Tj*deltaS(Tr)
                     + deltaCp(Tr) * [(Tj-Tr) - Tj*ln(Tj/Tr)]

    where deltaG(Tj) = counterion_factor * R * Tj * ln(X_cmc(Tj)) (same
    formula and counterion_factor convention as
    gibbs_free_energy_micellization), and Tr is a reference temperature.

    The source paper describes fitting deltaH(Ti)/deltaS(Ti)/deltaCp(Ti)
    by "nonlinear optimization" -- but for a FIXED reference temperature
    Tr, this equation is actually LINEAR in the three unknowns
    (deltaH, deltaS, deltaCp): the Tj-dependence lives entirely in the
    three regressor coefficients [1, -Tj, (Tj-Tr)-Tj*ln(Tj/Tr)], not in
    the parameters themselves. This implementation exploits that and
    solves the exact ordinary-least-squares closed form (via Gaussian
    elimination on the 3x3 normal-equations system) instead of iterative
    nonlinear optimization -- mathematically equivalent for a fixed Tr,
    verified here via exact round-trip recovery of synthetic data (see
    tests/test_thermodynamics.py), and avoids adding a nonlinear-solver
    dependency to a project that is zero-external-dependency by design.

    reference_temperature_K: defaults to the mean of temperatures_K,
    which keeps the fit well-conditioned (reduces the numerical
    correlation between the fitted deltaH and deltaCp) -- deltaH and
    deltaS come out "at" this reference temperature (a real, unavoidable
    consequence of deltaCp being nonzero: enthalpy and entropy genuinely
    differ at different reference temperatures, unlike deltaCp itself or
    the fitted deltaG(T) curve, both of which are reference-independent).
    Override to get deltaH/deltaS reported at a specific literature-
    comparable temperature (e.g. 298.15 K) instead.

    Needs at least 3 points (the fit has 3 free parameters) with
    genuinely different temperatures; 5+ is recommended for deltaCp to
    be reliably resolved from real (noisy) data, per the source paper's
    own experimental design (they use 5 temperature points).
    """
    if len(cmc_mole_fractions) != len(temperatures_K):
        raise ValueError("cmc_mole_fractions and temperatures_K must be the same length")
    if len(cmc_mole_fractions) < 3:
        raise ValueError("need at least 3 temperature points to fit deltaH, deltaS, and deltaCp simultaneously")
    if any(not (0.0 < x < 1.0) for x in cmc_mole_fractions):
        raise ValueError("all CMC mole fractions must be strictly between 0 and 1")
    if any(t <= 0 for t in temperatures_K):
        raise ValueError("all temperatures must be positive")

    n = len(temperatures_K)
    tr = reference_temperature_K if reference_temperature_K is not None else sum(temperatures_K) / n

    delta_g = [
        (counterion_factor * R_GAS * t * math.log(x)) / 1000.0  # kJ/mol, same as gibbs_free_energy_micellization
        for x, t in zip(cmc_mole_fractions, temperatures_K)
    ]

    # Regressors per point: r1=1 (coeff deltaH), r2=-Tj (coeff deltaS),
    # r3=(Tj-Tr)-Tj*ln(Tj/Tr) (coeff deltaCp); all in kJ-consistent units
    # (deltaS, deltaCp come out in kJ/(mol.K) here, converted to J after).
    r3 = [(t - tr) - t * math.log(t / tr) for t in temperatures_K]

    m00 = float(n)
    m01 = -sum(temperatures_K)
    m02 = sum(r3)
    m11 = sum(t * t for t in temperatures_K)
    m12 = -sum(t * g for t, g in zip(temperatures_K, r3))
    m22 = sum(g * g for g in r3)
    matrix = [[m00, m01, m02], [m01, m11, m12], [m02, m12, m22]]

    b0 = sum(delta_g)
    b1 = -sum(t * g for t, g in zip(temperatures_K, delta_g))
    b2 = sum(g * dg for g, dg in zip(r3, delta_g))
    vector = [b0, b1, b2]

    delta_h_kJ, delta_s_kJ, delta_cp_kJ = _solve_3x3(matrix, vector)

    fitted = [delta_h_kJ - t * delta_s_kJ + delta_cp_kJ * g for t, g in zip(temperatures_K, r3)]
    mean_g = sum(delta_g) / n
    ss_res = sum((y - f) ** 2 for y, f in zip(delta_g, fitted))
    ss_tot = sum((y - mean_g) ** 2 for y in delta_g)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

    return VantHoffMultiPointFitResult(
        delta_H_kJ_per_mol=delta_h_kJ,
        delta_S_J_per_mol_K=delta_s_kJ * 1000.0,
        delta_Cp_J_per_mol_K=delta_cp_kJ * 1000.0,
        reference_temperature_K=tr,
        r_squared=r_squared,
        n_points=n,
        delta_G_fitted_kJ_per_mol=fitted,
    )


def entropy_micellization(delta_g_mic_kJ_per_mol: float, delta_h_mic_kJ_per_mol: float, temperature_K: float) -> float:
    """Entropy of micellization, J/(mol.K): deltaS = (deltaH - deltaG) / T.
    Completes the deltaG/deltaH/deltaS triad given the other two."""
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    delta_s_kJ = (delta_h_mic_kJ_per_mol - delta_g_mic_kJ_per_mol) / temperature_K
    return delta_s_kJ * 1000.0  # kJ/(mol.K) -> J/(mol.K)
