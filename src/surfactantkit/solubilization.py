"""Micellar solubilization capacity.

Originally implemented only the Molar Solubilization Ratio (MSR), which
has one clean, unambiguous definition in the literature -- the
micelle-water partition coefficient (Km) was deliberately left out,
because a check of the literature found its exact formula is study-
specific (some define it on a mole fraction basis, others on a molar-
concentration-ratio basis, with no single dominant convention), and
implementing one silently risked reproducing the wrong convention for a
given paper.

UPDATE, 2026-09-10 (alternative-methods sweep): re-checked, and the
ambiguity is real and still unresolved -- but Km is also a real,
commonly co-reported COMPLEMENTARY metric to MSR (this project's own
literature review: apparent solubilization capacity depends heavily on
WHERE the solute sits, palisade/interface vs. core -- MSR alone doesn't
capture this). Rather than continue omitting it entirely, this module
now implements the MOLE-FRACTION-BASIS convention specifically (X_mic/
X_w, one of the two real conventions found, not invented) with the
ambiguity explicitly disclosed in micelle_water_partition_coefficient's
own docstring -- letting a researcher use it with eyes open, matching
this project's "let a researcher choose the method suited to their
system" goal, rather than a blanket omission. Same discipline as the
missing Davies HLB groups and the unimplemented Krafft point for
anything that remains genuinely unresolved.
"""

from __future__ import annotations

WATER_MOLARITY_M = 55.5  # mol/L, same dilute-aqueous-solution convention as thermodynamics.py


def molar_solubilization_ratio(total_solubilized_M: float, intrinsic_water_solubility_M: float, surfactant_concentration_M: float, cmc_M: float) -> float:
    """Molar Solubilization Ratio (MSR, dimensionless): moles of
    solubilizate taken up per mole of MICELLIZED surfactant (i.e. only
    the surfactant above the CMC, which is what actually forms micelles):
    MSR = (total_solubilized - intrinsic_water_solubility) / (surfactant_concentration - cmc).

    All four inputs must be in the same concentration unit (e.g. all mM
    or all mol/L). surfactant_concentration must exceed cmc -- below the
    CMC there are no micelles to solubilize anything into.
    """
    if surfactant_concentration_M <= cmc_M:
        raise ValueError("surfactant_concentration_M must exceed cmc_M (no micelles form below the CMC)")
    if total_solubilized_M < intrinsic_water_solubility_M:
        raise ValueError("total_solubilized_M cannot be less than intrinsic_water_solubility_M")
    if intrinsic_water_solubility_M < 0:
        raise ValueError("intrinsic_water_solubility_M must be non-negative")
    return (total_solubilized_M - intrinsic_water_solubility_M) / (surfactant_concentration_M - cmc_M)


def micelle_water_partition_coefficient(
    total_solubilized_M: float,
    intrinsic_water_solubility_M: float,
    surfactant_concentration_M: float,
    cmc_M: float,
    water_molarity_M: float = WATER_MOLARITY_M,
) -> float:
    """Micelle-water partition coefficient (Km, dimensionless), MOLE-
    FRACTION-BASIS convention: Km = X_solubilizate,micelle /
    X_solubilizate,water, a real, commonly co-reported complementary
    metric to molar_solubilization_ratio -- captures how STRONGLY the
    solubilizate favors the micellar pseudophase over the bulk aqueous
    phase (an intensive, concentration-normalized quantity), whereas MSR
    is an extensive uptake-per-mole-of-surfactant quantity. A high Km
    with a low MSR (or vice versa) is real, meaningful information MSR
    alone does not capture (this project's own literature review: Ianiro
    et al., Langmuir 2019, PMC6448116 -- apparent solubilization
    capacity depends heavily on WHERE the solute sits, palisade/
    interface vs. core).

    IMPORTANT, disclosed honestly: this is ONE of (at least) two real
    conventions found in the literature for Km -- some papers instead
    use a molar-CONCENTRATION-RATIO basis (Km = C_solubilizate,micelle /
    C_solubilizate,water, no mole-fraction normalization), which gives a
    numerically DIFFERENT value from this mole-fraction-basis
    calculation for the same raw data. There is no single dominant
    convention -- verify which one a given paper you are comparing
    against actually uses before treating numbers as directly
    comparable; do not assume they match.

    Same four raw-data inputs as molar_solubilization_ratio (see that
    function for their meaning), reused directly for consistency:

        X_micelle = (total_solubilized - intrinsic_water_solubility)
                    / [(total_solubilized - intrinsic_water_solubility)
                       + (surfactant_concentration - cmc)]
        X_water = intrinsic_water_solubility / water_molarity_M
        Km = X_micelle / X_water

    water_molarity_M: dilute-water-molarity approximation (default 55.5
    mol/L, same convention as thermodynamics.cmc_to_mole_fraction).
    intrinsic_water_solubility_M is taken as the actual aqueous-phase
    solubilizate concentration at equilibrium (valid when excess
    undissolved solute is present, saturating the water phase at its own
    intrinsic solubility -- the same assumption molar_solubilization_ratio
    already makes for this same input).
    """
    if surfactant_concentration_M <= cmc_M:
        raise ValueError("surfactant_concentration_M must exceed cmc_M (no micelles form below the CMC)")
    if total_solubilized_M < intrinsic_water_solubility_M:
        raise ValueError("total_solubilized_M cannot be less than intrinsic_water_solubility_M")
    if intrinsic_water_solubility_M <= 0:
        raise ValueError("intrinsic_water_solubility_M must be positive (X_water is undefined at zero)")
    if water_molarity_M <= 0:
        raise ValueError("water_molarity_M must be positive")

    n_solubilized_micelle = total_solubilized_M - intrinsic_water_solubility_M
    n_surfactant_micellized = surfactant_concentration_M - cmc_M
    x_micelle = n_solubilized_micelle / (n_solubilized_micelle + n_surfactant_micellized)
    x_water = intrinsic_water_solubility_M / water_molarity_M
    return x_micelle / x_water
