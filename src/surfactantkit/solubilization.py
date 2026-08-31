"""Micellar solubilization capacity.

Deliberately implements only the Molar Solubilization Ratio, which has
one clean, unambiguous definition in the literature. The micelle-water
partition coefficient (Km), despite being a commonly reported
companion quantity, is NOT implemented here: a check of the literature
found its exact formula is study-specific (some define it on a mole
fraction basis, others on a molar-concentration-ratio basis, with no
single dominant convention) -- implementing one silently would risk
reproducing the wrong convention for a given paper. Same discipline as
the missing Davies HLB groups and the unimplemented Krafft point.
"""

from __future__ import annotations


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
