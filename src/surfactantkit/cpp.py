"""Critical packing parameter (CPP) and Tanford's chain-length/volume
formulas.

Constants verified against the underlying per-CH2/CH3 volume formula
(Tanford, "The Hydrophobic Effect"): v(CH3) = 54.3 cubic-Angstrom,
v(CH2) = 26.9 cubic-Angstrom at 25C, giving v_total = v(CH3) + (nc-1)*v(CH2)
= 27.4 + 26.9*nc for an nc-carbon chain -- the commonly-cited shorthand
checks out exactly against the more detailed per-group formula, not
just taken on faith.
"""

from __future__ import annotations

_V_CH3 = 54.3   # cubic Angstrom, 25 C
_V_CH2 = 26.9   # cubic Angstrom, 25 C


def tanford_tail_volume(n_carbons: int) -> float:
    """Hydrophobic tail volume (cubic Angstrom) for a saturated,
    unbranched alkyl chain of n_carbons carbons. Equivalent to the
    shorthand 27.4 + 26.9*n_carbons.
    """
    if n_carbons < 1:
        raise ValueError("n_carbons must be at least 1")
    return _V_CH3 + (n_carbons - 1) * _V_CH2


def tanford_critical_length(n_carbons: int) -> float:
    """Maximum effective (critical) chain length (Angstrom) for a
    saturated, unbranched alkyl chain: lc = 1.5 + 1.265 * n_carbons."""
    if n_carbons < 1:
        raise ValueError("n_carbons must be at least 1")
    return 1.5 + 1.265 * n_carbons


def critical_packing_parameter(volume: float, head_area: float, length: float) -> float:
    """CPP = v / (a0 * lc). volume in cubic Angstrom, head_area in
    square Angstrom, length in Angstrom -- CPP itself is dimensionless."""
    if volume <= 0 or head_area <= 0 or length <= 0:
        raise ValueError("volume, head_area, and length must all be positive")
    return volume / (head_area * length)


def aggregation_number_spherical(tail_volume_A3: float, core_radius_A: float) -> float:
    """Estimated aggregation number for a spherical micelle from
    geometric packing: N_agg = V_core / v_tail, V_core = (4/3)*pi*r^3.

    tail_volume_A3: hydrophobic tail volume per surfactant (cubic
    Angstrom, from tanford_tail_volume). core_radius_A: micelle core
    radius (Angstrom) -- typically approximated as the critical chain
    length (tanford_critical_length) for a micelle at its maximum
    packing, but pass the actual value if known (e.g. from SANS/SAXS),
    since the two are not always equal in practice.

    This is a geometric estimate, not a substitute for aggregation
    numbers measured directly (e.g. by fluorescence quenching or
    scattering) -- treat it as a sanity-check order of magnitude.
    """
    import math

    if tail_volume_A3 <= 0:
        raise ValueError("tail_volume_A3 must be positive")
    if core_radius_A <= 0:
        raise ValueError("core_radius_A must be positive")
    core_volume = (4.0 / 3.0) * math.pi * core_radius_A ** 3
    return core_volume / tail_volume_A3


def classify_aggregate_morphology(cpp: float) -> str:
    """Expected aggregate morphology from the critical packing parameter.

    Thresholds: <=1/3 spherical micelles; (1/3, 1/2] cylindrical/rod
    micelles; (1/2, 1] vesicles/bilayers; >1 inverted structures.
    """
    if cpp <= 0:
        raise ValueError("cpp must be positive")
    if cpp <= 1.0 / 3.0:
        return "spherical micelle"
    if cpp <= 0.5:
        return "cylindrical/rodlike micelle"
    if cpp <= 1.0:
        return "vesicle/bilayer"
    return "inverted structure"
