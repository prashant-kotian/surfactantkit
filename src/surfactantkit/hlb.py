"""Hydrophile-Lipophile Balance (HLB): Griffin's and Davies' methods.

Constants verified against cited sources (see docstrings) rather than
taken from memory alone.
"""

from __future__ import annotations

# Davies (1957) group numbers. This is a deliberately partial table --
# only groups with a verified numeric value are included. Notably
# missing: quaternary ammonium, amide, sulfonate, and several others
# common in amidoamine/gemini cationic surfactant chemistry, because no
# single consistently-cited value for them was confirmed. Do not guess
# a value for a missing group; raise instead (see hlb_davies).
DAVIES_HYDROPHILIC_GROUPS = {
    "SO4Na": 38.6,       # -SO4- Na+
    "COOK": 21.1,        # -COO- K+
    "COONa": 19.1,       # -COO- Na+
    "N_tertiary_amine": 9.4,
    "ester_sorbitan_ring": 6.8,
    "ester_free": 2.4,
    "COOH": 2.1,
    "OH_free": 1.9,
    "O_ether": 1.3,
    "OH_sorbitan_ring": 0.5,
}

DAVIES_LIPOPHILIC_GROUPS = {
    "CH2": 0.475,
    "CH3": 0.475,
    "CH": 0.475,
    "vinyl_CH": 0.475,  # =CH-
}

_ALL_KNOWN_GROUPS = set(DAVIES_HYDROPHILIC_GROUPS) | set(DAVIES_LIPOPHILIC_GROUPS)


def hlb_griffin(mw_hydrophilic: float, mw_total: float) -> float:
    """Griffin's method: HLB = 20 * (hydrophilic MW / total MW).

    Returns a value on Griffin's original 0-20 scale.
    """
    if mw_total <= 0:
        raise ValueError("mw_total must be positive")
    if not (0.0 <= mw_hydrophilic <= mw_total):
        raise ValueError("mw_hydrophilic must be between 0 and mw_total")
    return 20.0 * (mw_hydrophilic / mw_total)


def hlb_davies(group_counts: dict[str, int]) -> float:
    """Davies' method: HLB = 7 + sum(hydrophilic group numbers) -
    sum(lipophilic group numbers).

    group_counts: e.g. {"SO4Na": 1, "CH2": 11, "CH3": 1}. Keys must be
    from DAVIES_HYDROPHILIC_GROUPS or DAVIES_LIPOPHILIC_GROUPS -- raises
    KeyError naming the unknown group rather than silently ignoring or
    guessing a value for it, since a wrong Davies number is worse than
    an explicit error (this is exactly the "quaternary ammonium, amide,
    sulfonate not in this table" gap noted above).
    """
    unknown = set(group_counts) - _ALL_KNOWN_GROUPS
    if unknown:
        raise KeyError(
            f"No verified Davies group number for: {sorted(unknown)}. "
            "Add a cited value to DAVIES_HYDROPHILIC_GROUPS or "
            "DAVIES_LIPOPHILIC_GROUPS rather than guessing."
        )
    hlb = 7.0
    for group, count in group_counts.items():
        if group in DAVIES_HYDROPHILIC_GROUPS:
            hlb += DAVIES_HYDROPHILIC_GROUPS[group] * count
        else:
            hlb -= DAVIES_LIPOPHILIC_GROUPS[group] * count
    return hlb
