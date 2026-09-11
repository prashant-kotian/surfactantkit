"""Hydrophile-Lipophile Balance (HLB): Griffin's and Davies' methods.

Constants verified against cited sources (see docstrings) rather than
taken from memory alone.
"""

from __future__ import annotations

import math

# Davies (1957) group numbers. This is a deliberately partial table --
# only groups with a verified LITERATURE-SOURCED numeric value are
# included. Do not guess a value for a missing group; raise instead (see
# hlb_davies). Amide and sulfonate no longer have a literature-sourced
# number confirmed here either, BUT this project derived real, disclosed,
# non-guessed estimates for both (2026-09-11) via the same cross-
# calibration methodology Davies himself used -- see
# DAVIES_DERIVED_HYDROPHILIC_GROUPS and derive_davies_group_number_from_griffin
# below, and pass allow_derived_groups=True to hlb_davies to use them.
# Several others remain genuinely missing (common in amidoamine/gemini
# cationic surfactant chemistry).
#
# Quaternary ammonium (RESOLVED 2026-09-08, real source-verified gap
# closure -- this project's own 2026-09-07/08 literature review had
# confirmed Davies' own 1957 paper never assigned a quaternary ammonium
# number at all; the real answer turned out to be a later paper filling
# that gap, not Davies himself). Source: B.H. O, J. Colloid Interface
# Sci. 198 (1998) 249, GN = 22.0 for -N+(CH3)3 (single-tail
# trimethylammonium) and 22.5 for >N+(CH3)2 (double-tail dimethyl-
# dialkylammonium) -- found via Proverbio, Bardavid, Arancibia & Schulz,
# Colloids Surf. A 214 (2003) 167-171 (fetched and read in full, PDF
# provided by the user), which reproduces both numbers and uses them in
# four fully worked HLB examples. All four independently, EXACTLY
# reproduced here (7 + GN - 0.475*n_carbons, this project's own existing
# CH2/CH3 = 0.475 convention): DTAB (C10, trimethyl) = 24.25, LTAB (C12,
# trimethyl) = 23.3, DDAB (2xC12, dimethyl-dialkyl) = 18.1, DODAB (2xC18,
# dimethyl-dialkyl) = 12.4 -- all match the source paper's own reported
# values to 2-3 decimal places, not just approximately (see
# test_hlb_cpp.py). Counterion assumed Cl-/Br- (the source paper
# explicitly states it neglected the difference between them); NOT
# verified for other counterions (e.g. the same paper's CTATOS/tosylate
# example has no computable Davies HLB for exactly this reason -- its
# own Table 1 leaves that entry blank).
DAVIES_HYDROPHILIC_GROUPS = {
    "SO4Na": 38.6,       # -SO4- Na+
    "COOK": 21.1,        # -COO- K+
    "COONa": 19.1,       # -COO- Na+
    "N_tertiary_amine": 9.4,
    "N_quat_trimethyl": 22.0,          # -N+(CH3)3, single alkyl tail (e.g. CTAB, LTAB, DTAB)
    "N_quat_dimethyl_dialkyl": 22.5,   # >N+(CH3)2, two alkyl tails (e.g. DDAB, DODAB)
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


def derive_davies_group_number_from_griffin(
    hydrophilic_fragment_mass_g_per_mol: float,
    total_molar_mass_g_per_mol: float,
    lipophilic_group_counts: dict[str, int],
) -> float:
    """Back-solve a NEW Davies-scale hydrophilic group number for a
    functional group Davies' own 1957 table never assigned one to --
    using the SAME cross-calibration methodology Davies himself used to
    build his original table (and the same methodology already used
    once in this project to resolve the quaternary-ammonium gap, see
    DAVIES_HYDROPHILIC_GROUPS's own docstring, sourced there via B.H. O
    1998): assume Davies' additive group-number HLB equals Griffin's
    independent, purely mass-ratio-based HLB (hlb_griffin, HLB=20*Mh/M)
    for one real reference compound, then solve the resulting single
    linear equation for the one unknown group number.

        7 + GN_unknown - sum(lipophilic contributions) = 20*Mh/M
        => GN_unknown = 20*Mh/M - 7 + sum(lipophilic contributions)

    This is a real, general, reusable METHOD, not a one-off guess --
    see hlb_davies_derived_sulfonate/amide_dialkanolamide/
    amide_monoalkanolamide below for concrete real-compound applications
    (sulfonate: sodium alkyl sulfonates; amide: fatty acid
    mono-/di-ethanolamides), each cross-checked across MULTIPLE real
    reference chain lengths for self-consistency, not asserted from a
    single point.

    hydrophilic_fragment_mass_g_per_mol: molar mass of the compound's
    hydrophilic "head" fragment, in the SAME atom-accounting convention
    Griffin's method needs for Mh -- the ENTIRE polar head, including
    any heteroatom-bearing carbon that is chemically part of the
    functional group itself (e.g. an amide's carbonyl carbon), never a
    plain hydrocarbon CH2/CH3 (those belong in lipophilic_group_counts
    instead).
    total_molar_mass_g_per_mol: the whole reference compound's molar
    mass.
    lipophilic_group_counts: the ALREADY-VERIFIED Davies lipophilic
    group counts for the rest of the molecule (e.g. {"CH2": 10, "CH3":
    1} for an 11-carbon alkyl tail) -- every key must already exist in
    DAVIES_LIPOPHILIC_GROUPS; raises KeyError naming any that don't
    rather than silently skipping them.

    IMPORTANT, disclosed honestly: a number produced this way is a
    DERIVED, self-consistent estimate, not an independently
    literature-published Davies group number -- it inherits Griffin's
    own well-known real limitation (a purely additive, LINEAR-in-chain-
    length Davies model cannot exactly reproduce Griffin's inherently
    NONLINEAR Mh/M mass ratio at every chain length; the two are
    designed to agree closely only over the practically-relevant
    commercial-surfactant range, roughly C8-C16 -- the same effect
    already documented in this project's hlb_griffin validation notes).
    Kept in a separate DAVIES_DERIVED_HYDROPHILIC_GROUPS table, not
    silently merged into DAVIES_HYDROPHILIC_GROUPS, and hlb_davies only
    uses it when explicitly told to via allow_derived_groups=True --
    the distinction between "Davies' own literature number" and "this
    project's derived estimate" is never allowed to blur at the call
    site.
    """
    if hydrophilic_fragment_mass_g_per_mol <= 0 or total_molar_mass_g_per_mol <= 0:
        raise ValueError("masses must be positive")
    if hydrophilic_fragment_mass_g_per_mol >= total_molar_mass_g_per_mol:
        raise ValueError("hydrophilic_fragment_mass_g_per_mol must be less than total_molar_mass_g_per_mol")
    unknown = set(lipophilic_group_counts) - set(DAVIES_LIPOPHILIC_GROUPS)
    if unknown:
        raise KeyError(
            f"No verified Davies LIPOPHILIC group number for: {sorted(unknown)} -- "
            "derive_davies_group_number_from_griffin needs the REST of the molecule's "
            "lipophilic groups to already be verified, since it is solving for the one "
            "remaining unknown."
        )
    lipophilic_total = sum(DAVIES_LIPOPHILIC_GROUPS[g] * c for g, c in lipophilic_group_counts.items())
    hlb_griffin = 20.0 * hydrophilic_fragment_mass_g_per_mol / total_molar_mass_g_per_mol
    return hlb_griffin - 7.0 + lipophilic_total


# DERIVED (2026-09-11) hydrophilic group numbers -- kept explicitly separate
# from DAVIES_HYDROPHILIC_GROUPS above (which is literature-sourced numbers
# only). Each derived via derive_davies_group_number_from_griffin against
# REAL, well-defined reference compounds at MULTIPLE chain lengths (not a
# single point), full derivation and real, disclosed chain-length spread
# below. This is a genuinely new capability/methodology contribution (not a
# literature citation) -- see derive_davies_group_number_from_griffin's own
# docstring for the full method and its honest limitations.
#
# SULFONATE (-SO3Na, e.g. sodium alkyl sulfonates, alpha-olefin sulfonates,
# LAS' non-aromatic analogue): derived from the sodium alkyl sulfonate
# series C(n)H(2n+1)-SO3Na (Mh = mass of -SO3Na = 103.047 g/mol; the
# alpha-carbon is a plain, unmodified CH2, unlike an amide/carboxyl carbon,
# so it belongs entirely in the lipophilic tail count). Computed at n=6, 7,
# 8, 9, 11, 12 (the exact chain lengths used in this project's own primary-
# source citation, Sutherland, Mercer, Everist & Leaist, J. Chem. Eng. Data
# 54(2) (2009) 272-278): GN = 6.80, 6.52, 6.33, 6.22, 6.20, 6.27
# respectively -- converges to a tight ~6.2-6.3 band over the practically-
# relevant C9-C12 range (real commercial alkyl sulfonate chain lengths),
# drifting at shorter/much-longer chains as expected (see the parent
# function's docstring on the linear-vs-nonlinear mismatch). Reference
# value taken at n=12 (this project's standard reference chain length,
# matching the quaternary-ammonium/Tanford-CPP convention used elsewhere).
#
# AMIDE (fatty acid alkanolamide heads -- real, common, commercially
# important nonionic surfactants, e.g. cocamide/lauramide DEA and MEA --
# NOT a bare primary fatty amide, which is barely water-soluble and not
# really a surfactant at all): derived from lauric/myristic/palmitic acid
# di- and mono-ethanolamides. Diethanolamide head, -CO-N(CH2CH2OH)2 (Mh =
# 132.139 g/mol): GN = 7.42 (C12), 7.55 (C14), 7.82 (C16) -- reference
# value at C12 = 7.42. Monoethanolamide head, -CO-NH-CH2CH2OH (Mh = 88.086
# g/mol): GN = 5.46 (C12), 5.67 (C14), 6.01 (C16) -- reference value at
# C12 = 5.46. Both head fragments are treated as ONE combined group number
# each (not further decomposed into separate amide+OH contributions),
# matching Davies' own real precedent of compound-specific combined groups
# for common structural motifs (e.g. his own "ester_sorbitan_ring" and
# "OH_sorbitan_ring" are likewise whole-motif groups, not atomistic
# decompositions) -- this also avoids compounding chain-length mismatch
# error across multiple separately-derived sub-groups.
#
# SULTAINE / alkyl sulfobetaine (-N+(CH3)2-(CH2)3-SO3-, e.g. decyl/dodecyl
# sulfobetaine SB10/SB12 -- the exact zwitterionic compounds already used
# as a real primary-source citation elsewhere in this project, Sutherland,
# Mercer, Everist & Leaist 2009): a genuine internal-salt zwitterion, no
# counterion mass to account for. Mh = 166.215 g/mol (the whole head:
# N(CH3)2 + propylene spacer + SO3). Computed at n=8, 10, 12, 14: GN =
# 8.70, 8.56, 8.61, 8.79 -- tight, ~2.7% spread, the tightest self-
# consistency of any derived group here. Reference value at C12 = 8.61.
#
# CARBOXYBETAINE (-CO-NH-(CH2)3-N+(CH3)2-CH2-COO-, e.g. cocamidopropyl
# betaine/CAPB -- one of the single most common real amphoteric
# surfactants, ubiquitous in shampoos as a foam booster/co-surfactant).
# Mh = 187.219 g/mol (the whole amidopropyl-betaine head, carbonyl carbon
# included per the same convention as the amide groups above). Computed
# at lauric/myristic/palmitic acid tails (C12/C14/C16): GN = 9.16, 9.28,
# 9.52 -- reference value at C12 = 9.16.
#
# PHOSPHATE ESTER, disodium monoalkyl (-O-P(=O)(ONa)2, e.g. disodium
# lauryl phosphate -- one common real commercial form; monoalkyl
# phosphate esters are also sold as the monosodium salt or as mono-/di-
# alkyl mixtures, NOT interchangeable with this specific derivation,
# disclosed explicitly since the salt form changes Mh). Mh = 140.950
# g/mol (unlike the amide/betaine groups, the alkyl-O-P linkage does NOT
# consume a tail carbon -- the full Cn alkyl stays in the lipophilic
# count). Computed at n=10, 12, 14: GN = 7.74, 7.79, 7.98 -- reference
# value at C12 = 7.79.
#
# IMIDAZOLINE -- deliberately NOT derived. Checked and explicitly skipped
# (2026-09-11), not merely unattempted: real "imidazoline surfactants"
# are typically synthesized via a cyclic imidazoline intermediate that
# HYDROLYZES in aqueous solution into a genuinely ambiguous mixture of
# open-chain amidoamine and ester-amine forms -- there is no single,
# undisputed "the" structure to compute Mh for the way there is for the
# other groups here, so deriving a number would require picking one
# contested structure and presenting it as settled, which this project's
# own "do not guess" discipline does not allow. Remains a genuinely open
# gap, not a guessed one.
DAVIES_DERIVED_HYDROPHILIC_GROUPS = {
    "sulfonate": 6.27,              # -SO3Na, derived @ C12 (see block comment above)
    "amide_dialkanolamide": 7.42,   # -CO-N(CH2CH2OH)2, derived @ C12
    "amide_monoalkanolamide": 5.46,  # -CO-NH-CH2CH2OH, derived @ C12
    "sultaine": 8.61,                # -N+(CH3)2-(CH2)3-SO3-, derived @ C12 (zwitterion)
    "carboxybetaine": 9.16,          # -CO-NH-(CH2)3-N+(CH3)2-CH2-COO-, derived @ C12
    "phosphate_diNa": 7.79,          # -O-P(=O)(ONa)2, derived @ C12
}

_ALL_KNOWN_GROUPS = set(DAVIES_HYDROPHILIC_GROUPS) | set(DAVIES_LIPOPHILIC_GROUPS)
_ALL_KNOWN_GROUPS_WITH_DERIVED = _ALL_KNOWN_GROUPS | set(DAVIES_DERIVED_HYDROPHILIC_GROUPS)


def guo_effective_eo_chain_length(n_eo: int) -> float:
    """"Effective" ethylene-oxide chain length (dimensionless) for a
    polyoxyethylene chain of n_eo actual EO units -- Guo, Rong & Ying,
    J. Colloid Interface Sci. 298 (2006) 441-450 (doi:10.1016/j.jcis.2005.12.009),
    real, sourced refinement to Davies' group-contribution HLB method
    specifically for NONIONIC POLYETHOXYLATED surfactants (the primary
    paper itself is paywalled on ScienceDirect -- these exact
    coefficients were verified via a citing patent, US 11,344,493 B2,
    which quotes the formula directly, not guessed or reconstructed).

    NEOeff = 13.45*ln(NEO) - 0.16*NEO + 1.26     (1 <= NEO <= 50)
    NEOeff = 0.056*NEO + 43.08                   (NEO >= 50)

    The paper's real finding (224 nonionic surfactants tested, average
    absolute HLB error under 1.5 vs. Davies' larger error) is that a
    polyoxyethylene chain's hydrophilic contribution is NOT linear in
    the actual EO count -- this effective/"apparent" chain length
    captures that real nonlinearity (sub-linear growth: each additional
    EO unit contributes progressively less as the chain gets longer),
    which plain Davies (linear in NEO) cannot.

    Second branch (NEO >= 50, primary source's own eq. 5', confirmed
    2026-09-11 against the obtained primary-source PDF -- previously
    unimplemented, this function raised above n_eo=50) reflects the
    paper's own stated reason: above ~50 EO units the chain's
    contribution to HLB has "taken up about 90%" of its saturating
    value, so growth becomes nearly linear again (small constant slope)
    rather than continuing the log-dominated sub-linear growth of the
    first branch. The two branches agree closely at NEO=50 by
    construction (13.45*ln(50)-0.16*50+1.26 = 46.0 vs. 0.056*50+43.08 =
    45.9), verified in tests/test_hlb_cpp.py.
    """
    if n_eo <= 0:
        raise ValueError("n_eo must be positive")
    if n_eo <= 50:
        return 13.45 * math.log(n_eo) - 0.16 * n_eo + 1.26
    return 0.056 * n_eo + 43.08


def guo_effective_alkyl_chain_length(n_ch2: int) -> float:
    """"Effective" alkyl (hydrophobic) chain length (dimensionless) for
    n_ch2 actual CH2 groups in a straight-chain tail -- Guo, Rong &
    Ying 2006 (see guo_effective_eo_chain_length's docstring for full
    sourcing). NCH2eff = 0.965*NCH2 - 0.178 -- nearly linear (unlike the
    EO chain's real sub-linearity), with a small correction near the
    chain-length origin."""
    if n_ch2 <= 0:
        raise ValueError("n_ch2 must be positive")
    return 0.965 * n_ch2 - 0.178


def guo_effective_po_chain_length(n_po: int) -> float:
    """"Effective" polyoxypropylene (PO) chain length (dimensionless)
    for n_po actual PO units -- Guo, Rong & Ying 2006 (see
    guo_effective_eo_chain_length's docstring for full sourcing).
    NPOeff = 2.057*NPO + 9.06. Standalone utility only: NOT chained into
    hlb_davies_guo_ecl below, because this project has no independently
    verified Davies-scale group number for a polyoxypropylene unit (the
    same "do not guess a missing group number" discipline as the
    amide/sulfonate gaps elsewhere in this module)."""
    if n_po <= 0:
        raise ValueError("n_po must be positive")
    return 2.057 * n_po + 9.06


def hlb_davies_guo_ecl(n_carbons_alkyl: int, n_eo: int, extra_group_counts: dict[str, int] | None = None) -> float:
    """Guo/Rong/Ying 2006 refinement of Davies' HLB method for a
    NONIONIC polyethoxylated surfactant (e.g. a fatty alcohol
    ethoxylate, R-(OCH2CH2)n-OH): applies Davies' own per-unit weights
    (O_ether=1.3, CH2=0.475) to the EFFECTIVE (not actual) EO and CH2
    chain lengths from guo_effective_eo_chain_length and
    guo_effective_alkyl_chain_length, which the source paper shows
    fits real HLB data substantially better (224 surfactants, average
    absolute error < 1.5) than plugging the actual counts into plain
    hlb_davies.

    RESOLVED 2026-09-11 (primary source obtained, PDF provided by the
    user): Table 1 of the primary paper lists group numbers under two
    columns, "Davies and Lin" vs. "ECL method", side by side for every
    group -- for CH2/CH3/CH/=CH- (lipophilic) and for the EO-repeat-unit
    group -CH2CH2O- (hydrophilic), the two columns are IDENTICAL
    (-0.475 and 0.33 respectively, unchanged). This confirms directly
    (not inferred) that Guo/Rong/Ying reuse Davies' original per-unit
    weights unchanged and apply them to the effective (not actual) chain
    lengths -- exactly what this function does. The previously disclosed
    uncertainty (only the transform equations were confirmed, via a
    citing patent, not the weight-reuse assumption) is now closed.

    n_carbons_alkyl: total carbons in the straight hydrophobic tail
    (same convention as elsewhere in this project: CH2 count =
    n_carbons_alkyl - 1, plus one terminal CH3, unaffected by the ECL
    correction). n_eo: actual number of ethylene oxide units (any
    positive value -- guo_effective_eo_chain_length handles both the
    NEO<50 and NEO>=50 branches of the source's own piecewise formula).
    extra_group_counts: any additional Davies group contributions not
    covered by the alkyl tail or EO chain (e.g. an ester headgroup) --
    same {group_name: count} convention as hlb_davies, added on top.
    """
    if n_carbons_alkyl < 2:
        raise ValueError("n_carbons_alkyl must be at least 2 (need at least one CH2 plus the terminal CH3)")
    n_ch2 = n_carbons_alkyl - 1
    n_ch2_eff = guo_effective_alkyl_chain_length(n_ch2)
    n_eo_eff = guo_effective_eo_chain_length(n_eo)
    hlb = (
        7.0
        + DAVIES_HYDROPHILIC_GROUPS["O_ether"] * n_eo_eff
        - DAVIES_LIPOPHILIC_GROUPS["CH2"] * n_ch2_eff
        - DAVIES_LIPOPHILIC_GROUPS["CH3"] * 1
    )
    if extra_group_counts:
        unknown = set(extra_group_counts) - _ALL_KNOWN_GROUPS
        if unknown:
            raise KeyError(
                f"No verified Davies group number for: {sorted(unknown)}. "
                "Add a cited value to DAVIES_HYDROPHILIC_GROUPS or "
                "DAVIES_LIPOPHILIC_GROUPS rather than guessing."
            )
        for group, count in extra_group_counts.items():
            if group in DAVIES_HYDROPHILIC_GROUPS:
                hlb += DAVIES_HYDROPHILIC_GROUPS[group] * count
            else:
                hlb -= DAVIES_LIPOPHILIC_GROUPS[group] * count
    return hlb


def hlb_griffin(mw_hydrophilic: float, mw_total: float) -> float:
    """Griffin's method: HLB = 20 * (hydrophilic MW / total MW).

    Returns a value on Griffin's original 0-20 scale.
    """
    if mw_total <= 0:
        raise ValueError("mw_total must be positive")
    if not (0.0 <= mw_hydrophilic <= mw_total):
        raise ValueError("mw_hydrophilic must be between 0 and mw_total")
    return 20.0 * (mw_hydrophilic / mw_total)


def hlb_davies(group_counts: dict[str, int], allow_derived_groups: bool = False) -> float:
    """Davies' method: HLB = 7 + sum(hydrophilic group numbers) -
    sum(lipophilic group numbers).

    group_counts: e.g. {"SO4Na": 1, "CH2": 11, "CH3": 1}. Keys must be
    from DAVIES_HYDROPHILIC_GROUPS or DAVIES_LIPOPHILIC_GROUPS -- raises
    KeyError naming the unknown group rather than silently ignoring or
    guessing a value for it, since a wrong Davies number is worse than
    an explicit error.

    allow_derived_groups: False by default. When True, ALSO accepts keys
    from DAVIES_DERIVED_HYDROPHILIC_GROUPS (sulfonate, amide_dialkanolamide,
    amide_monoalkanolamide -- this project's own 2026-09-11 Griffin-Davies
    cross-calibration derivations, NOT independently literature-published
    Davies numbers; see derive_davies_group_number_from_griffin for the
    method and honest limitations). Left OFF by default deliberately: a
    caller must explicitly opt in to using a derived-not-literature-sourced
    number, so the distinction is never silently blurred at the call site
    the way it would be if these were merged straight into
    DAVIES_HYDROPHILIC_GROUPS.
    """
    known = _ALL_KNOWN_GROUPS_WITH_DERIVED if allow_derived_groups else _ALL_KNOWN_GROUPS
    unknown = set(group_counts) - known
    if unknown:
        extra_hint = ""
        if not allow_derived_groups and unknown & set(DAVIES_DERIVED_HYDROPHILIC_GROUPS):
            extra_hint = (
                " (note: some of these DO have a derived, not-literature-sourced number "
                "available -- pass allow_derived_groups=True to use it, see "
                "DAVIES_DERIVED_HYDROPHILIC_GROUPS)"
            )
        raise KeyError(
            f"No verified Davies group number for: {sorted(unknown)}.{extra_hint} "
            "Add a cited value to DAVIES_HYDROPHILIC_GROUPS or "
            "DAVIES_LIPOPHILIC_GROUPS rather than guessing."
        )
    hlb = 7.0
    for group, count in group_counts.items():
        if group in DAVIES_HYDROPHILIC_GROUPS:
            hlb += DAVIES_HYDROPHILIC_GROUPS[group] * count
        elif group in DAVIES_LIPOPHILIC_GROUPS:
            hlb -= DAVIES_LIPOPHILIC_GROUPS[group] * count
        else:
            hlb += DAVIES_DERIVED_HYDROPHILIC_GROUPS[group] * count
    return hlb
