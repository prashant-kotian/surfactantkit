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

import math

_V_CH3 = 54.3   # cubic Angstrom, 25 C
_V_CH2 = 26.9   # cubic Angstrom, 25 C

# CGS-Gaussian constants for Nagarajan's electrostatic headgroup model below
# (the source paper works entirely in CGS-Gaussian units -- esu, cm, erg --
# not SI; kept local/self-contained to this one alternative-method function
# rather than reconciling unit systems with the rest of this SI-based module).
_E_ESU = 4.803204e-10       # elementary charge, esu (CGS-Gaussian)
_K_B_ERG = 1.380649e-16     # Boltzmann constant, erg/K
_N_AVOGADRO = 6.02214076e23  # /mol
_CM_PER_ANGSTROM = 1.0e-8


def tanford_tail_volume(n_carbons: int) -> float:
    """Hydrophobic tail volume (cubic Angstrom) for a saturated,
    unbranched alkyl chain of n_carbons carbons. Equivalent to the
    shorthand 27.4 + 26.9*n_carbons.

    CONFIRMED 2026-09-11 directly against a freely-hosted full primary
    source (Bales, Messina, Vidal, Peric & Nascimento, J. Phys. Chem. B
    1998, 102, 10347-10358, their eq. 9: Vtail=27.4+26.9*Nc, Nc=12 taken
    as the number of methylene-equivalent carbons for SDS -- exactly
    this project's own already-cited value and formula, previously
    confirmed only via a citing secondary paper's Table 2, now against
    the primary text itself, hosted free by the author at
    csun.edu/~vcphy00s/). The same paper's own illustrative example
    (SDS, 69 mM, salt-free, referencing Cabane's SANS data) uses
    NA=63 with this exact Vtail=350.2 A^3 -- see
    tests/test_hlb_cpp.py for both cross-checks.
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
    if tail_volume_A3 <= 0:
        raise ValueError("tail_volume_A3 must be positive")
    if core_radius_A <= 0:
        raise ValueError("core_radius_A must be positive")
    core_volume = (4.0 / 3.0) * math.pi * core_radius_A ** 3
    return core_volume / tail_volume_A3


def nagarajan_debye_huckel_kappa_inverse(cmc_M: float, temperature_K: float = 298.15, dielectric_constant: float = 80.0) -> float:
    """Debye-Huckel inverse screening length kappa^-1 (Angstrom) at the
    CMC of an IONIC surfactant with no added electrolyte -- Nagarajan,
    Langmuir 18 (2002) 31-38, eq. 7. Real, sourced alternative
    (headgroup-electrostatics-aware) input to the equilibrium area per
    molecule (see nagarajan_equilibrium_area_ionic), used to show that
    the packing parameter genuinely varies with tail length for ionic
    surfactants (contrary to the common simplifying assumption that only
    the headgroup, not the tail, sets the area per molecule).

    Uses the source paper's own explicit approximation: "near the cmc,
    the ionic strength is determined primarily by the concentration of
    the singly dispersed surfactant in solution, which is practically
    equivalent to the cmc" -- i.e. ionic strength ~= cmc_M for a 1:1
    ionic surfactant with NO added salt. Not valid with added
    electrolyte (use the real total ionic strength via
    electrostatics.debye_length instead in that case).

    Verified numerically against the source paper's own Table 2 (real
    worked example, sodium alkyl sulfates): for n_C=12, cmc=0.008 M,
    this reproduces kappa^-1 = 34.43 Angstrom to within ~0.2% (the small
    residual is rounding in the paper's own tabulated inputs, not a
    formula error) -- see tests/test_cpp.py.
    """
    if cmc_M <= 0:
        raise ValueError("cmc_M must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if dielectric_constant <= 0:
        raise ValueError("dielectric_constant must be positive")
    n0_per_cm3 = cmc_M * _N_AVOGADRO / 1000.0  # mol/L -> mol/cm^3 -> ions/cm^3
    kappa_per_cm = math.sqrt(
        8.0 * math.pi * n0_per_cm3 * (_E_ESU ** 2) / (dielectric_constant * _K_B_ERG * temperature_K)
    )
    kappa_per_angstrom = kappa_per_cm * _CM_PER_ANGSTROM
    return 1.0 / kappa_per_angstrom


def nagarajan_equilibrium_area_ionic(cmc_M: float, tail_length_A: float, headgroup_prefactor_A: float, temperature_K: float = 298.15, dielectric_constant: float = 80.0) -> float:
    """Equilibrium area per molecule a_e (square Angstrom) for an IONIC
    surfactant micelle, accounting for chain-length-dependent
    electrostatic headgroup repulsion -- Nagarajan, Langmuir 18 (2002)
    31-38, eqs. 6-8, combined with nagarajan_debye_huckel_kappa_inverse.

    a_e = headgroup_prefactor_A / sqrt(1 + tail_length_A/kappa_inverse_A)

    This is a real, sourced ALTERNATIVE to simply plugging an
    experimentally-measured or assumed-constant a0 into
    critical_packing_parameter: the paper's central finding is that a_e
    (and hence v0/(a_e*l0), the packing parameter itself) is NOT
    independent of tail length for ionic surfactants -- it depends on
    tail length through the Debye screening length kappa^-1, which
    itself depends on the cmc, which depends on tail length. This
    directly contradicts the common assumption (explicit in the
    literature this paper critiques) that a_e is fixed by the headgroup
    alone. Use the returned a_e as the head_area input to
    critical_packing_parameter for an ionic surfactant when this
    tail-length effect matters (e.g. comparing a homologous series).

    headgroup_prefactor_A: sqrt(2*pi*e^2*d / (epsilon*sigma)) in
    Angstrom, where d is the headgroup's effective charge-separation
    (capacitor) distance and sigma is the bare hydrocarbon-water
    interfacial free energy per unit area -- BOTH are real, system/
    headgroup-specific physical quantities the source paper does not
    give generic universal values for (it reports only their COMBINED
    value used in its own illustrative calculation for sodium alkyl
    sulfate-like headgroups: 82 Angstrom -- pass 82.0 to reproduce that
    specific illustration, but do not treat 82 as a universal constant
    for other headgroups; source it from d and sigma for the actual
    system, or do not use this function if unknown -- do not guess).

    Verified numerically against the source paper's own Table 2: with
    headgroup_prefactor_A=82.0 (the paper's own stated value for this
    table), reproduces a_e = 67.4 Angstrom^2 at n_C=12 (cmc=0.008 M,
    tail_length_A=16.5) to within the table's own reported precision --
    see tests/test_cpp.py.
    """
    if tail_length_A <= 0:
        raise ValueError("tail_length_A must be positive")
    if headgroup_prefactor_A <= 0:
        raise ValueError("headgroup_prefactor_A must be positive")
    kappa_inverse_A = nagarajan_debye_huckel_kappa_inverse(cmc_M, temperature_K, dielectric_constant)
    return headgroup_prefactor_A / math.sqrt(1.0 + tail_length_A / kappa_inverse_A)


def estimate_axial_ratio_from_cpp_geometry(cpp: float, aggregation_number: float, n_carbons: int) -> float:
    """Real, buildable estimate of the prolate-ellipsoid axial ratio
    (a/b) for a cylindrical/rodlike micelle, closing the [COMPUTE] path
    identified for dynamics.py's Perrin axial_ratio bottleneck (see
    benchmark/paper3_groundzero/BOTTLENECK_RESOLUTION_PLAN.md item 8) --
    dynamics.hydrodynamic_radius_perrin_corrected's own docstring already
    names "the critical packing parameter's predicted morphology" as one
    legitimate source of axial_ratio; this function makes that concrete
    and quantitative instead of leaving it as an unstated intention.

    Geometric model: the hydrophobic core of a growing rodlike micelle is
    treated as a prolate ellipsoid of revolution whose MINOR semi-axis b
    is pinned at the extended tail length lc (tanford_critical_length) --
    the same packing constraint that already caps a spherical micelle's
    core radius at lc, since the hydrophobic core can never exceed the
    fully-extended chain length in its shortest cross-sectional
    dimension. The MAJOR semi-axis a then follows from real volume
    conservation given a REAL (independently measured -- fluorescence
    quenching, SLS, etc., never estimated from the spherical-micelle
    formula, which would assume the very geometry this function exists
    to correct for) aggregation number:

        V_core = aggregation_number * tanford_tail_volume(n_carbons)
        V_core = (4/3) * pi * a * b^2   (prolate ellipsoid of revolution)
        axial_ratio = a / b

    Deliberately scoped to CPP in (1/3, 1/2] -- the classical sphere-to-
    cylinder growth regime this ellipsoid-of-revolution model actually
    describes (classify_aggregate_morphology's "cylindrical/rodlike
    micelle" bucket). CPP <= 1/3 returns axial_ratio=1.0 directly (a
    sphere needs no shape correction -- Perrin's own factor is exactly 1
    there). CPP > 1/2 (vesicle/bilayer or inverted) raises rather than
    guessing: a solid ellipsoid of revolution is the wrong shape model
    entirely for a hollow bilayer/vesicle shell, and extending this
    formula there would silently apply a model outside its own domain.

    Also raises if the resulting axial_ratio < 1.0 -- a real, honest
    inconsistency signal, not silently rounded to 1: it means the
    supplied aggregation_number is smaller than what the pinned-b
    geometry needs for even a marginal sphere-to-rod transition,
    contradicting the cylindrical/rodlike classification the caller's
    own cpp value implies. Report the disagreement rather than picking
    one input over the other.
    """
    if cpp <= 0:
        raise ValueError("cpp must be positive")
    if aggregation_number <= 0:
        raise ValueError("aggregation_number must be positive")
    if n_carbons < 1:
        raise ValueError("n_carbons must be at least 1")
    if cpp <= 1.0 / 3.0:
        return 1.0
    if cpp > 0.5:
        raise ValueError(
            f"cpp={cpp:.4f} is above the cylindrical/rodlike regime (0.333, 0.5] this "
            "prolate-ellipsoid-of-revolution model describes -- a solid ellipsoid is not "
            "the right shape for a vesicle/bilayer or inverted-structure core; do not use "
            "this function outside that CPP range"
        )
    b = tanford_critical_length(n_carbons)
    v_tail = tanford_tail_volume(n_carbons)
    v_core = aggregation_number * v_tail
    a = (3.0 * v_core) / (4.0 * math.pi * b ** 2)
    axial_ratio = a / b
    if axial_ratio < 1.0:
        if math.isclose(axial_ratio, 1.0, rel_tol=1e-9):
            # Real floating-point boundary case, not a genuine physical
            # inconsistency: an aggregation_number that lands EXACTLY at
            # the sphere/cylinder boundary (a==b) can round to axial_ratio
            # a hair below 1.0 through this division, purely from binary
            # floating-point representation, not because the inputs
            # actually disagree -- clamp to exactly 1.0 rather than raise.
            return 1.0
        raise ValueError(
            f"aggregation_number={aggregation_number:.1f} implies a core major semi-axis "
            f"({a:.1f} A) smaller than the pinned minor axis (b={b:.1f} A) -- this "
            f"aggregation number is more consistent with a spherical micelle than the "
            f"cylindrical/rodlike morphology implied by cpp={cpp:.4f}; the two inputs "
            "disagree, do not trust either without checking both"
        )
    return axial_ratio


def classify_aggregate_morphology(cpp: float) -> str:
    """Expected aggregate morphology from the critical packing parameter.

    Thresholds: <=1/3 spherical micelles; (1/3, 1/2] cylindrical/rod
    micelles; (1/2, 1] vesicles/bilayers; >1 inverted structures --
    Israelachvili's classical thresholds.

    Real-world caveat, checked against the primary source (2026-09-11,
    PDF provided by the user): Schafer, Kolli, Christensen, Bore,
    Diezemann, Gauss, Milano, Lund & Cascella, Angew. Chem. Int. Ed. 59
    (2020) 18591-18598 (SAXS/SANS + simulation study of SDS's real
    sphere-to-cylinder transition). Their own reported number is P=0.493
    for the effective packing parameter of the paired ("dimeric") SDS
    subunit that drives the transition -- close to the classical 1/2
    threshold used here for the cylinder/vesicle boundary, NOT a
    departure from it, and NOT a bulk-CPP value near 0.2 (an earlier,
    inaccurate secondary-summary paraphrase of this paper, corrected
    this pass now that the primary text is in hand). Their own real
    finding is more subtle than a shifted threshold: only ~10-20% of SDS
    monomers need to pair into these locally-P=0.493 dimers to drive the
    sphere-to-cylinder transition, "without any discontinuity in the
    overall packing parameter" (their own words) -- i.e. the BULK/
    monomeric CPP computed by this function's own formula does not
    itself jump at the transition; the local packing of a transient
    substructure does. Treat this function's classical thresholds as
    what they are (a monomer-level rule of thumb), not as a claim that
    real morphology transitions are driven by a single sharp bulk-CPP
    value crossing 1/3 or 1/2.
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
