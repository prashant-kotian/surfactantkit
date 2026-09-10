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
