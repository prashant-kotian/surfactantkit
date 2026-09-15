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


# --- Nagarajan & Ruckenstein's full molecular-thermodynamic headgroup model
# (Langmuir 7 (1991) 2934-2969) -- Tier 3 bottleneck-resolution addition,
# 2026-09-15 (see benchmark/paper3_groundzero/BOTTLENECK_RESOLUTION_PLAN.md
# item 7). PDF provided directly by the user 2026-09-15 (this exact primary
# source, not a secondary citation) -- a real, DIFFERENT, more complete
# electrostatic/dipole headgroup treatment than the simpler 2002-paper
# headgroup_prefactor_A approach above (that shortcut folds everything into
# one pre-combined constant; this earlier, more general 1991 theory instead
# needs each headgroup's own real geometric parameters).
#
# EXTENDED same day with Nagarajan's own later, more comprehensive
# review chapter -- "Theory of Micelle Formation," Ch. 1 in Structure-
# Performance Relationships in Surfactants, 2nd ed., Taylor & Francis
# (2003), Table 3 ("Molecular Constants for Surfactant Headgroups"),
# primary PDF read in full 2026-09-15. This REAL, independent later
# source reports IDENTICAL values (converted nm->Angstrom) for every
# headgroup already in this table from the 1991 paper -- sodium sulfate
# (17/17/5.45), sodium sulfonate (17/17/3.85), N-betaine (30/21/d=5.0),
# and glucoside (40/21) -- a genuine cross-confirmation across two
# independent real sources, not just one paper's own internal
# consistency. It ALSO closes the two previously-disclosed gaps
# directly: real constants now exist for CATIONIC quaternary ammonium
# (trimethyl ammonium bromide) and CARBOXYLATE (both Na+ and K+
# counterions) headgroups, plus 4 more nonionic classes and a SECOND
# real zwitterionic class (lecithin, a phospholipid -- note lecithin
# genuinely needs BOTH delta (for its own charged phosphate) AND d (for
# its zwitterionic dipole character), a real structural complexity this
# table discloses rather than picks one representation):
#
#   headgroup class          a_p (A^2)  a_o (A^2)  delta (A)  d (A)  class
#   sodium sulfate             17.0       17.0       5.45      -     anionic
#   sodium sulfonate            17.0       17.0       3.85      -     anionic
#   sodium carboxylate          11.0       11.0       5.55      -     anionic
#   potassium carboxylate       11.0       11.0       6.00      -     anionic
#   trimethyl ammonium bromide  54.0       21.0       3.45      -     cationic
#   pyridinium bromide          34.0       21.0       2.20      -     cationic
#   N-betaine                   30.0       21.0       0.70     5.0   zwitterionic
#   lecithin                    45.0       42.0       6.50     6.2   zwitterionic
#   glucoside                   40.0       21.0        -        -    nonionic
#   methyl sulfoxide             39.0       21.0        -        -    nonionic
#   dimethyl phosphene oxide     48.0       21.0        -        -    nonionic
#   beta-maltoside               43.0       21.0        -        -    nonionic
#   N-methyl glucamine           34.0       21.0        -        -    nonionic
#
# Real, disclosed remaining scope: no other cationic headgroup classes
# (e.g. imidazolium) or gemini/dimeric-specific molecular constants are
# in either source paper.
NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS = {
    "sodium_sulfate": {"a_p_A2": 17.0, "a_o_A2": 17.0, "delta_A": 5.45, "class": "anionic"},
    "sodium_sulfonate": {"a_p_A2": 17.0, "a_o_A2": 17.0, "delta_A": 3.85, "class": "anionic"},
    "sodium_carboxylate": {"a_p_A2": 11.0, "a_o_A2": 11.0, "delta_A": 5.55, "class": "anionic"},
    "potassium_carboxylate": {"a_p_A2": 11.0, "a_o_A2": 11.0, "delta_A": 6.00, "class": "anionic"},
    "trimethyl_ammonium_bromide": {"a_p_A2": 54.0, "a_o_A2": 21.0, "delta_A": 3.45, "class": "cationic"},
    "pyridinium_bromide": {"a_p_A2": 34.0, "a_o_A2": 21.0, "delta_A": 2.20, "class": "cationic"},
    "n_betaine": {"a_p_A2": 30.0, "a_o_A2": 21.0, "delta_A": 0.70, "d_A": 5.0, "class": "zwitterionic"},
    "lecithin": {"a_p_A2": 45.0, "a_o_A2": 42.0, "delta_A": 6.50, "d_A": 6.2, "class": "zwitterionic"},
    "glucoside": {"a_p_A2": 40.0, "a_o_A2": 21.0, "class": "nonionic"},
    "methyl_sulfoxide": {"a_p_A2": 39.0, "a_o_A2": 21.0, "class": "nonionic"},
    "dimethyl_phosphene_oxide": {"a_p_A2": 48.0, "a_o_A2": 21.0, "class": "nonionic"},
    "beta_maltoside": {"a_p_A2": 43.0, "a_o_A2": 21.0, "class": "nonionic"},
    "n_methyl_glucamine": {"a_p_A2": 34.0, "a_o_A2": 21.0, "class": "nonionic"},
}


def nagarajan_ruckenstein_steric_headgroup_free_energy(area_per_molecule_A2: float, a_p_A2: float) -> float:
    """Steric repulsion free energy (dimensionless, units of kT) between
    headgroups crowded at the aggregate surface -- Nagarajan &
    Ruckenstein, Langmuir 7 (1991) 2934-2969, eq. 66 (van der Waals hard-
    particle approximation):

        (delta_mu)_steric / kT = -ln(1 - a_p/a)

    area_per_molecule_A2: the actual area per molecule at the aggregate
    surface (a in the source paper -- this is the same quantity
    critical_packing_parameter's head_area represents, or a candidate
    value being evaluated in a free-energy search). a_p_A2: the
    headgroup's own real cross-sectional area (see
    NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS), NEVER guessed. Diverges
    (correctly, not a bug) as area_per_molecule_A2 approaches a_p_A2 from
    above -- headgroups physically cannot pack tighter than their own
    hard-core area; raises instead of returning a nonsensical value if
    area_per_molecule_A2 <= a_p_A2.
    """
    if a_p_A2 <= 0:
        raise ValueError("a_p_A2 must be positive")
    if area_per_molecule_A2 <= a_p_A2:
        raise ValueError(
            f"area_per_molecule_A2={area_per_molecule_A2} must exceed the headgroup's own hard-core "
            f"area a_p_A2={a_p_A2} -- headgroups cannot pack tighter than their own cross-sectional area"
        )
    return -math.log(1.0 - a_p_A2 / area_per_molecule_A2)


def nagarajan_ruckenstein_dipole_headgroup_free_energy(
    area_per_molecule_A2: float,
    radius_A: float,
    d_A: float,
    geometry: str = "sphere",
    temperature_K: float = 298.15,
    dielectric_constant: float = 80.0,
) -> float:
    """Dipole-dipole interaction free energy (dimensionless, kT units)
    between ZWITTERIONIC headgroups (e.g. N-betaine) packed at the
    aggregate surface -- Nagarajan & Ruckenstein 1991, eqs. 67-68:

        sphere/globular/endcap: (delta_mu)_dipole/kT = 2*pi*e^2*R/(eps*a*kT) * [d/(d+R)]
        cylinder (middle part): (delta_mu)_dipole/kT = 2*pi*e^2*R/(eps*a*kT) * ln(1 + d/R)

    All lengths converted explicitly to cm and areas to cm^2 before
    combining with e (in esu) -- the same explicit-CGS-conversion
    pattern already used by nagarajan_debye_huckel_kappa_inverse in this
    module, deliberately NOT a rescaled-charge shortcut (a real, similar
    unit bug -- Faraday's constant vs. bare elementary charge -- was
    found and fixed elsewhere in this project's own electrostatics work
    once already; this function is written to make the conversion
    auditable line by line instead of trusting algebra done once).

    area_per_molecule_A2: the actual area per molecule at the aggregate
    surface (a in the source paper). radius_A: sphere/globular/endcap
    radius, or cylinder radius, per geometry. d_A: real charge-
    separation distance for the zwitterionic headgroup (e.g. d_A=5.0 for
    N-betaine -- see NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS, never
    guessed). geometry: 'sphere' (or 'globular'/'endcap') or 'cylinder'.
    """
    if area_per_molecule_A2 <= 0:
        raise ValueError("area_per_molecule_A2 must be positive")
    if radius_A <= 0:
        raise ValueError("radius_A must be positive")
    if d_A <= 0:
        raise ValueError("d_A must be positive")

    area_cm2 = area_per_molecule_A2 * _CM_PER_ANGSTROM ** 2
    radius_cm = radius_A * _CM_PER_ANGSTROM
    kT_erg = _K_B_ERG * temperature_K
    prefactor = 2.0 * math.pi * (_E_ESU ** 2) * radius_cm / (dielectric_constant * area_cm2 * kT_erg)

    key = geometry.lower()
    if key in ("sphere", "globular", "endcap"):
        bracket = d_A / (d_A + radius_A)  # a ratio of two lengths in the same unit -- unit-independent
    elif key == "cylinder":
        bracket = math.log(1.0 + d_A / radius_A)
    else:
        raise ValueError("geometry must be 'sphere' (or 'globular'/'endcap') or 'cylinder'")
    return prefactor * bracket


def nagarajan_ruckenstein_ionic_headgroup_free_energy(
    area_per_molecule_A2: float,
    core_radius_A: float,
    delta_A: float,
    counterion_concentration_M: float,
    added_salt_M: float = 0.0,
    temperature_K: float = 298.15,
    dielectric_constant: float = 80.0,
    geometry: str = "sphere",
) -> float:
    """Ionic (charged-headgroup) interaction free energy (dimensionless,
    kT units) at the aggregate surface -- Nagarajan & Ruckenstein 1991,
    eqs. 70-73, a curvature-corrected analytical solution to the Poisson-
    Boltzmann equation (a real, more complete alternative to the plain
    Debye-Huckel approximation with an empirical 0.46 correction factor
    the SAME paper's own earlier work used, and to the simpler 2002-
    paper headgroup_prefactor_A shortcut in this module).

        s = 4*pi*e^2 / (eps*kappa*a_delta*kT)
        (delta_mu)_ionic/kT = 2*[ln(s/2 + sqrt(1+(s/2)^2))
                                 - (2/s)*(sqrt(1+(s/2)^2) - 1)
                                 - (2*C/(kappa*s))*ln(0.5*(1+sqrt(1+(s/2)^2)))]
        C = 2/(R+delta)   [sphere/globular/endcap]
          = 1/(R+delta)   [cylinder, middle part]
        kappa = sqrt(8*pi*n0*e^2 / (eps*kT)),  n0 = (C1+Cadd)*N_A/1000

    area_per_molecule_A2: a_delta in the source paper -- the area per
    molecule evaluated AT distance delta from the hydrophobic core
    surface (a real candidate value in a free-energy search, or an
    already-known equilibrium area; not the same as a_p, the headgroup's
    own bare cross-sectional area). core_radius_A: R in the source paper
    -- sphere/globular/endcap radius, or cylinder radius, per geometry.
    delta_A: real, headgroup-specific distance from the hydrophobic core
    surface at which ionic interactions are evaluated (see
    NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS -- e.g. 5.45 A for sodium
    sulfate, 3.85 A for sodium sulfonate -- never guessed).
    counterion_concentration_M: C1 in the source paper, the molar
    concentration of singly-dispersed surfactant (its own counterions);
    added_salt_M: Cadd, any additional electrolyte -- both real,
    system-specific inputs, never guessed (same discipline as
    nagarajan_debye_huckel_kappa_inverse's own cmc_M-as-ionic-strength
    approximation, generalized here to accept added salt explicitly).
    geometry: 'sphere' (or 'globular'/'endcap') or 'cylinder'.
    """
    if area_per_molecule_A2 <= 0:
        raise ValueError("area_per_molecule_A2 must be positive")
    if core_radius_A <= 0:
        raise ValueError("core_radius_A must be positive")
    if delta_A <= 0:
        raise ValueError("delta_A must be positive")
    if counterion_concentration_M < 0 or added_salt_M < 0:
        raise ValueError("counterion_concentration_M and added_salt_M must be non-negative")
    if counterion_concentration_M + added_salt_M <= 0:
        raise ValueError("counterion_concentration_M + added_salt_M must be positive (kappa is undefined at zero ionic strength)")

    n0_per_cm3 = (counterion_concentration_M + added_salt_M) * _N_AVOGADRO / 1000.0
    kappa_per_cm = math.sqrt(8.0 * math.pi * n0_per_cm3 * (_E_ESU ** 2) / (dielectric_constant * _K_B_ERG * temperature_K))
    kappa_per_A = kappa_per_cm * _CM_PER_ANGSTROM  # same conversion as nagarajan_debye_huckel_kappa_inverse

    # s = 4*pi*e^2/(eps*kappa*a_delta*kT) -- explicit cm conversion (area
    # back to cm^2, reusing kappa_per_cm already computed above), same
    # auditable pattern as the dipole function above, not a rescaled-
    # charge shortcut.
    area_cm2 = area_per_molecule_A2 * _CM_PER_ANGSTROM ** 2
    kT_erg = _K_B_ERG * temperature_K
    s = 4.0 * math.pi * (_E_ESU ** 2) / (dielectric_constant * kappa_per_cm * area_cm2 * kT_erg)

    key = geometry.lower()
    if key in ("sphere", "globular", "endcap"):
        C = 2.0 / (core_radius_A + delta_A)
    elif key == "cylinder":
        C = 1.0 / (core_radius_A + delta_A)
    else:
        raise ValueError("geometry must be 'sphere' (or 'globular'/'endcap') or 'cylinder'")

    sqrt_term = math.sqrt(1.0 + (s / 2.0) ** 2)
    bracket = (
        math.log(s / 2.0 + sqrt_term)
        - (2.0 / s) * (sqrt_term - 1.0)
        - (2.0 * C / (kappa_per_A * s)) * math.log(0.5 * (1.0 + sqrt_term))
    )
    return 2.0 * bracket


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
