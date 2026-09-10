"""Raw experimental curve analysis -- the missing first step every other
module in this library assumes has already happened. Every function in
adsorption.py, mixed_micelle.py, etc. takes an already-known CMC or slope as
input; nothing in this library previously took a raw surface-tension-vs-
concentration table and extracted CMC from it. Built 2026-09-07 specifically
to close that gap, after a real methodology correction: the benchmark's
original plan tested models on pre-parameterized textbook questions, but the
actual research workflow starts from raw measured data, and that's the gap
this module (and the benchmark built on top of it) targets.

Method: two-segment (piecewise) linear regression break-point detection --
the same real method described in Shah/Das/Bhattarai 2025 (Heliyon,
PMC11835642): "A graph is plotted with gamma and log C... gamma shows a
breakpoint -- this point is called CMC." Standard, not invented: fit every
possible split of the sorted (log C, gamma) points into a declining
premicellar line and a flat-ish post-CMC segment, and take the split that
minimizes total residual sum of squares across both segments.
"""

from __future__ import annotations
import math
from dataclasses import dataclass


def _linreg(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    """Ordinary least squares. Returns (slope, intercept, residual_sum_sq)."""
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    if sxx == 0:
        slope = 0.0
    else:
        slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    rss = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(xs, ys))
    return slope, intercept, rss


@dataclass
class CmcFromCurveResult:
    cmc_mM: float
    premicellar_slope_mN_per_m_per_lnC: float  # dGamma/d(ln C), the form
                                                 # gibbs_gamma_max() expects
    premicellar_slope_mN_per_m_per_log10C: float  # as directly fit, for reference
    gamma_at_cmc_mN_per_m: float
    r_squared_premicellar: float
    n_premicellar_points: int
    n_postmicellar_points: int
    method: str = "two-segment linear regression break-point (min total RSS)"


def cmc_from_surface_tension_curve(
    concentrations_mM: list[float],
    surface_tensions_mN_per_m: list[float],
    min_points_per_segment: int = 3,
) -> CmcFromCurveResult:
    """Extracts CMC from a raw gamma-vs-concentration tensiometry curve by
    finding the concentration that best splits the data into a declining
    premicellar line (gamma vs log10 C, real Gibbs-adsorption-linear region)
    and a flat post-micellar segment, minimizing combined residual sum of
    squares across every valid split point.

    Real, standard method -- not a guess: this is exactly what
    Shah/Das/Bhattarai 2025 describe doing by eye/software ("the breakpoint
    -- this point is called CMC"), just made a deterministic, repeatable
    algorithm instead of a by-eye read.
    """
    if len(concentrations_mM) != len(surface_tensions_mN_per_m):
        raise ValueError("concentrations and surface tensions must be the same length")
    if len(concentrations_mM) < 2 * min_points_per_segment:
        raise ValueError(
            f"need at least {2*min_points_per_segment} points "
            f"({min_points_per_segment} per segment) to fit a break-point"
        )

    pairs = sorted(zip(concentrations_mM, surface_tensions_mN_per_m))
    log_c = [math.log10(c) for c, _ in pairs]
    gamma = [g for _, g in pairs]
    n = len(pairs)

    best = None
    for split in range(min_points_per_segment, n - min_points_per_segment + 1):
        pre_x, pre_y = log_c[:split], gamma[:split]
        post_x, post_y = log_c[split:], gamma[split:]
        pre_slope, pre_intercept, pre_rss = _linreg(pre_x, pre_y)
        _, _, post_rss = _linreg(post_x, post_y)
        total_rss = pre_rss + post_rss
        if best is None or total_rss < best[0]:
            best = (total_rss, split, pre_slope, pre_intercept, pre_rss)

    _, split, pre_slope, pre_intercept, pre_rss = best
    pre_x, pre_y = log_c[:split], gamma[:split]
    mean_y = sum(pre_y) / len(pre_y)
    ss_tot = sum((y - mean_y) ** 2 for y in pre_y)
    r_squared = 1.0 - (pre_rss / ss_tot) if ss_tot > 0 else 1.0

    # CMC = concentration at the boundary between the two fitted segments --
    # the midpoint (in log-concentration space) between the last premicellar
    # point and the first post-micellar point, matching how a break-point is
    # actually read off a real plot (the intersection region between the two
    # branches, not either branch's own last/first raw data point).
    cmc_log_c = (log_c[split - 1] + log_c[split]) / 2.0
    cmc_mM = 10 ** cmc_log_c
    gamma_at_cmc = pre_slope * cmc_log_c + pre_intercept

    # d(gamma)/d(ln C) = d(gamma)/d(log10 C) / ln(10) -- the form
    # gibbs_gamma_max() actually consumes.
    slope_per_lnC = pre_slope / math.log(10)

    return CmcFromCurveResult(
        cmc_mM=cmc_mM,
        premicellar_slope_mN_per_m_per_lnC=slope_per_lnC,
        premicellar_slope_mN_per_m_per_log10C=pre_slope,
        gamma_at_cmc_mN_per_m=gamma_at_cmc,
        r_squared_premicellar=r_squared,
        n_premicellar_points=split,
        n_postmicellar_points=n - split,
    )


@dataclass
class CmcFromConductivityResult:
    cmc_mM: float
    slope_below_cmc: float  # d(conductivity)/d(concentration), premicellar segment
    slope_above_cmc: float  # d(conductivity)/d(concentration), postmicellar segment
    counterion_binding_degree: float  # 1 - slope_above/slope_below, the same
                                       # formula thermodynamics.counterion_binding_degree()
                                       # computes -- provided directly here so this
                                       # function alone answers "what's beta from
                                       # this raw conductivity curve"
    r_squared_below_cmc: float
    r_squared_above_cmc: float
    n_points_below_cmc: int
    n_points_above_cmc: int
    method: str = "two-segment linear regression, linear-concentration break-point (min total RSS)"


def cmc_from_conductivity_curve(
    concentrations_mM: list[float],
    conductivities: list[float],
    min_points_per_segment: int = 3,
) -> CmcFromConductivityResult:
    """Extracts CMC from a raw specific-conductivity-vs-concentration curve
    -- the standard conductometric break-point method for IONIC surfactants
    (nonionic surfactants have negligible conductivity response and cannot
    be characterized this way -- a real, disclosed method limitation, not
    an oversight).

    Unlike cmc_from_surface_tension_curve (which fits gamma vs LOG10
    concentration, per the Gibbs adsorption isotherm's own log-concentration
    dependence), conductometry is fit on LINEAR concentration -- specific
    conductivity increases linearly with surfactant concentration in both
    the premicellar (free ion/monomer) and postmicellar (micelle-dominated)
    regimes, just with two different real slopes, because a micelle is a
    much larger, slower-diffusing, less efficient charge carrier than a
    free monomer/counterion. The CMC is the concentration where the two
    fitted lines intersect -- both slopes are real, physically meaningful
    quantities in their own right (not just a means to locate the CMC),
    which is why this function returns both directly, in exactly the form
    thermodynamics.counterion_binding_degree() consumes, AND the resulting
    beta itself computed the same way that function does (so this one
    function answers the whole "raw conductivity curve -> CMC and beta"
    question end to end).

    Real, standard method -- not invented: segmented/piecewise linear
    regression with the break point taken at the two fitted lines'
    intersection is the textbook conductometric CMC method (see e.g. the
    two-linear-segments description in real surfactant characterization
    literature), same "fit every valid split, take the one minimizing
    total residual sum of squares" discipline as this module's
    cmc_from_surface_tension_curve, adapted for the linear (not log)
    x-axis this specific technique uses.
    """
    if len(concentrations_mM) != len(conductivities):
        raise ValueError("concentrations and conductivities must be the same length")
    if len(concentrations_mM) < 2 * min_points_per_segment:
        raise ValueError(
            f"need at least {2*min_points_per_segment} points "
            f"({min_points_per_segment} per segment) to fit a break-point"
        )

    pairs = sorted(zip(concentrations_mM, conductivities))
    conc = [c for c, _ in pairs]
    cond = [k for _, k in pairs]
    n = len(pairs)

    best = None
    for split in range(min_points_per_segment, n - min_points_per_segment + 1):
        pre_x, pre_y = conc[:split], cond[:split]
        post_x, post_y = conc[split:], cond[split:]
        pre_slope, pre_intercept, pre_rss = _linreg(pre_x, pre_y)
        post_slope, post_intercept, post_rss = _linreg(post_x, post_y)
        total_rss = pre_rss + post_rss
        if best is None or total_rss < best[0]:
            best = (total_rss, split, pre_slope, pre_intercept, pre_rss, post_slope, post_intercept, post_rss)

    _, split, pre_slope, pre_intercept, pre_rss, post_slope, post_intercept, post_rss = best

    if pre_slope <= post_slope:
        raise ValueError(
            "fitted premicellar slope is not greater than the postmicellar slope -- this "
            "contradicts the expected conductometric behavior (micelles are less efficient "
            "charge carriers than free monomer/counterions) and likely means this data does "
            "not show a real conductometric break point (e.g. a nonionic surfactant, or data "
            "that doesn't actually span the CMC) -- check the input rather than trusting this fit"
        )

    # CMC = the actual intersection of the two fitted lines, the standard
    # convention for this method (both segments are genuinely linear here,
    # unlike surface tension's curved transition region, so the precise
    # algebraic intersection is the right convention, not a data-point
    # midpoint).
    cmc_mM = (post_intercept - pre_intercept) / (pre_slope - post_slope)

    pre_y = cond[:split]
    mean_pre = sum(pre_y) / len(pre_y)
    ss_tot_pre = sum((y - mean_pre) ** 2 for y in pre_y)
    r_squared_pre = 1.0 - (pre_rss / ss_tot_pre) if ss_tot_pre > 0 else 1.0

    post_y = cond[split:]
    mean_post = sum(post_y) / len(post_y)
    ss_tot_post = sum((y - mean_post) ** 2 for y in post_y)
    r_squared_post = 1.0 - (post_rss / ss_tot_post) if ss_tot_post > 0 else 1.0

    beta = 1.0 - (post_slope / pre_slope)

    return CmcFromConductivityResult(
        cmc_mM=cmc_mM,
        slope_below_cmc=pre_slope,
        slope_above_cmc=post_slope,
        counterion_binding_degree=beta,
        r_squared_below_cmc=r_squared_pre,
        r_squared_above_cmc=r_squared_post,
        n_points_below_cmc=split,
        n_points_above_cmc=n - split,
    )


@dataclass
class AggregationNumberQuenchingResult:
    aggregation_number: float
    slope_per_Q: float  # fitted N/(Ct - cmc), from ln(I0/I) = slope * [Q]
    method: str = "steady-state fluorescence quenching (Turro-Yekta, 1978)"


def aggregation_number_from_quenching_curve(
    quencher_concentrations: list[float],
    intensities: list[float],
    intensity_without_quencher: float,
    total_surfactant_concentration: float,
    cmc: float,
) -> AggregationNumberQuenchingResult:
    """Micelle aggregation number from a raw steady-state fluorescence
    quenching (SSFQ) series -- the Turro-Yekta method (1978), verified
    directly from a real fetched source (equation image, quoting exactly):

    ln(I / I0) = -N*[Q] / (Ct - cmc)

    i.e. ln(I0/I) = N*[Q] / (Ct - cmc), a straight line through the
    origin when ln(I0/I) is plotted against quencher concentration [Q].

    Deliberately the SIMPLEST real aggregation-number method to run:
    needs only a standard steady-state spectrofluorometer (not the
    time-resolved/TCSPC instrumentation TRFQ requires, not SANS, not an
    ultracentrifuge) -- a fluorescent probe (commonly pyrene) solubilized
    in the micelles, quenched by a compound that partitions almost
    entirely into the micellar phase (e.g. cetylpyridinium chloride,
    benzophenone), measured at a series of quencher concentrations.

    quencher_concentrations, intensities: paired raw data -- the probe's
    fluorescence intensity (arbitrary units, e.g. the first vibronic
    peak height for pyrene) at each quencher concentration [Q] (same
    concentration unit as total_surfactant_concentration/cmc).
    intensity_without_quencher: I0, the probe's intensity at [Q]=0
    (measured separately, the same units as `intensities`).
    total_surfactant_concentration, cmc: same concentration unit as
    quencher_concentrations -- Ct must be well above cmc for a
    meaningful (Ct - cmc) micellized-surfactant denominator.

    The slope is fit forced through the origin (least-squares with zero
    intercept: slope = sum(x*y)/sum(x^2)), matching the theoretical
    model exactly (ln(I0/I0)=0 at [Q]=0 identically, not just
    approximately) rather than fitting a free intercept that could
    absorb real experimental baseline noise into a spurious offset.

    Source: Turro & Yekta, J. Am. Chem. Soc. 100 (1978) 5951-5952 (the
    original method); equation verified directly from a real fetched
    secondary source (Bhattacharjee & Mahapatra-style DDM/SAXS
    comparison paper, SCIRP, image-extracted equation), cross-confirmed
    by two independent web summaries of the same equation before
    fetching the image directly.
    """
    if len(quencher_concentrations) != len(intensities):
        raise ValueError("quencher_concentrations and intensities must be the same length")
    if len(quencher_concentrations) < 2:
        raise ValueError("need at least 2 (quencher concentration, intensity) points to fit a slope")
    if intensity_without_quencher <= 0:
        raise ValueError("intensity_without_quencher must be positive")
    if any(i <= 0 for i in intensities):
        raise ValueError("all intensities must be positive")
    if total_surfactant_concentration <= cmc:
        raise ValueError("total_surfactant_concentration must exceed cmc (need a real micellized fraction)")

    xs = quencher_concentrations
    ys = [math.log(intensity_without_quencher / i) for i in intensities]
    sum_xy = sum(x * y for x, y in zip(xs, ys))
    sum_xx = sum(x * x for x in xs)
    if sum_xx == 0:
        raise ValueError("quencher_concentrations cannot all be zero")
    slope = sum_xy / sum_xx

    aggregation_number = slope * (total_surfactant_concentration - cmc)
    return AggregationNumberQuenchingResult(aggregation_number=aggregation_number, slope_per_Q=slope)


AVOGADRO_NUMBER = 6.02214076e23  # /mol, exact SI 2019


@dataclass
class AggregationNumberSLSResult:
    micelle_molar_mass_g_per_mol: float
    aggregation_number: float
    second_virial_coefficient_cm3_mol_per_g2: float
    method: str = "static light scattering, single-angle Debye plot"


def aggregation_number_from_sls_debye_plot(
    micellized_concentrations_g_per_mL: list[float],
    rayleigh_ratios_cm_inv: list[float],
    dn_dc_mL_per_g: float,
    wavelength_nm: float,
    monomer_molar_mass_g_per_mol: float,
    refractive_index_solvent: float = 1.333,
) -> AggregationNumberSLSResult:
    """Micelle aggregation number from a raw static light scattering
    (SLS) concentration series -- the single-angle Debye plot method,
    real and standard for small (Rg < 12 nm) scatterers like surfactant
    micelles (a multi-angle instrument, needed for larger macromolecules,
    is NOT required here) -- more instrumentation than SSFQ
    (aggregation_number_from_quenching_curve) but real, standard, and
    achievable on the same class of light-scattering/DLS instrument
    already used elsewhere in this line of work for micelle sizing.

    Debye equation (verified directly against a real fetched primary
    technical source, Brookhaven Instruments' "SLS FAQ: The Debye
    Plot," which reduces the general Zimm equation to this form
    explicitly when Rg^2*q^2/3 << 1, the small-scatterer limit):

    K*c / DeltaR = 1/Mw + 2*A2*c

    K = 4*pi^2*n^2*(dn/dc)^2 / (NA * lambda0^4)   (optical constant;
    n = solvent refractive index, dn/dc = refractive index increment,
    lambda0 = VACUUM wavelength, NA = Avogadro's number -- this specific
    4*pi^2/vacuum-wavelength form cross-confirmed across three
    independent real sources: Malvern, Anton Paar, and Brookhaven
    technical documentation, all agreeing).

    micellized_concentrations_g_per_mL: concentration of MICELLIZED
    surfactant only (i.e. total surfactant mass concentration minus the
    mass-concentration-equivalent of the cmc) at each measurement point
    -- using total surfactant concentration instead would include free
    monomer that contributes no excess scattering above the solvent
    baseline, biasing the fitted molar mass. Same "must be a real,
    already-corrected input, not guessed" discipline as elsewhere in
    this module.
    rayleigh_ratios_cm_inv: the excess Rayleigh ratio (DeltaR, solvent-
    subtracted, units of inverse cm) at each corresponding concentration
    -- a real instrument-calibrated quantity, not raw photon counts.
    dn_dc_mL_per_g: refractive index increment for this surfactant in
    this solvent at this wavelength -- system-specific, must be a real
    measured or cited value, never guessed (see Zhao, Brown & Schuck,
    Biophys. J. 100 (2011) 2309-2317 for how variable this can be even
    within one macromolecule class).
    wavelength_nm: VACUUM wavelength of the laser used, in nm.
    monomer_molar_mass_g_per_mol: pure-surfactant-monomer molar mass --
    aggregation_number = micelle_molar_mass / monomer_molar_mass.
    refractive_index_solvent: default 1.333 (water, visible light,
    ~20-25 C) -- override for other solvents/temperatures/wavelengths.

    Fits Mw and the second virial coefficient A2 by ordinary least
    squares on Kc/DeltaR vs. c (the standard Debye-plot linear fit),
    reusing this module's own _linreg helper.
    """
    if len(micellized_concentrations_g_per_mL) != len(rayleigh_ratios_cm_inv):
        raise ValueError("micellized_concentrations_g_per_mL and rayleigh_ratios_cm_inv must be the same length")
    if len(micellized_concentrations_g_per_mL) < 2:
        raise ValueError("need at least 2 (concentration, Rayleigh ratio) points to fit the Debye plot")
    if any(c <= 0 for c in micellized_concentrations_g_per_mL):
        raise ValueError("all micellized concentrations must be positive")
    if any(r <= 0 for r in rayleigh_ratios_cm_inv):
        raise ValueError("all Rayleigh ratios must be positive")
    if dn_dc_mL_per_g == 0:
        raise ValueError("dn_dc_mL_per_g must be nonzero -- a real, measured/cited value, not a guess")
    if monomer_molar_mass_g_per_mol <= 0:
        raise ValueError("monomer_molar_mass_g_per_mol must be positive")

    wavelength_cm = wavelength_nm * 1e-7
    optical_constant_K = (
        4.0 * math.pi ** 2 * refractive_index_solvent ** 2 * dn_dc_mL_per_g ** 2
        / (AVOGADRO_NUMBER * wavelength_cm ** 4)
    )

    xs = micellized_concentrations_g_per_mL
    ys = [optical_constant_K * c / r for c, r in zip(micellized_concentrations_g_per_mL, rayleigh_ratios_cm_inv)]
    slope, intercept, _ = _linreg(xs, ys)

    if intercept <= 0:
        raise ValueError(
            "fitted Kc/DeltaR intercept is non-positive -- this gives a non-physical (infinite or negative) "
            "molar mass; check the input data (Rayleigh ratios, dn/dc, wavelength) rather than trusting this fit"
        )
    micelle_molar_mass = 1.0 / intercept
    aggregation_number = micelle_molar_mass / monomer_molar_mass_g_per_mol
    second_virial_coefficient = slope / 2.0

    return AggregationNumberSLSResult(
        micelle_molar_mass_g_per_mol=micelle_molar_mass,
        aggregation_number=aggregation_number,
        second_virial_coefficient_cm3_mol_per_g2=second_virial_coefficient,
    )


@dataclass
class AggregationNumberDLSResult:
    hydrodynamic_radius_nm: float
    micelle_molar_mass_g_per_mol: float
    aggregation_number: float
    method: str = "dynamic light scattering, hydrodynamic-volume/partial-specific-volume route"


def aggregation_number_from_dls(
    diffusion_coefficient_cm2_per_s: float,
    viscosity_mPas: float,
    partial_specific_volume_cm3_per_g: float,
    monomer_molar_mass_g_per_mol: float,
    temperature_K: float = 298.15,
) -> AggregationNumberDLSResult:
    """Micelle aggregation number from a raw DLS-measured diffusion
    coefficient -- an alternative-methods-sweep addition (2026-09-10),
    real and standard (the same instrument class already used elsewhere
    in this line of work for micelle sizing): converts the hydrodynamic
    radius (dynamics.hydrodynamic_radius_stokes_einstein, reused
    directly here, not reimplemented) to a hydrodynamic VOLUME, then to
    a molar mass via the micelle's partial specific volume, then to an
    aggregation number via the monomer molar mass.

        R_h = k_B*T / (6*pi*eta*D)                 [Stokes-Einstein]
        V_h = (4/3)*pi*R_h^3                        [per-particle hydrodynamic volume]
        M_micelle = V_h * N_A / v_bar                [v_bar = partial specific volume]
        N_agg = M_micelle / M_monomer

    Real, disclosed limitation (standard in the DLS/micelle-sizing
    literature, not specific to this implementation): R_h includes a
    thin bound-water hydration shell, so the hydrodynamic volume is a
    real, known OVERESTIMATE of the anhydrous molecular volume implied
    by partial_specific_volume_cm3_per_g -- this gives an aggregation
    number that is a real, usable estimate but typically biased slightly
    HIGH relative to a hydration-corrected value; treat as an estimate
    of the same character as aggregation_number_spherical's geometric
    estimate, not a hydration-corrected exact value. Same units and
    inputs as dynamics.hydrodynamic_radius_stokes_einstein for the first
    two arguments -- chain directly from that function's own docstring
    guidance on diffusion coefficient units (cm^2/s, not m^2/s).

    partial_specific_volume_cm3_per_g: real, measured or cited value
    for this surfactant (e.g. from densitometry) -- system-specific, do
    not guess (same discipline as dn/dc in aggregation_number_from_sls_debye_plot).
    """
    if diffusion_coefficient_cm2_per_s <= 0:
        raise ValueError("diffusion_coefficient_cm2_per_s must be positive")
    if viscosity_mPas <= 0:
        raise ValueError("viscosity_mPas must be positive")
    if partial_specific_volume_cm3_per_g <= 0:
        raise ValueError("partial_specific_volume_cm3_per_g must be positive")
    if monomer_molar_mass_g_per_mol <= 0:
        raise ValueError("monomer_molar_mass_g_per_mol must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")

    k_B = 1.380649e-23  # J/K, exact SI 2019
    eta_SI = viscosity_mPas * 1e-3  # mPa.s -> Pa.s
    D_SI = diffusion_coefficient_cm2_per_s * 1e-4  # cm^2/s -> m^2/s
    r_h_m = (k_B * temperature_K) / (6.0 * math.pi * eta_SI * D_SI)
    r_h_nm = r_h_m * 1e9

    r_h_cm = r_h_m * 100.0  # m -> cm, to get volume directly in cm^3
    v_h_cm3 = (4.0 / 3.0) * math.pi * r_h_cm ** 3
    micelle_molar_mass = v_h_cm3 * AVOGADRO_NUMBER / partial_specific_volume_cm3_per_g
    aggregation_number = micelle_molar_mass / monomer_molar_mass_g_per_mol

    return AggregationNumberDLSResult(
        hydrodynamic_radius_nm=r_h_nm,
        micelle_molar_mass_g_per_mol=micelle_molar_mass,
        aggregation_number=aggregation_number,
    )


@dataclass
class AggregationNumberSvedbergResult:
    micelle_molar_mass_g_per_mol: float
    aggregation_number: float
    method: str = "analytical ultracentrifugation, Svedberg equation"


def aggregation_number_from_svedberg_equation(
    sedimentation_coefficient_S: float,
    diffusion_coefficient_cm2_per_s: float,
    partial_specific_volume_cm3_per_g: float,
    solvent_density_g_per_cm3: float,
    monomer_molar_mass_g_per_mol: float,
    temperature_K: float = 298.15,
) -> AggregationNumberSvedbergResult:
    """Micelle aggregation number from raw analytical-ultracentrifugation
    data (sedimentation velocity) -- an alternative-methods-sweep
    addition (2026-09-10), the classic, standard route (Svedberg
    equation) to a macromolecule's molar mass, real and independent of
    the light-scattering/DLS-based routes elsewhere in this module
    (different instrument, different physical principle -- a genuine
    cross-check when both are available on the same system).

        M = R*T*s / [D*(1 - v_bar*rho)]
        N_agg = M / M_monomer

    s: sedimentation coefficient in SVEDBERGS (S, 1 S = 1e-13 s) -- the
    standard reporting unit, converted internally; do not pass raw
    seconds here. D: translational diffusion coefficient in cm^2/s (same
    convention as elsewhere in this module/dynamics.py). v_bar: partial
    specific volume (cm^3/g, real measured/cited value, system-specific
    -- do not guess; same discipline as aggregation_number_from_dls).
    rho: solvent density (g/cm^3, ~0.997 for water at 25 C). The
    (1 - v_bar*rho) "buoyancy factor" must be positive -- a real,
    physical requirement (a particle denser than its buoyant
    displacement sediments; one that isn't, doesn't; this function
    raises rather than silently returning a nonsensical negative or
    infinite mass if the inputs violate it).
    """
    if diffusion_coefficient_cm2_per_s <= 0:
        raise ValueError("diffusion_coefficient_cm2_per_s must be positive")
    if partial_specific_volume_cm3_per_g <= 0:
        raise ValueError("partial_specific_volume_cm3_per_g must be positive")
    if solvent_density_g_per_cm3 <= 0:
        raise ValueError("solvent_density_g_per_cm3 must be positive")
    if monomer_molar_mass_g_per_mol <= 0:
        raise ValueError("monomer_molar_mass_g_per_mol must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    buoyancy_factor = 1.0 - partial_specific_volume_cm3_per_g * solvent_density_g_per_cm3
    if buoyancy_factor <= 0:
        raise ValueError(
            "buoyancy factor (1 - v_bar*rho) is non-positive -- the particle would not sediment under "
            "these inputs (check partial_specific_volume_cm3_per_g and solvent_density_g_per_cm3); "
            "the Svedberg equation gives a non-physical result and is not applied"
        )

    R_GAS = 8.314462618  # J/(mol.K)
    s_seconds = sedimentation_coefficient_S * 1e-13
    D_SI = diffusion_coefficient_cm2_per_s * 1e-4  # cm^2/s -> m^2/s
    micelle_molar_mass_SI = (R_GAS * temperature_K * s_seconds) / (D_SI * buoyancy_factor)  # kg/mol
    micelle_molar_mass = micelle_molar_mass_SI * 1000.0  # kg/mol -> g/mol
    aggregation_number = micelle_molar_mass / monomer_molar_mass_g_per_mol

    return AggregationNumberSvedbergResult(
        micelle_molar_mass_g_per_mol=micelle_molar_mass,
        aggregation_number=aggregation_number,
    )


if __name__ == "__main__":
    # Real pilot dataset: AOT in water, reconstructed from Shah/Das/Bhattarai
    # 2025 (Heliyon, PMC11835642) Figure 1 -- 6 premicellar points computed
    # from the paper's own published regression line (y = -16.4x - 16.5,
    # r^2=0.994, x = log10[AOT mol/L]), 3 post-CMC plateau points estimated
    # by eye from the same figure's flat region. Paper's own reported CMC
    # for AOT in water: 2.51 mM.
    concentrations = [0.01585, 0.03981, 0.10000, 0.25119, 0.63096, 1.58489,
                       3.981, 6.31, 10.0]
    tensions = [62.22, 55.66, 49.10, 42.54, 35.98, 29.42,
                27.5, 27.8, 28.0]
    result = cmc_from_surface_tension_curve(concentrations, tensions)
    print(f"Extracted CMC: {result.cmc_mM:.3f} mM  (paper's reported value: 2.51 mM)")
    print(f"Premicellar slope: {result.premicellar_slope_mN_per_m_per_lnC:.3f} mN/m per ln(C)")
    print(f"Gamma at CMC: {result.gamma_at_cmc_mN_per_m:.2f} mN/m")
    print(f"R^2 (premicellar fit): {result.r_squared_premicellar:.4f}")
    print(f"Split: {result.n_premicellar_points} premicellar / {result.n_postmicellar_points} post-CMC points")
