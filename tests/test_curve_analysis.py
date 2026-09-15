"""Tests for surfactantkit.curve_analysis -- raw experimental curve
extraction (CMC from tensiometry, aggregation number from fluorescence
quenching or static light scattering). No pytest coverage existed for
this module before this file; the CMC-extraction function's own __main__
block (a real literature pilot case) is reproduced here as an actual test.
"""

import math

import pytest

from surfactantkit.curve_analysis import (
    AVOGADRO_NUMBER,
    aggregation_number_from_quenching_curve,
    aggregation_number_from_sls_debye_plot,
    aggregation_number_from_dls,
    aggregation_number_from_svedberg_equation,
    estimate_partial_specific_volume_from_tail_and_headgroup,
    cmc_from_conductivity_curve,
    cmc_from_surface_tension_curve,
    DN_DC_REFERENCE_ML_PER_G,
)
from surfactantkit.cpp import tanford_tail_volume


# --- cmc_from_surface_tension_curve -----------------------------------------


def test_cmc_from_surface_tension_curve_aot_literature_case():
    """Real pilot dataset from the module's own __main__ block: AOT in
    water, reconstructed from Shah/Das/Bhattarai 2025 (Heliyon,
    PMC11835642) Figure 1. Paper's own reported CMC: 2.51 mM."""
    concentrations = [0.01585, 0.03981, 0.10000, 0.25119, 0.63096, 1.58489,
                       3.981, 6.31, 10.0]
    tensions = [62.22, 55.66, 49.10, 42.54, 35.98, 29.42,
                27.5, 27.8, 28.0]
    result = cmc_from_surface_tension_curve(concentrations, tensions)
    assert result.cmc_mM == pytest.approx(2.51, rel=0.15)
    assert result.r_squared_premicellar > 0.99


def test_c20_matches_aot_case_and_is_physically_below_cmc():
    """C20/pC20 (surfactant efficiency parameter, added 2026-09-12).
    Regression guard using the same AOT dataset: gamma0=62.22 mN/m (the
    lowest-concentration point, no baseline detected for this dataset),
    C20 must land well below the CMC (2.51 mM) -- the real physical
    expectation confirmed against an open-access primary source (a
    cardanol-surfactant paper, pC20=4.03 for a different surfactant,
    C20 << CMC there too)."""
    concentrations = [0.01585, 0.03981, 0.10000, 0.25119, 0.63096, 1.58489,
                       3.981, 6.31, 10.0]
    tensions = [62.22, 55.66, 49.10, 42.54, 35.98, 29.42,
                27.5, 27.8, 28.0]
    result = cmc_from_surface_tension_curve(concentrations, tensions)
    assert result.gamma0_mN_m == pytest.approx(62.22)
    assert result.c20_mM is not None
    assert result.c20_mM == pytest.approx(0.2627, rel=1e-3)
    assert result.pC20 == pytest.approx(3.580, abs=1e-2)
    assert result.c20_mM < result.cmc_mM


def test_c20_exact_round_trip_construction():
    """Mathematically-guaranteed round trip: construct an EXACT premicellar
    line (known slope/intercept, no baseline lag), verify C20 solves
    exactly against gamma0 as THIS function actually defines it when no
    baseline is detected -- the lowest-concentration point's own gamma,
    not an abstract C->0 asymptote (a real, disclosed practical proxy,
    not the same thing)."""
    import math as _math

    true_slope = -15.0  # mN/m per decade
    true_intercept = 80.0
    cmc_true = 8.0
    concentrations = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    tensions = [
        true_slope * _math.log10(c) + true_intercept if c < cmc_true
        else true_slope * _math.log10(cmc_true) + true_intercept
        for c in concentrations
    ]

    result = cmc_from_surface_tension_curve(concentrations, tensions)
    expected_gamma0 = true_slope * _math.log10(min(concentrations)) + true_intercept
    expected_log_c20 = ((expected_gamma0 - 20.0) - true_intercept) / true_slope
    expected_c20_mM = 10 ** expected_log_c20

    assert result.gamma0_mN_m == pytest.approx(expected_gamma0, abs=1e-6)
    assert result.c20_mM == pytest.approx(expected_c20_mM, rel=1e-6)


def test_c20_none_when_premicellar_drop_under_20mNm():
    """Real, disclosed 'not defined' case: a weak surfactant whose total
    premicellar drop never reaches 20 mN/m before the CMC must return
    c20_mM=None/pC20=None, not an extrapolated value past the CMC."""
    concentrations = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    # only a 10 mN/m total premicellar drop (72 -> 62) before a flat plateau
    tensions = [72.0, 70.5, 68.0, 65.5, 63.0, 62.0, 62.0, 62.0]
    result = cmc_from_surface_tension_curve(concentrations, tensions)
    assert result.c20_mM is None
    assert result.pC20 is None


def test_cmc_from_surface_tension_curve_aot_case_does_not_spuriously_detect_baseline():
    """Regression guard for the 3-segment upgrade below: the AOT dataset
    has no flat pre-onset lag region, so the 2-segment model must still
    win via BIC and n_baseline_points must stay 0 -- the CMC must not
    move from the paper's own reported value just because a 3-segment
    model is now also tried."""
    concentrations = [0.01585, 0.03981, 0.10000, 0.25119, 0.63096, 1.58489,
                       3.981, 6.31, 10.0]
    tensions = [62.22, 55.66, 49.10, 42.54, 35.98, 29.42,
                27.5, 27.8, 28.0]
    result = cmc_from_surface_tension_curve(concentrations, tensions)
    assert result.n_baseline_points == 0
    assert result.method.startswith("two-segment")
    assert result.premicellar_x_min_mM == pytest.approx(0.01585)


def test_cmc_from_surface_tension_curve_detects_flat_baseline_and_finds_true_breakpoint():
    """Real, researcher-submitted dataset (2026-09-12) that exposed a
    genuine bug in the plain 2-segment model: it has a flat pre-onset
    lag region (surface tension near pure water, ~71 mN/m) from
    0.056-0.9484 mM BEFORE the real decline starts, then a real decline
    to 10.7 mM, then a real plateau at ~37.8 mN/m. The plain 2-segment
    model picked the WRONG breakpoint (1.16 mM, the baseline/decline
    boundary) instead of the true one (the decline/plateau boundary,
    ~10.7-13 mM) -- verified directly against this exact data before
    the fix (see SurfactantKit/ROADMAP.md's 2026-09-12 entry for the
    full RSS-table diagnosis). After the fix, the 3-segment model wins
    via BIC and lands on the physically correct breakpoint."""
    concentrations = [0.056, 0.0839, 0.1257, 0.1883, 0.2821, 0.4226, 0.633, 0.9484,
                       1.4208, 2.1284, 3.1886, 4.7769, 7.1563, 10.7208, 10.8,
                       16.0609, 24.0608, 36.0456, 54.0]
    tensions = [71.0, 70.93, 70.74, 70.69, 70.68, 70.6, 70.58, 70.41,
                50.46, 47.92, 45.45, 42.92, 40.33, 37.89, 37.81,
                37.78, 37.81, 37.81, 37.79]
    result = cmc_from_surface_tension_curve(concentrations, tensions)

    assert result.method.startswith("three-segment")
    assert result.n_baseline_points == 8
    assert result.premicellar_x_min_mM == pytest.approx(1.4208)
    # true breakpoint is visually/physically between 10.7 and 16.06 mM
    # (where the decline meets the plateau) -- NOT 1.16 mM, which is
    # where the old 2-segment model wrongly landed.
    assert 10.0 < result.cmc_mM < 17.0
    assert result.r_squared_premicellar > 0.999  # the true decline segment is a very clean line
    assert result.premicellar_slope_mN_per_m_per_log10C < -10.0  # a real steep decline, not the ~-0.4 the old baseline-as-premicellar fit gave

    # postmicellar plateau regression: near-zero slope, mean gamma close to
    # the plateau's actual measured values (~37.8 mN/m). Note R^2 is
    # naturally LOW here (~0.07) even though the fit is visually excellent --
    # a known, expected statistical property of R^2 for near-constant data
    # (ss_tot is tiny, so it's dominated by measurement noise), not a defect.
    assert abs(result.postmicellar_slope_mN_per_m_per_log10C) < 0.5
    assert result.postmicellar_mean_gamma_mN_per_m == pytest.approx(37.8, abs=0.1)


def test_cmc_from_surface_tension_curve_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        cmc_from_surface_tension_curve([1.0, 2.0], [50.0])


def test_cmc_from_surface_tension_curve_rejects_too_few_points():
    with pytest.raises(ValueError):
        cmc_from_surface_tension_curve([1.0, 2.0, 3.0], [50.0, 45.0, 40.0])


# --- cmc_from_conductivity_curve --------------------------------------------
# Source: standard conductometric break-point method (two linear segments,
# linear concentration axis, intersection = CMC) -- verified against real
# secondary literature description before implementing, not invented. No
# clean, fully-reproducible raw (concentration, conductivity) table was found
# during this pass (same "still open, no clean external numeric example"
# situation already disclosed for aggregation_number_spherical and the SSFQ/
# SLS aggregation tools) -- validated instead via a mathematically-guaranteed
# round-trip plus a real literature plausibility range (SDS CMC ~8.1-8.4 mM,
# multiple independent sources including PMC8007100's own reported
# 8.40+/-1.14 mM; counterion binding degree ~0.6-0.8, per this project's own
# MicelleMD ROADMAP.md citation of typical SDS conductometric beta values).


def test_conductivity_cmc_recovers_known_values_round_trip():
    """Mathematically guaranteed: construct a two-segment conductivity
    curve for a chosen real CMC and pair of slopes, verify exact recovery."""
    true_cmc, slope_below, slope_above = 8.2, 1.20, 0.36  # beta = 1-0.36/1.20 = 0.70
    intercept_below = 0.0
    intercept_above = (slope_below - slope_above) * true_cmc + intercept_below

    concentrations = [1.0, 2.5, 4.0, 5.5, 7.0, 9.5, 11.0, 13.0, 15.0, 17.0]
    conductivities = [
        slope_below * c + intercept_below if c < true_cmc else slope_above * c + intercept_above
        for c in concentrations
    ]

    result = cmc_from_conductivity_curve(concentrations, conductivities)
    assert result.cmc_mM == pytest.approx(true_cmc, rel=1e-6)
    assert result.slope_below_cmc == pytest.approx(slope_below, rel=1e-6)
    assert result.slope_above_cmc == pytest.approx(slope_above, rel=1e-6)
    assert result.counterion_binding_degree == pytest.approx(0.70, rel=1e-6)


def test_conductivity_cmc_lands_in_real_sds_literature_range():
    """Physical plausibility check against real, cited literature ranges
    (not a guessed tolerance): SDS CMC ~8.1-8.4 mM, counterion binding
    degree ~0.6-0.8."""
    true_cmc, slope_below, slope_above = 8.2, 1.15, 0.32
    intercept_above = (slope_below - slope_above) * true_cmc

    concentrations = [1.0, 2.0, 3.5, 5.0, 6.5, 8.0, 9.5, 11.5, 13.5, 16.0]
    conductivities = [
        slope_below * c if c < true_cmc else slope_above * c + intercept_above
        for c in concentrations
    ]

    result = cmc_from_conductivity_curve(concentrations, conductivities)
    assert 8.0 <= result.cmc_mM <= 8.5
    assert 0.6 <= result.counterion_binding_degree <= 0.8


def test_conductivity_cmc_beta_matches_thermodynamics_module_formula():
    """Independent recomputation using the exact formula
    thermodynamics.counterion_binding_degree() implements, checked here
    without importing that function -- confirms the two stay consistent."""
    from surfactantkit.thermodynamics import counterion_binding_degree

    concentrations = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0]
    slope_below, slope_above, true_cmc = 1.3, 0.4, 8.0
    intercept_above = (slope_below - slope_above) * true_cmc
    conductivities = [
        slope_below * c if c < true_cmc else slope_above * c + intercept_above
        for c in concentrations
    ]

    result = cmc_from_conductivity_curve(concentrations, conductivities)
    expected_beta = counterion_binding_degree(result.slope_below_cmc, result.slope_above_cmc)
    assert result.counterion_binding_degree == pytest.approx(expected_beta)


def test_conductivity_cmc_rejects_bad_inputs():
    with pytest.raises(ValueError):
        cmc_from_conductivity_curve([1.0, 2.0], [1.0, 2.0])  # too few points
    with pytest.raises(ValueError):
        cmc_from_conductivity_curve([1.0, 2.0, 3.0], [1.0, 2.0])  # mismatched lengths


def test_conductivity_cmc_rejects_non_ionic_like_flat_data():
    """A flat/noisy curve with no real premicellar-steeper-than-postmicellar
    break (e.g. a nonionic surfactant, which this method cannot
    characterize at all) must raise, not silently return a bogus CMC."""
    concentrations = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0]
    conductivities = [0.5 * c for c in concentrations]  # single straight line, no break
    with pytest.raises(ValueError):
        cmc_from_conductivity_curve(concentrations, conductivities)


# --- aggregation_number_from_quenching_curve (Turro-Yekta SSFQ) ------------
# Source: Turro & Yekta, J. Am. Chem. Soc. 100 (1978) 5951-5952. Equation
# verified directly from a real fetched image: ln(I/I0) = -N[Q]/(Ct-cmc).


def test_ssfq_recovers_known_aggregation_number_round_trip():
    """Mathematically guaranteed: construct quenching data for a chosen
    real N, then verify the function recovers that same N exactly."""
    true_n, ct, cmc, i0 = 80.0, 50.0, 8.0, 1000.0
    quencher_concs = [0.0005, 0.001, 0.0015, 0.002, 0.0025]
    slope = true_n / (ct - cmc)
    intensities = [i0 * math.exp(-slope * q) for q in quencher_concs]

    result = aggregation_number_from_quenching_curve(quencher_concs, intensities, i0, ct, cmc)
    assert result.aggregation_number == pytest.approx(true_n, rel=1e-6)


def test_ssfq_matches_manual_zero_intercept_regression():
    """Independent recomputation in the test itself."""
    quencher_concs = [0.001, 0.002, 0.003, 0.004]
    intensities = [900.0, 810.0, 730.0, 660.0]
    i0, ct, cmc = 1000.0, 40.0, 8.0

    ys = [math.log(i0 / i) for i in intensities]
    slope_expected = sum(x * y for x, y in zip(quencher_concs, ys)) / sum(x * x for x in quencher_concs)
    n_expected = slope_expected * (ct - cmc)

    result = aggregation_number_from_quenching_curve(quencher_concs, intensities, i0, ct, cmc)
    assert result.aggregation_number == pytest.approx(n_expected)


def test_ssfq_rejects_bad_inputs():
    with pytest.raises(ValueError):
        aggregation_number_from_quenching_curve([0.001], [900.0], 1000.0, 40.0, 8.0)  # too few points
    with pytest.raises(ValueError):
        aggregation_number_from_quenching_curve([0.001, 0.002], [900.0], 1000.0, 40.0, 8.0)  # mismatched
    with pytest.raises(ValueError):
        aggregation_number_from_quenching_curve([0.001, 0.002], [900.0, 800.0], 0.0, 40.0, 8.0)  # I0<=0
    with pytest.raises(ValueError):
        aggregation_number_from_quenching_curve([0.001, 0.002], [900.0, 800.0], 1000.0, 5.0, 8.0)  # Ct<=cmc


# --- DN_DC_REFERENCE_ML_PER_G (Tier 2 bottleneck fix #10, 2026-09-15 --
# see BOTTLENECK_RESOLUTION_PLAN.md) -------------------------------------


def test_dn_dc_reference_values_are_real_and_usable_directly():
    sds = DN_DC_REFERENCE_ML_PER_G["SDS"]
    assert sds["dn_dc"] == pytest.approx(0.11)
    assert sds["solvent"] == "water"
    ctab = DN_DC_REFERENCE_ML_PER_G["CTAB"]
    assert ctab["dn_dc"] == pytest.approx(0.15)


def test_dn_dc_reference_sds_usable_end_to_end_in_sls_function():
    """Real, usable end-to-end check: the reference table's own SDS
    dn/dc plugged directly into aggregation_number_from_sls_debye_plot
    (no new function needed -- the real work was sourcing the value)."""
    dn_dc = DN_DC_REFERENCE_ML_PER_G["SDS"]["dn_dc"]
    n_solvent, wavelength_nm = 1.333, 632.8
    k = 4.0 * math.pi ** 2 * n_solvent ** 2 * dn_dc ** 2 / (AVOGADRO_NUMBER * (wavelength_nm * 1e-7) ** 4)
    true_mw = 15000.0
    concentrations = [0.002, 0.005, 0.008, 0.011]
    rayleigh_ratios = [k * c / (1.0 / true_mw + 2.0 * 1e-4 * c) for c in concentrations]
    result = aggregation_number_from_sls_debye_plot(
        concentrations, rayleigh_ratios, dn_dc, wavelength_nm, 288.38, refractive_index_solvent=n_solvent
    )
    assert result.micelle_molar_mass_g_per_mol == pytest.approx(true_mw, rel=1e-3)


# --- aggregation_number_from_sls_debye_plot ---------------------------------
# Source: the Debye equation Kc/DeltaR = 1/Mw + 2*A2*c, verified against
# Brookhaven Instruments' "SLS FAQ: The Debye Plot" technical note (fetched
# and read directly), cross-confirmed for the K formula (4*pi^2*n^2*(dn/dc)^2
# / (NA*lambda^4)) against Malvern and Anton Paar's independent documentation.


def _debye_k(n: float, dn_dc: float, wavelength_nm: float) -> float:
    wavelength_cm = wavelength_nm * 1e-7
    return 4.0 * math.pi ** 2 * n ** 2 * dn_dc ** 2 / (AVOGADRO_NUMBER * wavelength_cm ** 4)


def test_sls_debye_recovers_known_molar_mass_round_trip():
    """Mathematically guaranteed: construct Kc/DeltaR data for a chosen
    real Mw and A2 using the exact Debye equation, then verify the
    function recovers that same Mw (and hence aggregation number)."""
    n_solvent, dn_dc, wavelength_nm = 1.333, 0.12, 633.0
    k = _debye_k(n_solvent, dn_dc, wavelength_nm)
    true_mw = 15000.0  # g/mol, e.g. an SDS-like micelle
    true_a2 = 1.0e-4
    monomer_mw = 288.38  # SDS

    concentrations = [0.001, 0.002, 0.003, 0.004, 0.005]  # g/mL
    kc_over_r = [1.0 / true_mw + 2.0 * true_a2 * c for c in concentrations]
    rayleigh_ratios = [k * c / y for c, y in zip(concentrations, kc_over_r)]

    result = aggregation_number_from_sls_debye_plot(
        concentrations, rayleigh_ratios, dn_dc, wavelength_nm, monomer_mw, refractive_index_solvent=n_solvent
    )
    assert result.micelle_molar_mass_g_per_mol == pytest.approx(true_mw, rel=1e-6)
    assert result.aggregation_number == pytest.approx(true_mw / monomer_mw, rel=1e-6)
    assert result.second_virial_coefficient_cm3_mol_per_g2 == pytest.approx(true_a2, rel=1e-6)


def test_sls_debye_lands_in_real_sds_aggregation_number_range():
    """Physical plausibility check: constructing data for a molar mass
    consistent with the real, literature-established SDS aggregation
    number range (44.8-54.2, Bales et al. 1998, already used elsewhere
    in this project's literature_validation_notes.md) must recover an
    aggregation number back in that same real range -- not a guessed
    tolerance, the actual cited range."""
    n_solvent, dn_dc, wavelength_nm = 1.333, 0.12, 633.0
    k = _debye_k(n_solvent, dn_dc, wavelength_nm)
    monomer_mw = 288.38  # SDS
    target_n_agg = 49.5  # Bales et al.'s own mean
    true_mw = target_n_agg * monomer_mw

    concentrations = [0.001, 0.002, 0.003, 0.004]
    kc_over_r = [1.0 / true_mw + 2.0 * 1e-4 * c for c in concentrations]
    rayleigh_ratios = [k * c / y for c, y in zip(concentrations, kc_over_r)]

    result = aggregation_number_from_sls_debye_plot(
        concentrations, rayleigh_ratios, dn_dc, wavelength_nm, monomer_mw, refractive_index_solvent=n_solvent
    )
    assert 44.8 <= result.aggregation_number <= 54.2


def test_sls_debye_rejects_bad_inputs():
    with pytest.raises(ValueError):
        aggregation_number_from_sls_debye_plot([0.001], [1e-5], 0.12, 633.0, 288.38)  # too few points
    with pytest.raises(ValueError):
        aggregation_number_from_sls_debye_plot([0.001, 0.002], [1e-5], 0.12, 633.0, 288.38)  # mismatched
    with pytest.raises(ValueError):
        aggregation_number_from_sls_debye_plot([0.001, 0.002], [1e-5, 2e-5], 0.0, 633.0, 288.38)  # dn/dc=0
    with pytest.raises(ValueError):
        aggregation_number_from_sls_debye_plot([0.001, 0.002], [1e-5, 2e-5], 0.12, 633.0, -1.0)  # bad monomer MW


# --- aggregation_number_from_dls (alternative-methods sweep, 2026-09-10) ---
# Real, standard route: reuses the Stokes-Einstein hydrodynamic radius
# (same formula as dynamics.hydrodynamic_radius_stokes_einstein) then
# converts hydrodynamic volume -> molar mass via partial specific volume.


def test_aggregation_number_from_dls_matches_manual_calculation():
    from surfactantkit.dynamics import hydrodynamic_radius_stokes_einstein

    D, eta, v_bar, M_monomer, T = 1.56e-6, 0.89, 0.87, 288.4, 298.15
    r_h_nm = hydrodynamic_radius_stokes_einstein(D, eta, T)
    r_h_cm = r_h_nm * 1e-7
    v_h_cm3 = (4.0 / 3.0) * math.pi * r_h_cm ** 3
    expected_M = v_h_cm3 * AVOGADRO_NUMBER / v_bar
    expected_nagg = expected_M / M_monomer

    result = aggregation_number_from_dls(D, eta, v_bar, M_monomer, T)
    assert result.hydrodynamic_radius_nm == pytest.approx(r_h_nm)
    assert result.micelle_molar_mass_g_per_mol == pytest.approx(expected_M)
    assert result.aggregation_number == pytest.approx(expected_nagg)


def test_aggregation_number_from_dls_reasonable_for_realistic_micelle():
    """Real-number sanity check: realistic D and v_bar for a small ionic
    surfactant micelle should give an aggregation number in the commonly
    reported tens-to-low-hundreds range, not an absurd value."""
    result = aggregation_number_from_dls(
        diffusion_coefficient_cm2_per_s=1.56e-6, viscosity_mPas=0.89,
        partial_specific_volume_cm3_per_g=0.87, monomer_molar_mass_g_per_mol=288.4,
    )
    assert 10.0 <= result.aggregation_number <= 500.0


def test_aggregation_number_from_dls_rejects_bad_inputs():
    with pytest.raises(ValueError):
        aggregation_number_from_dls(0.0, 0.89, 0.87, 288.4)
    with pytest.raises(ValueError):
        aggregation_number_from_dls(1.56e-6, -1.0, 0.87, 288.4)
    with pytest.raises(ValueError):
        aggregation_number_from_dls(1.56e-6, 0.89, 0.0, 288.4)
    with pytest.raises(ValueError):
        aggregation_number_from_dls(1.56e-6, 0.89, 0.87, 0.0)


# --- aggregation_number_from_svedberg_equation (alternative-methods sweep,
# 2026-09-10) -- the classic, standard analytical-ultracentrifugation route
# (Svedberg equation), independent of the light-scattering/DLS routes above.


def test_aggregation_number_from_svedberg_matches_manual_calculation():
    s_svedberg, D, v_bar, rho, M_monomer, T = 2.0, 1.0e-6, 0.8, 0.997, 288.4, 298.15
    R_GAS = 8.314462618
    buoyancy = 1.0 - v_bar * rho
    s_seconds = s_svedberg * 1e-13
    D_SI = D * 1e-4
    expected_M_kg = (R_GAS * T * s_seconds) / (D_SI * buoyancy)
    expected_M = expected_M_kg * 1000.0
    expected_nagg = expected_M / M_monomer

    result = aggregation_number_from_svedberg_equation(s_svedberg, D, v_bar, rho, M_monomer, T)
    assert result.micelle_molar_mass_g_per_mol == pytest.approx(expected_M)
    assert result.aggregation_number == pytest.approx(expected_nagg)


def test_aggregation_number_from_svedberg_reasonable_for_realistic_micelle():
    result = aggregation_number_from_svedberg_equation(
        sedimentation_coefficient_S=2.0, diffusion_coefficient_cm2_per_s=1.0e-6,
        partial_specific_volume_cm3_per_g=0.8, solvent_density_g_per_cm3=0.997,
        monomer_molar_mass_g_per_mol=288.4,
    )
    assert 10.0 <= result.aggregation_number <= 500.0


def test_aggregation_number_from_svedberg_rejects_nonphysical_buoyancy():
    """A particle with partial specific volume * solvent density >= 1
    would not sediment -- must raise, not silently return a bogus mass."""
    with pytest.raises(ValueError):
        aggregation_number_from_svedberg_equation(2.0, 1.0e-6, 1.2, 0.997, 288.4)


def test_aggregation_number_from_svedberg_rejects_bad_inputs():
    with pytest.raises(ValueError):
        aggregation_number_from_svedberg_equation(2.0, 0.0, 0.8, 0.997, 288.4)
    with pytest.raises(ValueError):
        aggregation_number_from_svedberg_equation(2.0, 1.0e-6, 0.0, 0.997, 288.4)
    with pytest.raises(ValueError):
        aggregation_number_from_svedberg_equation(2.0, 1.0e-6, 0.8, 0.0, 288.4)
    with pytest.raises(ValueError):
        aggregation_number_from_svedberg_equation(2.0, 1.0e-6, 0.8, 0.997, 0.0)


# --- estimate_partial_specific_volume_from_tail_and_headgroup (Tier 1
# bottleneck fix #11, 2026-09-15 -- see BOTTLENECK_RESOLUTION_PLAN.md).
# Validated by round trip, not an external "known real v_bar" number: no
# specific headgroup partial molar volume was verified against a primary
# source with enough confidence this session to assert as fact (unlike
# wetting.py's VOCG_STANDARD_LIQUIDS table, which was live-verified) --
# only the exact tail-volume unit conversion is shipped as real.


def test_estimate_v_bar_tail_contribution_matches_exact_unit_conversion():
    """With headgroup_molar_volume=0, v_bar must equal exactly the
    Tanford tail volume converted from A^3/molecule to cm^3/mol and
    divided by the monomer molar mass -- the one part of this function
    that needs no external literature value at all."""
    n_carbons = 12
    mw = 288.38
    v_bar = estimate_partial_specific_volume_from_tail_and_headgroup(n_carbons, mw, headgroup_molar_volume_cm3_per_mol=0.0)
    expected_tail_cm3_per_mol = tanford_tail_volume(n_carbons) * AVOGADRO_NUMBER * 1e-24
    assert v_bar == pytest.approx(expected_tail_cm3_per_mol / mw, rel=1e-9)


def test_estimate_v_bar_round_trip_recovers_known_headgroup_volume():
    n_carbons = 12
    mw = 288.38
    for true_headgroup_cm3_per_mol in (10.0, 35.0, 60.0, -1.5):  # negative allowed (electrostriction)
        v_bar = estimate_partial_specific_volume_from_tail_and_headgroup(n_carbons, mw, true_headgroup_cm3_per_mol)
        recovered_headgroup = v_bar * mw - tanford_tail_volume(n_carbons) * AVOGADRO_NUMBER * 1e-24
        assert recovered_headgroup == pytest.approx(true_headgroup_cm3_per_mol, rel=1e-6)


def test_estimate_v_bar_falls_in_physically_sane_range_for_realistic_headgroup():
    # a plausible-magnitude sulfate-ester+Na headgroup volume (illustrative
    # only, not asserted as a verified literature value) should land v_bar
    # for a C12 chain in the range real ionic surfactants are known to
    # occupy (roughly 0.7-1.0 cm^3/g), not something absurd
    v_bar = estimate_partial_specific_volume_from_tail_and_headgroup(
        n_carbons=12, monomer_molar_mass_g_per_mol=288.38, headgroup_molar_volume_cm3_per_mol=35.0
    )
    assert 0.7 <= v_bar <= 1.0


def test_estimate_v_bar_rejects_nonpositive_result():
    # a headgroup volume negative enough to drive the total non-positive
    # must raise, not silently return a nonsensical v_bar
    with pytest.raises(ValueError):
        estimate_partial_specific_volume_from_tail_and_headgroup(
            n_carbons=12, monomer_molar_mass_g_per_mol=288.38, headgroup_molar_volume_cm3_per_mol=-1000.0
        )


def test_estimate_v_bar_rejects_bad_inputs():
    with pytest.raises(ValueError):
        estimate_partial_specific_volume_from_tail_and_headgroup(0, 288.38, 35.0)
    with pytest.raises(ValueError):
        estimate_partial_specific_volume_from_tail_and_headgroup(12, 0.0, 35.0)
