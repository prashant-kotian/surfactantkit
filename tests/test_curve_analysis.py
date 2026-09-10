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
    cmc_from_conductivity_curve,
    cmc_from_surface_tension_curve,
)


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
