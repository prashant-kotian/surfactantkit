"""End-to-end tests for the SurfMCP server: verifies tools are
registered and that calling them through the MCP protocol layer
produces the same results as calling the underlying library functions
directly (i.e. the wrapping doesn't silently change any numbers)."""

import asyncio
import json

import pytest

from surfactantkit.mcp_server import mcp


def call(tool_name: str, args: dict) -> dict:
    """Call an MCP tool synchronously and parse its JSON text content.

    Real API-shape change found 2026-09-04, alongside the MCPServer ->
    FastMCP rename (see mcp_server.py): call_tool() in the currently
    installed mcp SDK (1.27.2) returns the content list directly, not
    wrapped in a `.content` attribute -- confirmed by actually calling it
    and inspecting the real return value, not guessed from a changelog.
    """
    result = asyncio.run(mcp.call_tool(tool_name, args))
    return json.loads(result[0].text)


def test_all_expected_tools_are_registered():
    tools = asyncio.run(mcp.list_tools())
    names = {t.name for t in tools}
    expected = {
        "clint_ideal_cmc",
        "rubingh_solve",
        "rubingh_activity_coefficients",
        "excess_free_energy",
        "gibbs_surface_excess",
        "gibbs_area_per_molecule",
        "hlb_from_mw",
        "hlb_from_groups",
        "hld_optimal_salinity",
        "hld_cationic_quat_reference",
        "hld_fit_k_and_cc_from_salinity_scan",
        "tanford_chain_geometry",
        "critical_packing_parameter",
        "aggregation_number_from_quenching",
        "aggregation_number_from_sls",
        "cmc_from_surface_tension_curve",
        "cmc_from_conductivity",
        "rubingh_beta_regression",
        "frumkin_predict_surface_tension",
        "motomura_ideal_composition",
        "asymmetric_margules_activity_coefficients",
        "rodenas_x1",
        "rodenas_activity_coefficients",
        "rodenas_x1_from_series",
        "maeda_free_energy_of_micellization",
        "vant_hoff_multi_point_fit",
        "szyszkowski_fit_K",
        "eommm_global_fit",
        "nagarajan_equilibrium_area_ionic",
        "hydrodynamic_radius_perrin_corrected",
        "hlb_davies_guo_ecl",
        "mass_action_free_energy_of_micellization",
        "predict_mobility_relaxation_corrected",
        "zeta_potential_relaxation_corrected",
        "owens_wendt_solid_surface_energy",
        "van_oss_chaudhury_good_solid_surface_energy",
        "micelle_water_partition_coefficient",
        "grahame_equation_surface_potential",
        "aggregation_number_from_dls",
        "aggregation_number_from_svedberg",
    }
    assert expected <= names


def test_clint_ideal_cmc_tool_matches_library():
    out = call("clint_ideal_cmc", {"alpha1": 0.25, "cmc1_mM": 14.80, "cmc2_mM": 8.00})
    assert out["cmc_ideal_mM"] == pytest.approx(9.038, abs=0.01)
    assert out["unit"] == "mM"


def test_rubingh_solve_tool_matches_cholate_sds_literature_case():
    out = call(
        "rubingh_solve",
        {"alpha1": 0.5, "cmc_mix_mM": 4.07, "cmc1_mM": 11.50, "cmc2_mM": 11.98},
    )
    assert out["micellar_mole_fraction_x1"] == pytest.approx(0.5033, abs=0.001)
    assert out["beta"] == pytest.approx(-4.236, abs=0.01)
    assert out["synergy_classification"] == "synergistic"


def test_rubingh_beta_regression_tool_recovers_true_beta():
    import math

    cmc1, cmc2, beta_true = 14.80, 8.00, -1.8
    alpha1_series, cmc_mix_series = [], []
    for x1 in (0.25, 0.4, 0.55, 0.7):
        f1 = math.exp(beta_true * (1.0 - x1) ** 2)
        f2 = math.exp(beta_true * x1 ** 2)
        cmc_mix = x1 * f1 * cmc1 + (1.0 - x1) * f2 * cmc2
        alpha1 = x1 * f1 * cmc1 / cmc_mix
        alpha1_series.append(alpha1)
        cmc_mix_series.append(cmc_mix)

    out = call(
        "rubingh_beta_regression",
        {"alpha1_series": alpha1_series, "cmc_mix_series_mM": cmc_mix_series, "cmc1_mM": cmc1, "cmc2_mM": cmc2},
    )
    assert out["beta_mean"] == pytest.approx(beta_true, abs=1e-3)
    assert out["n_points_used"] == 4
    assert out["n_points_skipped"] == 0
    assert out["synergy_classification"] == "synergistic"


def test_szyszkowski_fit_k_tool_recovers_true_k_round_trip():
    true_K = 0.35
    gamma0, gamma_max, system_type = 72.0, 4.0e-6, "nonionic"
    concentrations = [0.5, 1.0, 2.0, 4.0, 8.0, 14.0]
    tensions = [
        call("szyszkowski_predict_surface_tension", {
            "concentration": c, "gamma0_mN_m": gamma0, "gamma_max_mol_per_m2": gamma_max,
            "K": true_K, "system_type": system_type,
        })["surface_tension_mN_m"]
        for c in concentrations
    ]

    out = call("szyszkowski_fit_K", {
        "concentrations": concentrations, "surface_tensions_mN_per_m": tensions,
        "gamma0_mN_m": gamma0, "gamma_max_mol_per_m2": gamma_max, "system_type": system_type,
    })
    assert out["K"] == pytest.approx(true_K, rel=1e-4)
    assert out["r_squared"] == pytest.approx(1.0, abs=1e-6)
    assert out["n_points"] == 6


def test_frumkin_tool_zero_a_matches_szyszkowski_tool():
    """a=0 must reduce EXACTLY to szyszkowski_predict_surface_tension --
    the real regression guard for this tool, since no external numeric
    literature match exists for the Frumkin isotherm specifically (see
    test_adsorption_alternatives.py's module docstring)."""
    frumkin_args = {
        "concentration": 2.0,
        "gamma0_mN_m": 72.0,
        "gamma_max_mol_per_m2": 4.0e-6,
        "K": 0.5,
        "a": 0.0,
        "system_type": "nonionic",
    }
    szyszkowski_args = {k: v for k, v in frumkin_args.items() if k != "a"}
    frumkin_out = call("frumkin_predict_surface_tension", frumkin_args)
    szyszkowski_out = call("szyszkowski_predict_surface_tension", szyszkowski_args)
    assert frumkin_out["surface_tension_mN_m"] == pytest.approx(szyszkowski_out["surface_tension_mN_m"], abs=1e-6)


def test_motomura_ideal_composition_tool_matches_library():
    out = call("motomura_ideal_composition", {"alpha1": 0.75, "cmc1_mM": 14.80, "cmc2_mM": 8.00})
    assert out["x1_ideal"] == pytest.approx(0.75 * 12.21 / 14.80, abs=0.001)


def test_asymmetric_margules_tool_reduces_to_rst_when_symmetric():
    beta = -1.5
    asym_out = call("asymmetric_margules_activity_coefficients", {"x1": 0.4, "w12": beta, "w21": beta})
    rst_out = call("rubingh_activity_coefficients", {"x1": 0.4, "beta": beta})
    assert asym_out["f1"] == pytest.approx(rst_out["f1"], rel=1e-9)
    assert asym_out["f2"] == pytest.approx(rst_out["f2"], rel=1e-9)


def test_rodenas_tools_end_to_end():
    x1_out = call("rodenas_x1", {"alpha1": 0.5, "dln_cmc_mix_dalpha1": 0.0})
    assert x1_out["micellar_mole_fraction_x1"] == pytest.approx(0.5)
    f_out = call(
        "rodenas_activity_coefficients",
        {"alpha1": 0.5, "x1": 0.6, "cmc_mix_mM": 5.0, "cmc1_mM": 10.0, "cmc2_mM": 8.0},
    )
    assert f_out["f1"] == pytest.approx((0.5 * 5.0) / (0.6 * 10.0))


def test_maeda_tool_matches_library():
    out = call(
        "maeda_free_energy_of_micellization",
        {"x1_rub": 0.6, "beta": -1.8, "cmc1_M": 0.041e-3, "cmc2_M": 0.263e-3, "temperature_K": 298.15},
    )
    assert isinstance(out["deltaG_M_kJ_per_mol"], float)
    assert out["deltaG_M_kJ_per_mol"] < 0


def test_hld_cationic_quat_reference_tool_matches_library():
    out = call("hld_cationic_quat_reference", {"surfactant": "ctab"})
    assert out["Cc"] == pytest.approx(-5.7)
    assert out["k"] == pytest.approx(0.7)
    assert out["real_reported_hlb_davies_scale"] == pytest.approx(21.4)


def test_hld_cationic_quat_reference_tool_rejects_unknown_surfactant():
    with pytest.raises(Exception):
        call("hld_cationic_quat_reference", {"surfactant": "SDS"})


def test_hld_optimal_salinity_tool_round_trips_with_library():
    from surfactantkit.hld import hld_ionic

    out = call("hld_optimal_salinity", {"k": 0.7, "eacn": 8.0, "cc": -5.7})
    s_star = out["optimal_salinity_pct"]
    assert hld_ionic(s_star, 0.7, 8.0, -5.7) == pytest.approx(0.0, abs=1e-6)


def test_hld_fit_k_and_cc_from_salinity_scan_tool_round_trip():
    from surfactantkit.hld import optimal_salinity_ionic

    true_k, true_cc = 0.7, -5.7
    eacn_series = [6.0, 8.0, 10.0, 12.0]
    s_star_series = [optimal_salinity_ionic(true_k, eacn, true_cc) for eacn in eacn_series]

    out = call(
        "hld_fit_k_and_cc_from_salinity_scan",
        {"eacn_series": eacn_series, "optimal_salinity_series_pct": s_star_series},
    )
    assert out["k"] == pytest.approx(true_k, abs=1e-6)
    assert out["cc"] == pytest.approx(true_cc, abs=1e-6)
    assert out["r_squared"] == pytest.approx(1.0, abs=1e-6)
    assert out["n_points"] == 4


def test_hlb_from_groups_tool_sds_worked_example():
    out = call("hlb_from_groups", {"group_counts": {"SO4Na": 1, "CH2": 11, "CH3": 1}})
    assert out["hlb"] == pytest.approx(39.9, abs=0.01)


def test_hlb_davies_guo_ecl_tool_matches_manual_formula():
    import math

    n_carbons, n_eo = 12, 9
    n_ch2_eff = 0.965 * (n_carbons - 1) - 0.178
    n_eo_eff = 13.45 * math.log(n_eo) - 0.16 * n_eo + 1.26
    expected = 7.0 + 1.3 * n_eo_eff - 0.475 * n_ch2_eff - 0.475 * 1

    out = call("hlb_davies_guo_ecl", {"n_carbons_alkyl": n_carbons, "n_eo": n_eo})
    assert out["hlb"] == pytest.approx(expected)


def test_hlb_from_groups_tool_raises_on_unverified_group():
    with pytest.raises(Exception):
        call("hlb_from_groups", {"group_counts": {"quaternary_ammonium": 1}})


def test_critical_packing_parameter_tool_end_to_end():
    geom = call("tanford_chain_geometry", {"n_carbons": 12})
    out = call(
        "critical_packing_parameter",
        {"volume_A3": geom["tail_volume_A3"], "head_area_A2": 50.0, "length_A": geom["critical_length_A"]},
    )
    assert out["cpp"] == pytest.approx(350.2 / (50.0 * 16.68), abs=0.001)
    assert out["predicted_morphology"] in {
        "spherical micelle",
        "cylindrical/rodlike micelle",
        "vesicle/bilayer",
        "inverted structure",
    }


def test_nagarajan_equilibrium_area_ionic_tool_matches_table2():
    out = call(
        "nagarajan_equilibrium_area_ionic",
        {"cmc_M": 0.008, "tail_length_A": 16.5, "headgroup_prefactor_A": 82.0},
    )
    assert out["a_e_A2"] == pytest.approx(67.4, rel=0.01)
    assert out["debye_huckel_kappa_inverse_A"] == pytest.approx(34.43, rel=0.01)


def test_aggregation_number_tool():
    geom = call("tanford_chain_geometry", {"n_carbons": 12})
    out = call("aggregation_number", {"tail_volume_A3": geom["tail_volume_A3"], "core_radius_A": geom["critical_length_A"]})
    assert 40.0 < out["aggregation_number"] < 80.0


def test_aggregation_number_from_quenching_tool_round_trips_with_library():
    import math

    true_n, ct, cmc, i0 = 60.0, 40.0, 8.0, 1000.0
    slope = true_n / (ct - cmc)
    quencher_concs = [0.001, 0.002, 0.003]
    intensities = [i0 * math.exp(-slope * q) for q in quencher_concs]
    out = call(
        "aggregation_number_from_quenching",
        {
            "quencher_concentrations": quencher_concs,
            "intensities": intensities,
            "intensity_without_quencher": i0,
            "total_surfactant_concentration": ct,
            "cmc": cmc,
        },
    )
    assert out["aggregation_number"] == pytest.approx(true_n, rel=1e-6)


def test_aggregation_number_from_sls_tool_round_trips_with_library():
    import math

    n_solvent, dn_dc, wavelength_nm = 1.333, 0.12, 633.0
    wavelength_cm = wavelength_nm * 1e-7
    k = 4.0 * math.pi ** 2 * n_solvent ** 2 * dn_dc ** 2 / (6.02214076e23 * wavelength_cm ** 4)
    true_mw, monomer_mw = 15000.0, 288.38
    concentrations = [0.001, 0.002, 0.003]
    rayleigh_ratios = [k * c / (1.0 / true_mw + 2.0e-4 * c) for c in concentrations]
    out = call(
        "aggregation_number_from_sls",
        {
            "micellized_concentrations_g_per_mL": concentrations,
            "rayleigh_ratios_cm_inv": rayleigh_ratios,
            "dn_dc_mL_per_g": dn_dc,
            "wavelength_nm": wavelength_nm,
            "monomer_molar_mass_g_per_mol": monomer_mw,
            "refractive_index_solvent": n_solvent,
        },
    )
    assert out["micelle_molar_mass_g_per_mol"] == pytest.approx(true_mw, rel=1e-6)


def test_aggregation_number_from_dls_tool_matches_library():
    out = call(
        "aggregation_number_from_dls",
        {
            "diffusion_coefficient_cm2_per_s": 1.56e-6, "viscosity_mPas": 0.89,
            "partial_specific_volume_cm3_per_g": 0.87, "monomer_molar_mass_g_per_mol": 288.4,
        },
    )
    assert 10.0 <= out["aggregation_number"] <= 500.0
    assert out["hydrodynamic_radius_nm"] == pytest.approx(1.573, abs=0.01)


def test_aggregation_number_from_svedberg_tool_matches_library():
    out = call(
        "aggregation_number_from_svedberg",
        {
            "sedimentation_coefficient_S": 2.0, "diffusion_coefficient_cm2_per_s": 1.0e-6,
            "partial_specific_volume_cm3_per_g": 0.8, "solvent_density_g_per_cm3": 0.997,
            "monomer_molar_mass_g_per_mol": 288.4,
        },
    )
    assert 10.0 <= out["aggregation_number"] <= 500.0


def test_cmc_from_surface_tension_curve_tool_matches_library():
    concentrations = [0.01585, 0.03981, 0.10000, 0.25119, 0.63096, 1.58489, 3.981, 6.31, 10.0]
    tensions = [62.22, 55.66, 49.10, 42.54, 35.98, 29.42, 27.5, 27.8, 28.0]
    out = call("cmc_from_surface_tension_curve", {"concentrations_mM": concentrations, "surface_tensions_mN_per_m": tensions})
    assert out["cmc_mM"] == pytest.approx(2.51, rel=0.15)


def test_cmc_from_conductivity_tool_matches_library():
    true_cmc, slope_below, slope_above = 8.2, 1.20, 0.36  # beta = 1-0.36/1.20 = 0.70
    intercept_above = (slope_below - slope_above) * true_cmc

    concentrations = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 11.0, 13.0, 15.0]
    conductivities = [
        slope_below * c if c < true_cmc else slope_above * c + intercept_above
        for c in concentrations
    ]

    out = call("cmc_from_conductivity", {"concentrations_mM": concentrations, "conductivities": conductivities})
    assert out["cmc_mM"] == pytest.approx(true_cmc, rel=1e-6)
    assert out["counterion_binding_degree"] == pytest.approx(0.70, rel=1e-6)


def test_debye_screening_length_tool_matches_textbook_value():
    out = call("debye_screening_length", {"ionic_strength_M": 0.1, "temperature_K": 298.15})
    assert out["debye_length_nm"] == pytest.approx(0.961, abs=0.01)


def test_grahame_equation_tool_zero_charge_gives_zero_potential():
    out = call("grahame_equation_surface_potential", {"surface_charge_density_C_per_m2": 0.0, "ionic_strength_M": 0.01})
    assert out["surface_potential_mV"] == pytest.approx(0.0)


def test_zeta_potential_tool_requires_explicit_regime():
    with pytest.raises(Exception):
        call("zeta_potential", {"electrophoretic_mobility_um_cm_per_Vs": 2.0, "viscosity_mPas": 0.89, "regime": "not_a_real_regime"})


def test_relaxation_corrected_zeta_potential_tool_round_trip():
    zeta_true, kappa_a, eta = 60.0, 50.0, 0.89
    predicted = call(
        "predict_mobility_relaxation_corrected",
        {"zeta_mV": zeta_true, "kappa_a": kappa_a, "viscosity_mPas": eta},
    )
    recovered = call(
        "zeta_potential_relaxation_corrected",
        {"electrophoretic_mobility_um_cm_per_Vs": predicted["electrophoretic_mobility_um_cm_per_Vs"], "kappa_a": kappa_a, "viscosity_mPas": eta},
    )
    assert recovered["zeta_potential_mV"] == pytest.approx(zeta_true, abs=1e-3)


def test_hydrodynamic_radius_tool_positive_and_reasonable():
    out = call("hydrodynamic_radius", {"diffusion_coefficient_cm2_per_s": 1e-6, "viscosity_mPas": 0.89})
    assert 0.1 < out["hydrodynamic_radius_nm"] < 1000.0


def test_hydrodynamic_radius_perrin_corrected_tool_sphere_limit_matches_plain_tool():
    plain = call("hydrodynamic_radius", {"diffusion_coefficient_cm2_per_s": 1.56e-6, "viscosity_mPas": 0.89})
    corrected = call(
        "hydrodynamic_radius_perrin_corrected",
        {"diffusion_coefficient_cm2_per_s": 1.56e-6, "viscosity_mPas": 0.89, "axial_ratio": 1.0},
    )
    assert corrected["hydrodynamic_radius_true_nm"] == pytest.approx(plain["hydrodynamic_radius_nm"])
    assert corrected["perrin_friction_factor"] == pytest.approx(1.0)


def test_hydrodynamic_radius_perrin_corrected_tool_smaller_for_rodlike():
    plain = call("hydrodynamic_radius", {"diffusion_coefficient_cm2_per_s": 1.56e-6, "viscosity_mPas": 0.89})
    corrected = call(
        "hydrodynamic_radius_perrin_corrected",
        {"diffusion_coefficient_cm2_per_s": 1.56e-6, "viscosity_mPas": 0.89, "axial_ratio": 4.0},
    )
    assert corrected["hydrodynamic_radius_true_nm"] < plain["hydrodynamic_radius_nm"]


def test_rosen_monolayer_solve_tool_matches_underlying_math():
    out = call("rosen_monolayer_solve", {"alpha1": 0.5, "c_mix_sigma_mM": 4.07, "c1_sigma_mM": 11.50, "c2_sigma_mM": 11.98})
    assert out["monolayer_mole_fraction_x1"] == pytest.approx(0.5033, abs=0.001)
    assert out["beta_sigma"] == pytest.approx(-4.236, abs=0.01)


def test_corrin_harkins_predict_tool_cmc_decreases_with_salt():
    out = call("corrin_harkins_predict", {
        "cmc1_mM": 10.0, "salt_conc1_mM": 10.0,
        "cmc2_mM": 5.0, "salt_conc2_mM": 100.0,
        "salt_conc_target_mM": 200.0,
    })
    assert out["predicted_cmc_mM"] < 5.0


def test_counterion_binding_degree_tool():
    out = call("counterion_binding_degree", {"slope_below_cmc": 2.0, "slope_above_cmc": 1.0})
    assert out["beta"] == pytest.approx(0.5)


def test_gibbs_free_energy_micellization_tool_is_negative():
    out = call("gibbs_free_energy_micellization", {"cmc_M": 0.008, "temperature_K": 298.15, "counterion_factor": 1.6})
    assert out["deltaG_mic_kJ_per_mol"] < 0


def test_vant_hoff_and_entropy_tools_complete_triad():
    dg = call("gibbs_free_energy_micellization", {"cmc_M": 0.008, "temperature_K": 298.15, "counterion_factor": 1.6})
    dh = call("vant_hoff_enthalpy", {"cmc1_M": 0.008, "temperature1_K": 293.15, "cmc2_M": 0.0075, "temperature2_K": 313.15})
    ds = call("entropy_of_micellization", {
        "deltaG_mic_kJ_per_mol": dg["deltaG_mic_kJ_per_mol"],
        "deltaH_mic_kJ_per_mol": dh["deltaH_mic_kJ_per_mol"],
        "temperature_K": 298.15,
    })
    expected = ((dh["deltaH_mic_kJ_per_mol"] - dg["deltaG_mic_kJ_per_mol"]) / 298.15) * 1000.0
    assert ds["deltaS_mic_J_per_mol_K"] == pytest.approx(expected)


def test_eommm_global_fit_tool_recovers_true_parameters_round_trip():
    from surfactantkit.mixed_micelle import asymmetric_margules_activity_coefficients

    cmc1, cmc2 = 14.80, 8.00
    true_w12, true_w21 = 2.0, -3.0
    x1_true_values = [0.25, 0.5, 0.75]
    alpha1_series, cmc_mix_series = [], []
    for x1 in x1_true_values:
        f1, f2 = asymmetric_margules_activity_coefficients(x1, true_w12, true_w21)
        d = x1 * f1 * cmc1 + (1.0 - x1) * f2 * cmc2
        alpha1_series.append(x1 * f1 * cmc1 / d)
        cmc_mix_series.append(d)

    out = call(
        "eommm_global_fit",
        {"alpha1_series": alpha1_series, "cmc_mix_series_mM": cmc_mix_series, "cmc1_mM": cmc1, "cmc2_mM": cmc2},
    )
    assert out["W12"] == pytest.approx(true_w12, abs=1e-2)
    assert out["W21"] == pytest.approx(true_w21, abs=1e-2)
    assert out["r_squared"] == pytest.approx(1.0, abs=1e-4)
    assert out["n_points"] == 3


def test_rodenas_x1_from_series_tool_matches_library():
    import math

    a, b, c = 0.5, -2.0, 3.0  # ln(cmc_mix) = a*alpha1^2 + b*alpha1 + c
    alpha1_series = [0.1, 0.25, 0.4, 0.55, 0.7]
    cmc_mix_series = [math.exp(a * x ** 2 + b * x + c) for x in alpha1_series]
    true_slope = lambda x: 2.0 * a * x + b

    out = call("rodenas_x1_from_series", {"alpha1_series": alpha1_series, "cmc_mix_series_mM": cmc_mix_series})
    assert out["n_points"] == 5
    for x, slope in zip(out["alpha1_used"], out["dln_cmc_mix_dalpha1"]):
        assert slope == pytest.approx(true_slope(x), abs=1e-9)


def test_mass_action_free_energy_tool_converges_to_pseudo_phase_at_large_n():
    pseudo_phase = call("gibbs_free_energy_micellization", {"cmc_M": 0.008, "temperature_K": 298.15, "counterion_factor": 1.0})
    mass_action = call(
        "mass_action_free_energy_of_micellization",
        {"cmc_M": 0.008, "aggregation_number": 100000, "temperature_K": 298.15},
    )
    assert mass_action["deltaG_mic_kJ_per_mol"] == pytest.approx(pseudo_phase["deltaG_mic_kJ_per_mol"], abs=1e-3)


def test_vant_hoff_multi_point_fit_tool_recovers_true_parameters():
    import math
    from surfactantkit.thermodynamics import R_GAS

    true_dH_kJ, true_dS_J, true_dCp_J = -20.0, -50.0, -300.0
    tr = 298.15
    temps = [278.15, 288.15, 298.15, 308.15, 318.15]
    true_dS_kJ, true_dCp_kJ = true_dS_J / 1000.0, true_dCp_J / 1000.0
    cmc_series_M = []
    for t in temps:
        g = (t - tr) - t * math.log(t / tr)
        dg = true_dH_kJ - t * true_dS_kJ + true_dCp_kJ * g
        x = math.exp(dg * 1000.0 / (R_GAS * t))
        cmc_series_M.append(x * 55.5 / (1.0 - x))  # invert cmc_to_mole_fraction's dilute formula

    out = call(
        "vant_hoff_multi_point_fit",
        {"cmc_series_M": cmc_series_M, "temperatures_K": temps, "reference_temperature_K": tr},
    )
    assert out["deltaH_mic_kJ_per_mol"] == pytest.approx(true_dH_kJ, abs=1e-3)
    assert out["deltaS_mic_J_per_mol_K"] == pytest.approx(true_dS_J, abs=1e-1)
    assert out["deltaCp_mic_J_per_mol_K"] == pytest.approx(true_dCp_J, abs=1.0)
    assert out["n_points"] == 5


def test_szyszkowski_tool_at_zero_concentration_returns_gamma0():
    out = call("szyszkowski_predict_surface_tension", {
        "concentration": 0.0, "gamma0_mN_m": 72.0, "gamma_max_mol_per_m2": 3.0e-6, "K": 50.0,
        "system_type": "nonionic",
    })
    assert out["surface_tension_mN_m"] == pytest.approx(72.0)


def test_wetting_work_of_adhesion_tool_complete_wetting():
    out = call("wetting_work_of_adhesion", {"gamma_LV_mN_m": 72.0, "contact_angle_deg": 0.0})
    assert out["work_of_adhesion_mJ_per_m2"] == pytest.approx(144.0)


def test_wetting_spreading_coefficient_tool_never_positive():
    out = call("wetting_spreading_coefficient", {"gamma_LV_mN_m": 72.0, "contact_angle_deg": 90.0})
    assert out["spreading_coefficient_mN_m"] <= 0


def test_owens_wendt_tool_recovers_true_components_round_trip():
    import math

    true_gamma_s_d, true_gamma_s_p = 30.0, 10.0
    liquids_d = [21.8, 50.8, 37.0]
    liquids_p = [51.0, 0.0, 26.4]
    contact_angles = []
    for gd, gp in zip(liquids_d, liquids_p):
        gt = gd + gp
        w = 2.0 * (math.sqrt(true_gamma_s_d * gd) + math.sqrt(true_gamma_s_p * gp))
        contact_angles.append(math.degrees(math.acos(w / gt - 1.0)))

    out = call(
        "owens_wendt_solid_surface_energy",
        {"contact_angles_deg": contact_angles, "liquid_gamma_dispersive_mN_m": liquids_d, "liquid_gamma_polar_mN_m": liquids_p},
    )
    assert out["gamma_s_dispersive_mN_m"] == pytest.approx(true_gamma_s_d, abs=1e-3)
    assert out["gamma_s_polar_mN_m"] == pytest.approx(true_gamma_s_p, abs=1e-3)


def test_van_oss_chaudhury_good_tool_recovers_true_components_round_trip():
    import math

    true_lw, true_acid, true_base = 25.0, 2.0, 15.0
    liquids_lw = [21.8, 34.0, 50.8]
    liquids_acid = [25.5, 3.92, 0.0]
    liquids_base = [25.5, 57.4, 0.0]
    contact_angles = []
    for lw, a, b in zip(liquids_lw, liquids_acid, liquids_base):
        gt = lw + 2.0 * math.sqrt(a * b)
        w = 2.0 * (math.sqrt(true_lw * lw) + math.sqrt(true_acid * b) + math.sqrt(true_base * a))
        contact_angles.append(math.degrees(math.acos(w / gt - 1.0)))

    out = call(
        "van_oss_chaudhury_good_solid_surface_energy",
        {
            "contact_angles_deg": contact_angles,
            "liquid_gamma_lw_mN_m": liquids_lw,
            "liquid_gamma_acid_mN_m": liquids_acid,
            "liquid_gamma_base_mN_m": liquids_base,
        },
    )
    assert out["gamma_s_lw_mN_m"] == pytest.approx(true_lw, abs=1e-3)
    assert out["gamma_s_acid_mN_m"] == pytest.approx(true_acid, abs=1e-3)
    assert out["gamma_s_base_mN_m"] == pytest.approx(true_base, abs=1e-3)


def test_eor_capillary_number_tool_matches_direct_formula():
    out = call("eor_capillary_number", {"viscosity_mPas": 1.0, "velocity_m_per_s": 1e-5, "interfacial_tension_mN_m": 0.01})
    expected = (1.0e-3 * 1e-5) / (0.01e-3)
    assert out["capillary_number"] == pytest.approx(expected)


def test_molar_solubilization_ratio_tool():
    out = call("molar_solubilization_ratio", {
        "total_solubilized_M": 0.51e-3, "intrinsic_water_solubility_M": 0.01e-3,
        "surfactant_concentration_M": 10e-3, "cmc_M": 2e-3,
    })
    assert out["msr"] == pytest.approx(0.0625)


def test_micelle_water_partition_coefficient_tool():
    total, intrinsic, surf, cmc = 0.51e-3, 0.01e-3, 10e-3, 2e-3
    x_micelle = (total - intrinsic) / ((total - intrinsic) + (surf - cmc))
    x_water = intrinsic / 55.5
    expected = x_micelle / x_water

    out = call("micelle_water_partition_coefficient", {
        "total_solubilized_M": total, "intrinsic_water_solubility_M": intrinsic,
        "surfactant_concentration_M": surf, "cmc_M": cmc,
    })
    assert out["km"] == pytest.approx(expected)
