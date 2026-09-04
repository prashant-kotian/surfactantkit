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
        "tanford_chain_geometry",
        "critical_packing_parameter",
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


def test_hlb_from_groups_tool_sds_worked_example():
    out = call("hlb_from_groups", {"group_counts": {"SO4Na": 1, "CH2": 11, "CH3": 1}})
    assert out["hlb"] == pytest.approx(39.9, abs=0.01)


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


def test_aggregation_number_tool():
    geom = call("tanford_chain_geometry", {"n_carbons": 12})
    out = call("aggregation_number", {"tail_volume_A3": geom["tail_volume_A3"], "core_radius_A": geom["critical_length_A"]})
    assert 40.0 < out["aggregation_number"] < 80.0


def test_debye_screening_length_tool_matches_textbook_value():
    out = call("debye_screening_length", {"ionic_strength_M": 0.1, "temperature_K": 298.15})
    assert out["debye_length_nm"] == pytest.approx(0.961, abs=0.01)


def test_zeta_potential_tool_requires_explicit_regime():
    with pytest.raises(Exception):
        call("zeta_potential", {"electrophoretic_mobility_um_cm_per_Vs": 2.0, "viscosity_mPas": 0.89, "regime": "not_a_real_regime"})


def test_hydrodynamic_radius_tool_positive_and_reasonable():
    out = call("hydrodynamic_radius", {"diffusion_coefficient_cm2_per_s": 1e-6, "viscosity_mPas": 0.89})
    assert 0.1 < out["hydrodynamic_radius_nm"] < 1000.0


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
