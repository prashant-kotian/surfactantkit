"""End-to-end tests for the SurfMCP server: verifies tools are
registered and that calling them through the MCP protocol layer
produces the same results as calling the underlying library functions
directly (i.e. the wrapping doesn't silently change any numbers)."""

import asyncio
import json

import pytest

from surfactantkit.mcp_server import mcp


def call(tool_name: str, args: dict) -> dict:
    """Call an MCP tool synchronously and parse its JSON text content."""
    result = asyncio.run(mcp.call_tool(tool_name, args))
    return json.loads(result.content[0].text)


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
