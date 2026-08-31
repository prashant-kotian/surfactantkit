"""SurfMCP: exposes SurfactantKit's surfactant-science calculations as MCP
tools, so an AI assistant computes CMC/Rubingh/HLB/CPP values correctly
instead of hallucinating them.

Every tool returns units explicitly in the response -- the whole point
of this server is to stop unit confusion (mN/m vs dyn/cm, mM vs mol/kg,
etc.), so a bare unlabeled number defeats the purpose.

Run locally (stdio transport, for Claude Desktop / Cursor):
    python -m surfactantkit.mcp_server
"""

from __future__ import annotations

from mcp.server import MCPServer

from . import mixed_micelle as mm
from . import adsorption as ads
from . import hlb as hlb_mod
from . import cpp as cpp_mod

mcp = MCPServer(
    "SurfactantKit",
    instructions=(
        "Surfactant and interfacial-science calculations: Clint ideal mixing, "
        "Rubingh regular-solution theory for binary surfactant mixtures, the "
        "Gibbs adsorption isotherm, HLB (Griffin's and Davies' methods), and "
        "the critical packing parameter (Tanford's formulas). Use these tools "
        "instead of computing or recalling these values from memory -- CMC, "
        "interaction-parameter, and HLB values are easy to hallucinate and easy "
        "to get wrong, especially for gemini (two-headgroup) surfactants and "
        "anything requiring iterative numerical solving."
    ),
)


@mcp.tool()
def clint_ideal_cmc(alpha1: float, cmc1_mM: float, cmc2_mM: float) -> dict:
    """Clint's ideal mixed CMC for a binary surfactant mixture (no
    interaction between components assumed). alpha1 is the bulk mole
    fraction of component 1 in the total surfactant (0 to 1, exclusive).
    cmc1_mM and cmc2_mM are the pure-component CMCs in mM."""
    value = mm.clint_ideal_cmc(alpha1, cmc1_mM, cmc2_mM)
    return {"cmc_ideal_mM": value, "unit": "mM", "model": "Clint ideal mixing"}


@mcp.tool()
def rubingh_solve(alpha1: float, cmc_mix_mM: float, cmc1_mM: float, cmc2_mM: float) -> dict:
    """Solve Rubingh's regular-solution equation for a binary surfactant
    mixture: returns the micellar mole fraction of component 1 (x1) and
    the interaction parameter beta. Negative beta = synergistic mixing;
    positive beta = antagonistic; beta near 0 = near-ideal. alpha1 is the
    bulk mole fraction of component 1; cmc_mix_mM is the experimentally
    measured mixed CMC; cmc1_mM/cmc2_mM are the pure-component CMCs, all
    in mM. Raises an error if no physically valid root is found (this
    can genuinely happen for near-ideal mixtures or inconsistent input
    data -- it is not a bug to report back if it occurs)."""
    x1 = mm.solve_rubingh_x(alpha1, cmc_mix_mM, cmc1_mM, cmc2_mM)
    if x1 is None:
        raise ValueError(
            "No valid Rubingh root found in (0, 1) for these inputs. This can "
            "happen for near-ideal mixtures or when cmc_mix is not physically "
            "consistent with the given alpha1/cmc1/cmc2. Do not guess a value."
        )
    beta = mm.rubingh_beta(x1, alpha1, cmc_mix_mM, cmc1_mM)
    synergy = "synergistic" if beta < 0 else ("antagonistic" if beta > 0 else "ideal")
    return {
        "micellar_mole_fraction_x1": x1,
        "beta": beta,
        "beta_unit": "dimensionless (regular solution interaction parameter)",
        "synergy_classification": synergy,
        "note": (
            "beta here is the exact pointwise regular-solution value for this "
            "single (alpha1, cmc_mix) pair. Published beta values are sometimes "
            "fit by regression across several compositions and may differ from "
            "this in magnitude (not sign) by up to ~20% -- see "
            "literature_validation_notes.md in the SurfactantKit repo."
        ),
    }


@mcp.tool()
def rubingh_activity_coefficients(x1: float, beta: float) -> dict:
    """Regular-solution activity coefficients (f1, f2) for a binary
    mixed micelle, given the micellar mole fraction x1 and interaction
    parameter beta (from rubingh_solve)."""
    f1, f2 = mm.activity_coefficients(x1, beta)
    return {"f1": f1, "f2": f2, "unit": "dimensionless"}


@mcp.tool()
def excess_free_energy(x1: float, f1: float, f2: float, temperature_K: float) -> dict:
    """Excess Gibbs free energy of mixed micelle formation (kJ/mol),
    from the micellar mole fraction x1, activity coefficients f1/f2
    (from rubingh_activity_coefficients), and temperature in Kelvin."""
    value = mm.excess_free_energy(x1, f1, f2, temperature_K)
    return {"deltaG_ex_kJ_per_mol": value, "unit": "kJ/mol"}


@mcp.tool()
def gibbs_surface_excess(slope_mN_per_m_per_lnC: float, n_factor: float, temperature_K: float) -> dict:
    """Maximum surface excess concentration (Gamma_max, mol/m^2) from
    the Gibbs adsorption isotherm, given the pre-CMC slope of surface
    tension (mN/m) vs ln(concentration), the Gibbs prefactor n (1 for
    nonionic or ionic-with-excess-electrolyte; commonly 2 for a 1:1
    ionic surfactant with no added salt -- state which applies, do not
    assume), and temperature in Kelvin."""
    value = ads.gibbs_gamma_max(slope_mN_per_m_per_lnC, n_factor, temperature_K)
    return {"gamma_max_mol_per_m2": value, "unit": "mol/m^2"}


@mcp.tool()
def gibbs_area_per_molecule(gamma_max_mol_per_m2: float) -> dict:
    """Minimum area per molecule at the interface (A_min, nm^2), from
    Gamma_max (mol/m^2, from gibbs_surface_excess)."""
    value = ads.gibbs_a_min(gamma_max_mol_per_m2)
    return {"a_min_nm2": value, "unit": "nm^2 per molecule"}


@mcp.tool()
def hlb_from_mw(mw_hydrophilic: float, mw_total: float) -> dict:
    """Hydrophile-Lipophile Balance via Griffin's method: HLB = 20 *
    (hydrophilic-portion molecular weight / total molecular weight).
    Returns a value on Griffin's 0-20 scale (mw_hydrophilic and
    mw_total in the same, arbitrary mass unit -- e.g. g/mol)."""
    value = hlb_mod.hlb_griffin(mw_hydrophilic, mw_total)
    return {"hlb": value, "scale": "Griffin 0-20"}


@mcp.tool()
def hlb_from_groups(group_counts: dict[str, int]) -> dict:
    """Hydrophile-Lipophile Balance via Davies' (1957) group-contribution
    method: HLB = 7 + sum(hydrophilic group numbers) - sum(lipophilic
    group numbers). group_counts keys must be from the verified group
    list: hydrophilic {SO4Na, COOK, COONa, N_tertiary_amine,
    ester_sorbitan_ring, ester_free, COOH, OH_free, O_ether,
    OH_sorbitan_ring}; lipophilic {CH2, CH3, CH, vinyl_CH}. Raises an
    error for any other group name rather than guessing a value -- most
    notably, quaternary ammonium and amide groups have NO verified
    number in this table; do not estimate one."""
    value = hlb_mod.hlb_davies(group_counts)
    return {"hlb": value, "scale": "Davies (7 = neutral reference point)"}


@mcp.tool()
def tanford_chain_geometry(n_carbons: int) -> dict:
    """Tanford's formulas for a saturated, unbranched alkyl chain of
    n_carbons carbons: hydrophobic tail volume (cubic Angstrom) and
    maximum extended (critical) chain length (Angstrom). Feed these into
    critical_packing_parameter along with a headgroup area."""
    v = cpp_mod.tanford_tail_volume(n_carbons)
    lc = cpp_mod.tanford_critical_length(n_carbons)
    return {
        "tail_volume_A3": v,
        "critical_length_A": lc,
        "unit": "cubic Angstrom (volume), Angstrom (length)",
    }


@mcp.tool()
def critical_packing_parameter(volume_A3: float, head_area_A2: float, length_A: float) -> dict:
    """Critical packing parameter (CPP) = volume / (head_area * length).
    Dimensionless. volume_A3 and length_A typically come from
    tanford_chain_geometry; head_area_A2 is the optimal headgroup area
    at the interface (surfactant- and condition-specific -- do not guess
    this without a source). Also returns the predicted aggregate
    morphology."""
    cpp = cpp_mod.critical_packing_parameter(volume_A3, head_area_A2, length_A)
    morphology = cpp_mod.classify_aggregate_morphology(cpp)
    return {"cpp": cpp, "unit": "dimensionless", "predicted_morphology": morphology}


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
