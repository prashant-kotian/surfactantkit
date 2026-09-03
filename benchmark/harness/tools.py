"""Tool registry for the SurfBench harness.

The 25 SurfactantKit tools, defined once as canonical JSON schemas backed by
the real library, plus converters to each provider's function-calling format
(OpenAI / Anthropic-style / Gemini) and a single dispatch() entry point.

This is deliberately independent of the MCP transport layer (surfactantkit.
mcp_server): the benchmark exposes the same tools to every model through that
provider's native function-calling interface, for uniform cross-vendor control.
The MCP server serves the identical tools to real end-users over MCP.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from surfactantkit import mixed_micelle as mm
from surfactantkit import adsorption as ads
from surfactantkit import hlb as hlb_mod
from surfactantkit import cpp as cpp_mod
from surfactantkit import electrostatics as elec
from surfactantkit import dynamics as dyn
from surfactantkit import thermodynamics as thermo
from surfactantkit import wetting
from surfactantkit import solubilization as solub


# --- dispatch functions (mirror the MCP wrappers' compute logic) ------------
# Each takes a kwargs dict and returns a JSON-serialisable dict with the value
# and an explicit unit -- a bare number would defeat the point of the tool.

def _clint(a):
    v = mm.clint_ideal_cmc(a["alpha1"], a["cmc1_mM"], a["cmc2_mM"])
    return {"cmc_ideal_mM": v, "unit": "mM"}

def _rubingh(a):
    x1 = mm.solve_rubingh_x(a["alpha1"], a["cmc_mix_mM"], a["cmc1_mM"], a["cmc2_mM"])
    if x1 is None:
        return {"error": "No valid Rubingh root in (0,1); cmc_mix not physically "
                "consistent with alpha1/cmc1/cmc2. Do not guess a value."}
    beta = mm.rubingh_beta(x1, a["alpha1"], a["cmc_mix_mM"], a["cmc1_mM"])
    syn = "synergistic" if beta < 0 else ("antagonistic" if beta > 0 else "ideal")
    return {"micellar_mole_fraction_x1": x1, "beta": beta, "synergy_classification": syn}

def _rubingh_activity(a):
    f1, f2 = mm.activity_coefficients(a["x1"], a["beta"])
    return {"f1": f1, "f2": f2, "unit": "dimensionless"}

def _excess_g(a):
    v = mm.excess_free_energy(a["x1"], a["f1"], a["f2"], a["temperature_K"])
    return {"delta_G_excess_kJ_per_mol": v, "unit": "kJ/mol"}

def _rosen(a):
    x1 = mm.solve_rosen_monolayer_x(a["alpha1"], a["c_mix_mM"], a["c1_mM"], a["c2_mM"])
    if x1 is None:
        return {"error": "No valid Rosen monolayer root in (0,1) for these inputs."}
    beta = mm.rosen_beta_sigma(x1, a["alpha1"], a["c_mix_mM"], a["c1_mM"])
    return {"monolayer_mole_fraction_x1": x1, "beta_sigma": beta, "unit": "dimensionless"}

def _corrin(a):
    pred, g = mm.corrin_harkins_predict_cmc(
        a["cmc1_mM"], a["salt_conc1_mM"], a["cmc2_mM"],
        a["salt_conc2_mM"], a["salt_conc_target_mM"])
    return {"predicted_cmc_mM": pred, "corrin_harkins_slope": g, "unit": "mM"}

def _gibbs_gamma(a):
    v = ads.gibbs_gamma_max(a["slope_mN_per_lnC"], a["n_factor"], a["temperature_K"])
    return {"gamma_max_mol_per_m2": v, "unit": "mol/m^2"}

def _gibbs_amin(a):
    v = ads.gibbs_a_min(a["gamma_max_mol_per_m2"])
    return {"a_min_nm2": v, "unit": "nm^2"}

def _szysz(a):
    v = ads.szyszkowski_surface_tension(
        a["concentration"], a["gamma0_mN_m"], a["gamma_max_mol_per_m2"],
        a["K"], a.get("temperature_K", 298.15), a.get("n_factor", 1.0))
    return {"surface_tension_mN_m": v, "unit": "mN/m"}

def _hlb_mw(a):
    v = hlb_mod.hlb_griffin(a["mw_hydrophilic"], a["mw_total"])
    return {"hlb": v, "method": "Griffin", "scale": "0-20"}

def _hlb_groups(a):
    v = hlb_mod.hlb_davies(a["group_counts"])
    return {"hlb": v, "method": "Davies"}

def _tanford(a):
    nc = a["n_carbons"]
    return {"tail_volume_A3": cpp_mod.tanford_tail_volume(nc),
            "critical_length_A": cpp_mod.tanford_critical_length(nc),
            "unit_volume": "cubic Angstrom", "unit_length": "Angstrom"}

def _cpp(a):
    v = cpp_mod.critical_packing_parameter(a["volume_A3"], a["head_area_A2"], a["length_A"])
    return {"cpp": v, "morphology": cpp_mod.classify_aggregate_morphology(v),
            "unit": "dimensionless"}

def _nagg(a):
    v = cpp_mod.aggregation_number_spherical(a["tail_volume_A3"], a["core_radius_A"])
    return {"aggregation_number": v, "note": "geometric estimate"}

def _debye(a):
    if a.get("ion_concentrations"):
        # keys arrive as strings over JSON; coerce charge numbers back to int
        conc = {int(z): c for z, c in a["ion_concentrations"].items()}
        I = elec.ionic_strength(conc)
    else:
        I = a["ionic_strength_M"]
    v = elec.debye_length(I, a.get("temperature_K", 298.15))
    return {"debye_length_nm": v, "ionic_strength_M": I, "unit": "nm"}

def _zeta(a):
    v = elec.zeta_potential_henry(
        a["electrophoretic_mobility_um_cm_per_Vs"], a["viscosity_mPas"], a["regime"])
    return {"zeta_potential_mV": v, "unit": "mV"}

def _rh(a):
    v = dyn.hydrodynamic_radius_stokes_einstein(
        a["diffusion_coefficient_cm2_per_s"], a["viscosity_mPas"],
        a.get("temperature_K", 298.15))
    return {"hydrodynamic_radius_nm": v, "unit": "nm"}

def _beta_bind(a):
    v = thermo.counterion_binding_degree(a["slope_below_cmc"], a["slope_above_cmc"])
    return {"counterion_binding_degree_beta": v, "unit": "dimensionless"}

def _dg_mic(a):
    x = thermo.cmc_to_mole_fraction(a["cmc_M"])
    v = thermo.gibbs_free_energy_micellization(x, a["temperature_K"], a.get("counterion_factor", 1.0))
    return {"delta_G_mic_kJ_per_mol": v, "cmc_mole_fraction": x, "unit": "kJ/mol"}

def _dh_mic(a):
    x1 = thermo.cmc_to_mole_fraction(a["cmc1_M"])
    x2 = thermo.cmc_to_mole_fraction(a["cmc2_M"])
    v = thermo.vant_hoff_enthalpy(x1, a["temperature1_K"], x2, a["temperature2_K"])
    return {"delta_H_mic_kJ_per_mol": v, "unit": "kJ/mol"}

def _ds_mic(a):
    v = thermo.entropy_micellization(
        a["deltaG_mic_kJ_per_mol"], a["deltaH_mic_kJ_per_mol"], a["temperature_K"])
    return {"delta_S_mic_J_per_mol_K": v, "unit": "J/(mol.K)"}

def _wa(a):
    v = wetting.work_of_adhesion(a["gamma_LV_mN_m"], a["contact_angle_deg"])
    return {"work_of_adhesion_mJ_per_m2": v, "unit": "mJ/m^2"}

def _spread(a):
    v = wetting.spreading_coefficient(a["gamma_LV_mN_m"], a["contact_angle_deg"])
    return {"spreading_coefficient_mN_per_m": v, "unit": "mN/m"}

def _ca(a):
    v = wetting.capillary_number(a["viscosity_mPas"], a["velocity_m_per_s"], a["interfacial_tension_mN_m"])
    regime = "capillary-dominated" if v < 1e-5 else "transitional/viscous-dominated"
    return {"capillary_number": v, "regime": regime, "unit": "dimensionless"}

def _msr(a):
    v = solub.molar_solubilization_ratio(
        a["total_solubilized_M"], a["intrinsic_water_solubility_M"],
        a["surfactant_concentration_M"], a["cmc_M"])
    return {"molar_solubilization_ratio": v, "unit": "dimensionless"}


def _num(desc):
    return {"type": "number", "description": desc}

def _spec(name, description, properties, required, fn):
    return {"name": name, "description": description,
            "parameters": {"type": "object", "properties": properties, "required": required},
            "fn": fn}


# --- canonical tool registry (name -> spec) ---------------------------------
_TOOL_LIST = [
    _spec("clint_ideal_cmc",
          "Clint ideal mixed CMC (mM) for a binary surfactant mixture, no interaction assumed.",
          {"alpha1": _num("bulk mole fraction of component 1 (0-1)"),
           "cmc1_mM": _num("pure CMC of component 1 (mM)"),
           "cmc2_mM": _num("pure CMC of component 2 (mM)")},
          ["alpha1", "cmc1_mM", "cmc2_mM"], _clint),
    _spec("rubingh_solve",
          "Solve Rubingh regular-solution theory: micellar mole fraction x1 and interaction "
          "parameter beta (negative=synergistic, positive=antagonistic). Returns an error if "
          "no physically valid root exists.",
          {"alpha1": _num("bulk mole fraction of component 1"),
           "cmc_mix_mM": _num("experimentally measured mixed CMC (mM)"),
           "cmc1_mM": _num("pure CMC of component 1 (mM)"),
           "cmc2_mM": _num("pure CMC of component 2 (mM)")},
          ["alpha1", "cmc_mix_mM", "cmc1_mM", "cmc2_mM"], _rubingh),
    _spec("rubingh_activity_coefficients",
          "Regular-solution activity coefficients f1, f2 from micellar mole fraction x1 and beta.",
          {"x1": _num("micellar mole fraction of component 1"),
           "beta": _num("interaction parameter")},
          ["x1", "beta"], _rubingh_activity),
    _spec("excess_free_energy",
          "Excess Gibbs free energy of mixed micelle formation (kJ/mol).",
          {"x1": _num("micellar mole fraction of component 1"),
           "f1": _num("activity coefficient of component 1"),
           "f2": _num("activity coefficient of component 2"),
           "temperature_K": _num("temperature in Kelvin")},
          ["x1", "f1", "f2", "temperature_K"], _excess_g),
    _spec("rosen_monolayer_solve",
          "Rosen mixed-monolayer theory (surface-tension data, NOT CMC): monolayer mole fraction "
          "x1 and interaction parameter beta^sigma.",
          {"alpha1": _num("bulk mole fraction of component 1"),
           "c_mix_mM": _num("mixture concentration reaching the reference surface tension (mM)"),
           "c1_mM": _num("component-1 concentration reaching that reference tension (mM)"),
           "c2_mM": _num("component-2 concentration reaching that reference tension (mM)")},
          ["alpha1", "c_mix_mM", "c1_mM", "c2_mM"], _rosen),
    _spec("corrin_harkins_predict",
          "Corrin-Harkins log-linear salt-CMC relation: predict CMC (mM) at a target counterion "
          "concentration from two (CMC, salt) points.",
          {"cmc1_mM": _num("CMC at salt_conc1 (mM)"),
           "salt_conc1_mM": _num("counterion concentration 1 (mM)"),
           "cmc2_mM": _num("CMC at salt_conc2 (mM)"),
           "salt_conc2_mM": _num("counterion concentration 2 (mM)"),
           "salt_conc_target_mM": _num("target counterion concentration (mM)")},
          ["cmc1_mM", "salt_conc1_mM", "cmc2_mM", "salt_conc2_mM", "salt_conc_target_mM"], _corrin),
    _spec("gibbs_surface_excess",
          "Maximum surface excess Gamma_max (mol/m^2) from the pre-CMC slope of surface tension "
          "vs ln(concentration).",
          {"slope_mN_per_lnC": _num("d(gamma)/d(lnC) in mN/m (negative for a surfactant)"),
           "n_factor": _num("Gibbs prefactor: 1 nonionic/excess salt, 2 for 1:1 ionic no salt"),
           "temperature_K": _num("temperature in Kelvin")},
          ["slope_mN_per_lnC", "n_factor", "temperature_K"], _gibbs_gamma),
    _spec("gibbs_area_per_molecule",
          "Minimum area per molecule A_min (nm^2) from Gamma_max (mol/m^2).",
          {"gamma_max_mol_per_m2": _num("maximum surface excess (mol/m^2)")},
          ["gamma_max_mol_per_m2"], _gibbs_amin),
    _spec("szyszkowski_predict_surface_tension",
          "Predict surface tension (mN/m) via the Szyszkowski/Langmuir equation.",
          {"concentration": _num("surfactant concentration (same unit as 1/K)"),
           "gamma0_mN_m": _num("pure-solvent surface tension (mN/m)"),
           "gamma_max_mol_per_m2": _num("saturation surface excess (mol/m^2)"),
           "K": _num("Szyszkowski/Langmuir adsorption constant"),
           "temperature_K": _num("temperature in Kelvin"),
           "n_factor": _num("Gibbs prefactor (1 or 2)")},
          ["concentration", "gamma0_mN_m", "gamma_max_mol_per_m2", "K"], _szysz),
    _spec("hlb_from_mw",
          "HLB by Griffin's method (0-20 scale) from hydrophilic and total molecular weight.",
          {"mw_hydrophilic": _num("hydrophilic-portion MW (g/mol)"),
           "mw_total": _num("total MW (g/mol)")},
          ["mw_hydrophilic", "mw_total"], _hlb_mw),
    _spec("hlb_from_groups",
          "HLB by Davies' group-contribution method. Raises if a group has no verified Davies number.",
          {"group_counts": {"type": "object",
                            "description": "map of Davies group name -> count, e.g. {\"SO4Na\":1,\"CH2\":11,\"CH3\":1}",
                            "additionalProperties": {"type": "integer"}}},
          ["group_counts"], _hlb_groups),
    _spec("tanford_chain_geometry",
          "Tanford tail volume (cubic Angstrom) and critical/extended chain length (Angstrom) for "
          "a saturated unbranched alkyl chain.",
          {"n_carbons": {"type": "integer", "description": "number of carbons in the tail"}},
          ["n_carbons"], _tanford),
    _spec("critical_packing_parameter",
          "Critical packing parameter CPP = v/(a0*lc) and predicted aggregate morphology.",
          {"volume_A3": _num("tail volume (cubic Angstrom)"),
           "head_area_A2": _num("optimal headgroup area a0 (square Angstrom)"),
           "length_A": _num("critical chain length lc (Angstrom)")},
          ["volume_A3", "head_area_A2", "length_A"], _cpp),
    _spec("aggregation_number",
          "Geometric aggregation number for a spherical micelle from tail volume and core radius.",
          {"tail_volume_A3": _num("tail volume per surfactant (cubic Angstrom)"),
           "core_radius_A": _num("micelle core radius (Angstrom)")},
          ["tail_volume_A3", "core_radius_A"], _nagg),
    _spec("debye_screening_length",
          "Debye screening length (nm). Provide either ionic_strength_M directly, or an "
          "ion_concentrations map {charge_number: molar_concentration} to compute I first.",
          {"ionic_strength_M": _num("ionic strength (mol/L), if known directly"),
           "ion_concentrations": {"type": "object",
                                  "description": "map of signed charge number -> concentration (mol/L), "
                                  "e.g. {\"1\":0.1,\"-1\":0.1} for 0.1 M NaCl",
                                  "additionalProperties": {"type": "number"}},
           "temperature_K": _num("temperature in Kelvin (default 298.15)")},
          [], _debye),
    _spec("zeta_potential",
          "Zeta potential (mV) via the Henry equation from electrophoretic mobility.",
          {"electrophoretic_mobility_um_cm_per_Vs": _num("mobility in (um*cm)/(V*s)"),
           "viscosity_mPas": _num("solvent viscosity (mPa.s)"),
           "regime": {"type": "string", "enum": ["huckel", "smoluchowski"],
                      "description": "huckel (small particles/low I) or smoluchowski (large/high I)"}},
          ["electrophoretic_mobility_um_cm_per_Vs", "viscosity_mPas", "regime"], _zeta),
    _spec("hydrodynamic_radius",
          "Hydrodynamic radius (nm) via Stokes-Einstein from a diffusion coefficient.",
          {"diffusion_coefficient_cm2_per_s": _num("translational D in cm^2/s (NOT m^2/s)"),
           "viscosity_mPas": _num("solvent viscosity (mPa.s)"),
           "temperature_K": _num("temperature in Kelvin (default 298.15)")},
          ["diffusion_coefficient_cm2_per_s", "viscosity_mPas"], _rh),
    _spec("counterion_binding_degree",
          "Degree of counterion binding beta from the conductometric slope-ratio method.",
          {"slope_below_cmc": _num("conductivity-vs-concentration slope below the CMC"),
           "slope_above_cmc": _num("slope above the CMC (must be smaller)")},
          ["slope_below_cmc", "slope_above_cmc"], _beta_bind),
    _spec("gibbs_free_energy_micellization",
          "Standard Gibbs free energy of micellization deltaG_mic (kJ/mol). Converts CMC (mol/L) "
          "to mole fraction internally.",
          {"cmc_M": _num("CMC in mol/L"),
           "temperature_K": _num("temperature in Kelvin"),
           "counterion_factor": _num("1.0 nonionic; (2-beta) for an ionic surfactant")},
          ["cmc_M", "temperature_K"], _dg_mic),
    _spec("vant_hoff_enthalpy",
          "Enthalpy of micellization deltaH_mic (kJ/mol) via van't Hoff, from CMC at two temperatures.",
          {"cmc1_M": _num("CMC at temperature1 (mol/L)"),
           "temperature1_K": _num("temperature 1 (K)"),
           "cmc2_M": _num("CMC at temperature2 (mol/L)"),
           "temperature2_K": _num("temperature 2 (K)")},
          ["cmc1_M", "temperature1_K", "cmc2_M", "temperature2_K"], _dh_mic),
    _spec("entropy_of_micellization",
          "Entropy of micellization deltaS_mic (J/mol.K) from deltaG, deltaH, T.",
          {"deltaG_mic_kJ_per_mol": _num("deltaG_mic (kJ/mol)"),
           "deltaH_mic_kJ_per_mol": _num("deltaH_mic (kJ/mol)"),
           "temperature_K": _num("temperature in Kelvin")},
          ["deltaG_mic_kJ_per_mol", "deltaH_mic_kJ_per_mol", "temperature_K"], _ds_mic),
    _spec("wetting_work_of_adhesion",
          "Young-Dupre work of adhesion (mJ/m^2) from surface tension and contact angle.",
          {"gamma_LV_mN_m": _num("liquid-vapor surface tension (mN/m)"),
           "contact_angle_deg": _num("contact angle (degrees)")},
          ["gamma_LV_mN_m", "contact_angle_deg"], _wa),
    _spec("wetting_spreading_coefficient",
          "Spreading coefficient S (mN/m) from surface tension and contact angle.",
          {"gamma_LV_mN_m": _num("liquid-vapor surface tension (mN/m)"),
           "contact_angle_deg": _num("contact angle (degrees)")},
          ["gamma_LV_mN_m", "contact_angle_deg"], _spread),
    _spec("eor_capillary_number",
          "Capillary number Ca (dimensionless) and flow regime for enhanced oil recovery.",
          {"viscosity_mPas": _num("displacing-fluid viscosity (mPa.s)"),
           "velocity_m_per_s": _num("velocity (m/s)"),
           "interfacial_tension_mN_m": _num("oil-water interfacial tension (mN/m)")},
          ["viscosity_mPas", "velocity_m_per_s", "interfacial_tension_mN_m"], _ca),
    _spec("molar_solubilization_ratio",
          "Molar Solubilization Ratio (dimensionless): moles solubilizate per mole micellized "
          "surfactant. Requires surfactant concentration above the CMC.",
          {"total_solubilized_M": _num("total solubilized solute (mol/L)"),
           "intrinsic_water_solubility_M": _num("intrinsic water solubility of the solute (mol/L)"),
           "surfactant_concentration_M": _num("surfactant concentration (mol/L), must exceed CMC"),
           "cmc_M": _num("CMC (mol/L)")},
          ["total_solubilized_M", "intrinsic_water_solubility_M", "surfactant_concentration_M", "cmc_M"], _msr),
]

TOOLS = {t["name"]: t for t in _TOOL_LIST}
assert len(TOOLS) == 25, f"expected 25 tools, got {len(TOOLS)}"


def dispatch(name: str, args: dict) -> dict:
    """Execute a tool call against the real library. Returns a JSON-serialisable
    dict; on any library error returns {'error': <message>} rather than raising,
    so the calling model can react to it the way a real tool user would."""
    if name not in TOOLS:
        return {"error": f"unknown tool '{name}'"}
    try:
        return TOOLS[name]["fn"](args)
    except Exception as e:  # library raised (bad input, no-solution guard, etc.)
        return {"error": f"{type(e).__name__}: {e}"}


# --- per-provider schema converters -----------------------------------------

def to_openai_tools() -> list[dict]:
    return [{"type": "function",
             "function": {"name": t["name"], "description": t["description"],
                          "parameters": t["parameters"]}}
            for t in _TOOL_LIST]

def to_anthropic_tools() -> list[dict]:
    return [{"name": t["name"], "description": t["description"],
             "input_schema": t["parameters"]}
            for t in _TOOL_LIST]

def to_gemini_tools() -> list[dict]:
    # Gemini function declarations; it does not accept additionalProperties, so strip it.
    def clean(schema):
        props = {}
        for k, v in schema["properties"].items():
            v2 = {kk: vv for kk, vv in v.items() if kk != "additionalProperties"}
            props[k] = v2
        return {"type": "object", "properties": props, "required": schema.get("required", [])}
    return [{"function_declarations": [
        {"name": t["name"], "description": t["description"], "parameters": clean(t["parameters"])}
        for t in _TOOL_LIST]}]


if __name__ == "__main__":
    # quick self-check: every tool dispatches on a representative input
    import json
    print(f"{len(TOOLS)} tools registered")
    samples = {
        "clint_ideal_cmc": {"alpha1": 0.5, "cmc1_mM": 1.0, "cmc2_mM": 5.0},
        "rubingh_solve": {"alpha1": 0.5, "cmc_mix_mM": 0.8, "cmc1_mM": 1.0, "cmc2_mM": 5.0},
        "debye_screening_length": {"ionic_strength_M": 0.1},
        "hlb_from_groups": {"group_counts": {"SO4Na": 1, "CH2": 11, "CH3": 1}},
        "tanford_chain_geometry": {"n_carbons": 12},
    }
    for name, args in samples.items():
        print(f"  {name}: {json.dumps(dispatch(name, args))}")
