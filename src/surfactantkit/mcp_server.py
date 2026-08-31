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
from . import electrostatics as elec
from . import dynamics as dyn
from . import thermodynamics as thermo
from . import wetting
from . import solubilization as solub

mcp = MCPServer(
    "SurfactantKit",
    instructions=(
        "Surfactant and interfacial-science calculations: Clint ideal mixing, "
        "Rubingh regular-solution theory for binary surfactant mixtures (and its "
        "Rosen extension to mixed adsorbed monolayers), the Gibbs adsorption "
        "isotherm, HLB (Griffin's and Davies' methods), the critical packing "
        "parameter and aggregation number (Tanford's formulas), full "
        "micellization thermodynamics (deltaG/deltaH/deltaS, counterion "
        "binding), the Corrin-Harkins salt-CMC relation, Debye screening "
        "length, zeta potential (Henry equation), and hydrodynamic radius "
        "(Stokes-Einstein). Use these tools instead of computing or recalling "
        "these values from memory -- CMC, interaction-parameter, HLB, and "
        "electrolyte/electrostatic values are easy to hallucinate and easy to "
        "get wrong, especially for gemini (two-headgroup) surfactants, mixed "
        "salt/ionic-strength systems, and anything requiring iterative "
        "numerical solving or careful unit conversion (e.g. cP vs Pa.s, "
        "cm^2/s vs m^2/s -- these are exactly the kind of silent unit errors "
        "this server exists to prevent). Also covers the Szyszkowski surface "
        "tension equation, wetting/adhesion (Young-Dupre, spreading "
        "coefficient), the capillary number (enhanced oil recovery), and "
        "micellar solubilization (Molar Solubilization Ratio)."
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
def szyszkowski_predict_surface_tension(concentration: float, gamma0_mN_m: float, gamma_max_mol_per_m2: float, K: float, temperature_K: float = 298.15, n_factor: float = 1.0) -> dict:
    """Predict surface tension (mN/m) at a given surfactant concentration
    via the Szyszkowski/Langmuir equation: gamma = gamma0 -
    n_factor*R*T*Gamma_max*ln(1+K*concentration). K is a fitted
    Szyszkowski/Langmuir constant specific to the surfactant + condition
    (not a universal constant -- must come from real data, e.g. a fit to
    measured surface-tension-vs-concentration points). concentration and
    K must use consistent units (K*concentration must be dimensionless).
    Only valid below the CMC, where surfactant is present as free
    monomer at the interface."""
    value = ads.szyszkowski_surface_tension(concentration, gamma0_mN_m, gamma_max_mol_per_m2, K, temperature_K, n_factor)
    return {"surface_tension_mN_m": value, "unit": "mN/m"}


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


@mcp.tool()
def aggregation_number(tail_volume_A3: float, core_radius_A: float) -> dict:
    """Estimated spherical-micelle aggregation number from geometric
    packing: N_agg = core_volume / tail_volume. tail_volume_A3 typically
    comes from tanford_chain_geometry; core_radius_A is usually
    approximated by the critical chain length from that same tool, but
    pass a measured value (e.g. from SANS/SAXS) if you have one. This is
    a geometric estimate, not a substitute for a directly measured
    aggregation number (e.g. fluorescence quenching)."""
    value = cpp_mod.aggregation_number_spherical(tail_volume_A3, core_radius_A)
    return {"aggregation_number": value, "unit": "dimensionless (molecules per micelle)",
            "note": "geometric estimate assuming a spherical micelle core"}


@mcp.tool()
def rosen_monolayer_solve(alpha1: float, c_mix_sigma_mM: float, c1_sigma_mM: float, c2_sigma_mM: float) -> dict:
    """Rosen's extension of Rubingh theory to the mixed ADSORBED
    MONOLAYER at an interface, from surface-tension data -- NOT CMC
    data. c1_sigma_mM/c2_sigma_mM are the pure-component concentrations
    needed to reach a chosen reference surface tension (or pressure);
    c_mix_sigma_mM is the mixed concentration needed to reach that SAME
    reference surface tension at bulk mole fraction alpha1. Using CMC
    values here instead is a common category error -- if you have CMC
    data, use rubingh_solve instead."""
    x1 = mm.solve_rosen_monolayer_x(alpha1, c_mix_sigma_mM, c1_sigma_mM, c2_sigma_mM)
    if x1 is None:
        raise ValueError("No valid root found in (0, 1) for these inputs.")
    beta_sigma = mm.rosen_beta_sigma(x1, alpha1, c_mix_sigma_mM, c1_sigma_mM)
    return {
        "monolayer_mole_fraction_x1": x1,
        "beta_sigma": beta_sigma,
        "synergy_classification": "synergistic" if beta_sigma < 0 else ("antagonistic" if beta_sigma > 0 else "ideal"),
    }


@mcp.tool()
def corrin_harkins_predict(cmc1_mM: float, salt_conc1_mM: float, cmc2_mM: float, salt_conc2_mM: float, salt_conc_target_mM: float) -> dict:
    """Predict CMC at a new counterion (salt) concentration via the
    Corrin-Harkins log-linear relation, fit from two known (CMC, salt
    concentration) data points for the SAME surfactant + salt system.
    The fitted slope g and intercept are system-specific -- there is no
    universal lookup table for them, which is why two real data points
    are required as input. Only reliable as interpolation between (or a
    close extrapolation beyond) the two given concentrations."""
    predicted_cmc, g = mm.corrin_harkins_predict_cmc(cmc1_mM, salt_conc1_mM, cmc2_mM, salt_conc2_mM, salt_conc_target_mM)
    return {"predicted_cmc_mM": predicted_cmc, "fitted_slope_g": g}


@mcp.tool()
def debye_screening_length(ionic_strength_M: float, temperature_K: float = 298.15) -> dict:
    """Debye screening length (nm) for an electrolyte solution.
    ionic_strength_M must be the IONIC STRENGTH (0.5*sum(c_i*z_i^2)),
    not the raw salt concentration -- for a simple 1:1 salt like NaCl
    these are numerically equal, but for anything else (CaCl2, MgSO4,
    etc.) they are not. Use the ionic_strength helper concept: for each
    ion, multiply its molar concentration by its charge squared, sum,
    and halve."""
    value = elec.debye_length(ionic_strength_M, temperature_K)
    return {"debye_length_nm": value, "unit": "nm"}


@mcp.tool()
def zeta_potential(electrophoretic_mobility_um_cm_per_Vs: float, viscosity_mPas: float, regime: str) -> dict:
    """Zeta potential (mV) from electrophoretic mobility via the Henry
    equation. regime must be explicitly stated as 'huckel' (small
    particles / low ionic strength, kappa*a << 1) or 'smoluchowski'
    (larger colloids, ionic strength >= ~10 mM, kappa*a >> 1) -- there is
    no default, because silently picking the wrong regime is exactly the
    kind of error this tool exists to prevent. viscosity_mPas: solvent
    viscosity in mPa.s (=cP; water is ~0.89 at 25C -- do not pass Pa.s
    here, a 1000x unit error)."""
    value = elec.zeta_potential_henry(electrophoretic_mobility_um_cm_per_Vs, viscosity_mPas, regime)
    return {"zeta_potential_mV": value, "unit": "mV", "regime_used": regime.lower()}


@mcp.tool()
def hydrodynamic_radius(diffusion_coefficient_cm2_per_s: float, viscosity_mPas: float, temperature_K: float = 298.15) -> dict:
    """Hydrodynamic radius (nm) from a DLS-measured diffusion
    coefficient via the Stokes-Einstein equation. diffusion_coefficient
    must be in cm^2/s (the commonly-reported DLS unit) -- NOT m^2/s, a
    classic 1e4x unit error. viscosity_mPas: solvent viscosity in mPa.s
    (=cP; water is ~0.89 at 25C) -- NOT Pa.s."""
    value = dyn.hydrodynamic_radius_stokes_einstein(diffusion_coefficient_cm2_per_s, viscosity_mPas, temperature_K)
    return {"hydrodynamic_radius_nm": value, "unit": "nm"}


@mcp.tool()
def counterion_binding_degree(slope_below_cmc: float, slope_above_cmc: float) -> dict:
    """Degree of counterion binding to the micelle (beta) from the
    slope-ratio method on a conductivity-vs-concentration plot: beta = 1
    - (slope_above_cmc / slope_below_cmc). Both slopes must be positive,
    with slope_above_cmc < slope_below_cmc. Note: this method's physical
    interpretation has been questioned for some systems in the
    literature -- the arithmetic is standard and widely used, but treat
    the result as an approximate, commonly-reported value."""
    value = thermo.counterion_binding_degree(slope_below_cmc, slope_above_cmc)
    return {"beta": value, "unit": "dimensionless (fraction, 0 to 1)"}


@mcp.tool()
def gibbs_free_energy_micellization(cmc_M: float, temperature_K: float, counterion_factor: float = 1.0) -> dict:
    """Standard Gibbs free energy of micellization (kJ/mol) from CMC
    (molarity, mol/L) and temperature. counterion_factor must be stated
    explicitly: use 1.0 for a nonionic surfactant, or (2 - beta) for an
    ionic surfactant given its counterion binding degree beta (see
    counterion_binding_degree) -- this is never auto-detected, since
    silently assuming nonionic behavior for an ionic surfactant (or vice
    versa) is exactly the kind of hidden-assumption error this tool
    exists to prevent."""
    x_cmc = thermo.cmc_to_mole_fraction(cmc_M)
    value = thermo.gibbs_free_energy_micellization(x_cmc, temperature_K, counterion_factor)
    return {"deltaG_mic_kJ_per_mol": value, "unit": "kJ/mol", "cmc_mole_fraction_used": x_cmc}


@mcp.tool()
def vant_hoff_enthalpy(cmc1_M: float, temperature1_K: float, cmc2_M: float, temperature2_K: float) -> dict:
    """Van't Hoff enthalpy of micellization (kJ/mol) from CMC (molarity)
    measured at two temperatures. Assumes deltaH is constant over the
    temperature interval -- treat results from a wide T range with
    appropriate caution, and note this breaks down if the aggregation
    number itself varies significantly with temperature."""
    x1 = thermo.cmc_to_mole_fraction(cmc1_M)
    x2 = thermo.cmc_to_mole_fraction(cmc2_M)
    value = thermo.vant_hoff_enthalpy(x1, temperature1_K, x2, temperature2_K)
    return {"deltaH_mic_kJ_per_mol": value, "unit": "kJ/mol"}


@mcp.tool()
def entropy_of_micellization(deltaG_mic_kJ_per_mol: float, deltaH_mic_kJ_per_mol: float, temperature_K: float) -> dict:
    """Entropy of micellization (J/(mol.K)): deltaS = (deltaH - deltaG) /
    T. Completes the deltaG/deltaH/deltaS thermodynamic triad given the
    other two (from gibbs_free_energy_micellization and
    vant_hoff_enthalpy)."""
    value = thermo.entropy_micellization(deltaG_mic_kJ_per_mol, deltaH_mic_kJ_per_mol, temperature_K)
    return {"deltaS_mic_J_per_mol_K": value, "unit": "J/(mol.K)"}


@mcp.tool()
def wetting_work_of_adhesion(gamma_LV_mN_m: float, contact_angle_deg: float) -> dict:
    """Young-Dupre work of adhesion (mJ/m^2, numerically = mN/m):
    W_a = gamma_LV * (1 + cos(theta)). gamma_LV_mN_m is the liquid-vapor
    surface tension; contact_angle_deg is the measured equilibrium
    contact angle in degrees (0-180)."""
    value = wetting.work_of_adhesion(gamma_LV_mN_m, contact_angle_deg)
    return {"work_of_adhesion_mJ_per_m2": value, "unit": "mJ/m^2"}


@mcp.tool()
def wetting_spreading_coefficient(gamma_LV_mN_m: float, contact_angle_deg: float) -> dict:
    """Spreading coefficient (mN/m) from Young's equation:
    S = gamma_LV * (cos(theta) - 1). Always <= 0 in this formulation
    (equality only at complete wetting, theta=0); more negative means
    poorer spreading."""
    value = wetting.spreading_coefficient(gamma_LV_mN_m, contact_angle_deg)
    return {"spreading_coefficient_mN_m": value, "unit": "mN/m"}


@mcp.tool()
def eor_capillary_number(viscosity_mPas: float, velocity_m_per_s: float, interfacial_tension_mN_m: float) -> dict:
    """Capillary number (dimensionless): Ca = (viscosity*velocity) /
    interfacial_tension -- governs oil-displacement efficiency in
    enhanced oil recovery. Below Ca ~ 1e-5, flow is capillary-dominated
    and residual oil stays trapped; surfactants raise Ca mainly by
    driving interfacial tension toward ultra-low values. viscosity_mPas
    in mPa.s (=cP); velocity_m_per_s in m/s; interfacial_tension_mN_m in mN/m."""
    value = wetting.capillary_number(viscosity_mPas, velocity_m_per_s, interfacial_tension_mN_m)
    return {"capillary_number": value, "unit": "dimensionless",
            "flow_regime": "capillary-dominated (residual oil trapped)" if value < 1e-5 else "transitional/viscous-dominated"}


@mcp.tool()
def molar_solubilization_ratio(total_solubilized_M: float, intrinsic_water_solubility_M: float, surfactant_concentration_M: float, cmc_M: float) -> dict:
    """Molar Solubilization Ratio (MSR, dimensionless): moles of
    solubilizate taken up per mole of MICELLIZED surfactant (only the
    surfactant above the CMC, since that's what forms micelles):
    MSR = (total_solubilized - intrinsic_water_solubility) /
    (surfactant_concentration - cmc). All four inputs must be in the
    same concentration unit. Note: the companion quantity Km
    (micelle-water partition coefficient) is deliberately NOT provided
    by this server -- its exact formula is study-specific in the
    literature with no single dominant convention, and guessing one
    risks reproducing the wrong definition for a given paper."""
    value = solub.molar_solubilization_ratio(total_solubilized_M, intrinsic_water_solubility_M, surfactant_concentration_M, cmc_M)
    return {"msr": value, "unit": "dimensionless (mol solubilizate per mol micellized surfactant)"}


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
