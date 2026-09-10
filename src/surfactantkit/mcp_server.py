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

from mcp.server.fastmcp import FastMCP

from . import mixed_micelle as mm
from . import adsorption as ads
from . import hlb as hlb_mod
from . import hld as hld_mod
from . import cpp as cpp_mod
from . import electrostatics as elec
from . import dynamics as dyn
from . import thermodynamics as thermo
from . import wetting
from . import solubilization as solub
from . import curve_analysis as curve

mcp = FastMCP(
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
            "literature_validation_notes.md in the SurfactantKit repo. For a "
            "whole composition series, use rubingh_beta_regression instead."
        ),
    }


@mcp.tool()
def rubingh_beta_regression(alpha1_series: list[float], cmc_mix_series_mM: list[float], cmc1_mM: float, cmc2_mM: float) -> dict:
    """Fit a single representative Rubingh beta across a WHOLE composition
    series -- the real standard practice for reporting one interaction
    parameter for a binary system, rather than rubingh_solve's single-point
    value. Computes the pointwise x1/beta independently at each
    (alpha1, cmc_mix) point in the series (same math as rubingh_solve),
    then reports the arithmetic mean beta plus its standard deviation
    across the series (a large spread signals the system does not
    actually have one constant beta, i.e. regular-solution theory is only
    an approximation for it -- a real diagnostic, not noise). Points
    where no valid root exists (near-ideal or self-inconsistent data) are
    skipped, not errors, unless fewer than 2 of the series solve at all.
    alpha1_series and cmc_mix_series_mM must be the same length and in
    matching order; cmc1_mM/cmc2_mM are the pure-component CMCs."""
    result = mm.rubingh_beta_regression(alpha1_series, cmc_mix_series_mM, cmc1_mM, cmc2_mM)
    synergy = "synergistic" if result.beta_mean < 0 else ("antagonistic" if result.beta_mean > 0 else "ideal")
    return {
        "beta_mean": result.beta_mean,
        "beta_std": result.beta_std,
        "beta_unit": "dimensionless (regular solution interaction parameter)",
        "synergy_classification": synergy,
        "betas_per_point": result.betas,
        "x1_per_point": result.x1_values,
        "alpha1_used": result.alpha1_used,
        "n_points_used": result.n_points_used,
        "n_points_skipped": result.n_points_skipped,
    }


@mcp.tool()
def motomura_ideal_composition(alpha1: float, cmc1_mM: float, cmc2_mM: float) -> dict:
    """Motomura & Aratono's ideal micellar composition (X1) for a binary
    surfactant mixture -- the zero-interaction-assumption prediction for
    how much of component 1 ends up IN the micelle, as opposed to
    clint_ideal_cmc which predicts the ideal mixed CMC itself. Same
    zero-parameter ideal-mixing assumption as Clint's model (no
    interaction parameter, unlike rubingh_solve's beta); use this when
    the ideal-mixing micellar composition itself is asked for, not just
    the ideal CMC. alpha1 is the bulk mole fraction of component 1;
    cmc1_mM/cmc2_mM are the pure-component CMCs, in mM."""
    value = mm.motomura_ideal_composition(alpha1, cmc1_mM, cmc2_mM)
    return {"x1_ideal": value, "unit": "dimensionless (micellar mole fraction)", "model": "Motomura ideal mixing"}


@mcp.tool()
def rodenas_x1(alpha1: float, dln_cmc_mix_dalpha1: float) -> dict:
    """Rodenas' Gibbs-Duhem-based, "model-independent" micellar mole
    fraction of component 1 -- an alternative to rubingh_solve that
    does NOT assume any particular regular-solution/activity-
    coefficient model, at the cost of needing the LOCAL SLOPE of the
    real, multi-point experimental ln(cmc_mix) vs. alpha1 curve
    (dln_cmc_mix_dalpha1) rather than a single (alpha1, cmc_mix) pair.
    This slope MUST come from real neighboring experimental data
    (a fit or numerical derivative) -- do not guess or assume a value.
    Follow with rodenas_activity_coefficients using the returned x1."""
    value = mm.rodenas_x1(alpha1, dln_cmc_mix_dalpha1)
    return {"micellar_mole_fraction_x1": value, "model": "Rodenas (Gibbs-Duhem, model-independent)"}


@mcp.tool()
def rodenas_x1_from_series(alpha1_series: list[float], cmc_mix_series_mM: list[float]) -> dict:
    """Compute Rodenas' d[ln(cmc_mix)]/d(alpha1) AND the resulting x1 at
    every composition directly from a RAW multi-alpha1 cmc_mix(alpha1)
    series -- closes the gap rodenas_x1 itself flags (it needs the local
    slope as an already-known input). At each composition, fits the
    unique quadratic through it and its 2 nearest neighbors along
    alpha1, then differentiates analytically at that exact alpha1 --
    this reduces exactly to the standard central/forward/backward finite-
    difference formulas for evenly-spaced data, generalized to real
    (unevenly-spaced) experimental series. alpha1_series and
    cmc_mix_series_mM must be the same length, at least 3 points, all
    alpha1 distinct and strictly between 0 and 1. Points are sorted by
    alpha1 internally; results are returned in that sorted order."""
    result = mm.rodenas_x1_series(alpha1_series, cmc_mix_series_mM)
    return {
        "alpha1_used": result.alpha1_used,
        "dln_cmc_mix_dalpha1": result.dln_cmc_mix_dalpha1,
        "x1_rodenas": result.x1_rodenas,
        "n_points": result.n_points,
        "method": result.method,
    }


@mcp.tool()
def rodenas_activity_coefficients(alpha1: float, x1: float, cmc_mix_mM: float, cmc1_mM: float, cmc2_mM: float) -> dict:
    """Activity coefficients (f1, f2) for Rodenas' model, given the
    already-determined x1 (from rodenas_x1). Same mass-balance
    definition solve_rubingh_x uses internally -- Rodenas and Rubingh
    differ only in how x1 is determined, not in how f1/f2 follow from
    a known x1."""
    f1, f2 = mm.rodenas_activity_coefficients(alpha1, x1, cmc_mix_mM, cmc1_mM, cmc2_mM)
    return {"f1": f1, "f2": f2, "unit": "dimensionless", "model": "Rodenas"}


@mcp.tool()
def rubingh_activity_coefficients(x1: float, beta: float) -> dict:
    """Regular-solution activity coefficients (f1, f2) for a binary
    mixed micelle, given the micellar mole fraction x1 and interaction
    parameter beta (from rubingh_solve)."""
    f1, f2 = mm.activity_coefficients(x1, beta)
    return {"f1": f1, "f2": f2, "unit": "dimensionless"}


@mcp.tool()
def asymmetric_margules_activity_coefficients(x1: float, w12: float, w21: float) -> dict:
    """Activity coefficients (f1, f2) from the EOMMM binary asymmetric
    Margules model -- generalizes rubingh_activity_coefficients with
    TWO independent interaction parameters (w12 != w21) instead of
    RST's single symmetric beta. w12/w21 are dimensionless (RT units,
    same convention as rubingh_solve's beta) -- w12 is the energy of
    introducing a molecule of component 1 into a pure micelle of
    component 2 (and vice versa for w21). At w12 == w21 this reduces
    EXACTLY to rubingh_activity_coefficients (RST is the provable
    symmetric special case). ONLY computes activity coefficients for
    an ALREADY-KNOWN w12/w21 pair (e.g. quoted directly from a paper) --
    does NOT fit w12/w21 from raw CMC/alpha1 data the way rubingh_solve
    does for RST; for that, use eommm_global_fit."""
    f1, f2 = mm.asymmetric_margules_activity_coefficients(x1, w12, w21)
    return {"f1": f1, "f2": f2, "unit": "dimensionless", "model": "EOMMM asymmetric Margules (binary)"}


@mcp.tool()
def eommm_global_fit(alpha1_series: list[float], cmc_mix_series_mM: list[float], cmc1_mM: float, cmc2_mM: float) -> dict:
    """Fit EOMMM's two independent Margules interaction parameters
    (W12, W21) -- and each point's own micellar mole fraction x1 --
    SIMULTANEOUSLY across a whole composition series, closing the gap
    asymmetric_margules_activity_coefficients's own docstring flags
    (it only computes activity coefficients for an already-known
    w12/w21 pair). IMPORTANT, disclosed honestly: this is a FIRST-
    PRINCIPLES construction (least-squares over the same mass-balance
    equations rubingh_solve uses, generalized to asymmetric activity
    coefficients, solved via variable projection with multi-start
    coordinate descent), NOT a verified transcription of Schulz &
    Durand 2016's exact published procedure (their Eq. 3.2), which was
    unavailable (paywalled; see this project's own ROADMAP.md for the
    full disclosure). Validated via mathematically-guaranteed round
    trips only -- no external published raw dataset was available to
    check against. This is a genuine multi-parameter nonlinear
    optimization with no numpy/scipy dependency, so it is noticeably
    slower than other tools in this server (a few seconds for a small
    dataset) -- expect that, it is not an error. alpha1_series and
    cmc_mix_series_mM must be the same length, at least 3 points (5+
    recommended); cmc1_mM/cmc2_mM are the pure-component CMCs. Returns
    W12, W21, the fitted x1 at each point, the fit's r_squared (against
    ln(cmc_mix)), and the raw sse objective value."""
    result = mm.eommm_global_fit(alpha1_series, cmc_mix_series_mM, cmc1_mM, cmc2_mM)
    return {
        "W12": result.W12,
        "W21": result.W21,
        "x1_values": result.x1_values,
        "alpha1_used": result.alpha1_used,
        "sse": result.sse,
        "r_squared": result.r_squared,
        "n_points": result.n_points,
        "method": result.method,
    }


@mcp.tool()
def maeda_free_energy_of_micellization(x1_rub: float, beta: float, cmc1_M: float, cmc2_M: float, temperature_K: float) -> dict:
    """Maeda's free energy of micellization (kJ/mol) for a binary
    ionic/nonionic mixed micelle -- built specifically for ionic-
    nonionic systems, accounting for chain-chain (not just head-group)
    interaction. Reuses Rubingh's own already-solved x1/beta (from
    rubingh_solve) rather than fitting independently: deltaG_M =
    RT*(B0 + B1*X1_Rub + B2*X1_Rub^2), B0=ln(Xcmc2), B1+B2=ln(Xcmc1/
    Xcmc2), B2=-beta, where Xcmc are the pure-component CMCs on the
    mole-fraction scale. cmc1_M/cmc2_M must be in mol/L (molarity),
    not mM -- convert before calling."""
    value = mm.maeda_free_energy_of_micellization(x1_rub, beta, cmc1_M, cmc2_M, temperature_K)
    return {"deltaG_M_kJ_per_mol": value, "unit": "kJ/mol", "model": "Maeda"}


@mcp.tool()
def excess_free_energy(x1: float, f1: float, f2: float, temperature_K: float) -> dict:
    """Excess Gibbs free energy of mixed micelle formation (kJ/mol),
    from the micellar mole fraction x1, activity coefficients f1/f2
    (from rubingh_activity_coefficients), and temperature in Kelvin."""
    value = mm.excess_free_energy(x1, f1, f2, temperature_K)
    return {"deltaG_ex_kJ_per_mol": value, "unit": "kJ/mol"}


@mcp.tool()
def gibbs_surface_excess(slope_mN_per_m_per_lnC: float, system_type: str, temperature_K: float) -> dict:
    """Maximum surface excess concentration (Gamma_max, mol/m^2) from
    the Gibbs adsorption isotherm, given the pre-CMC slope of surface
    tension (mN/m) vs ln(concentration), temperature in Kelvin, and
    system_type -- one of 'nonionic', 'ionic_excess_electrolyte', or
    'ionic_no_added_salt' -- which fixes the Gibbs prefactor n (1 for
    the first two, 2 for the last). REQUIRED, no default: if the
    question does not state the surfactant's ionic character AND
    (for an ionic surfactant) whether excess inert electrolyte is
    present, DO NOT GUESS a system_type -- this is exactly the kind of
    silent assumption that produces a real, confidently wrong answer;
    report that Gamma_max cannot be determined from the given
    information instead of calling this tool with an assumed value."""
    value = ads.gibbs_gamma_max(slope_mN_per_m_per_lnC, system_type, temperature_K)
    return {"gamma_max_mol_per_m2": value, "unit": "mol/m^2"}


@mcp.tool()
def gibbs_area_per_molecule(gamma_max_mol_per_m2: float) -> dict:
    """Minimum area per molecule at the interface (A_min, nm^2), from
    Gamma_max (mol/m^2, from gibbs_surface_excess)."""
    value = ads.gibbs_a_min(gamma_max_mol_per_m2)
    return {"a_min_nm2": value, "unit": "nm^2 per molecule"}


@mcp.tool()
def szyszkowski_predict_surface_tension(concentration: float, gamma0_mN_m: float, gamma_max_mol_per_m2: float, K: float, system_type: str, temperature_K: float = 298.15) -> dict:
    """Predict surface tension (mN/m) at a given surfactant concentration
    via the Szyszkowski/Langmuir equation: gamma = gamma0 -
    n_factor*R*T*Gamma_max*ln(1+K*concentration). K is a fitted
    Szyszkowski/Langmuir constant specific to the surfactant + condition
    (not a universal constant -- must come from real data, e.g. a fit to
    measured surface-tension-vs-concentration points). concentration and
    K must use consistent units (K*concentration must be dimensionless).
    system_type: same required, no-default Gibbs-prefactor selector as
    gibbs_surface_excess ('nonionic', 'ionic_excess_electrolyte', or
    'ionic_no_added_salt') -- DO NOT GUESS this if the surfactant's ionic
    character or electrolyte condition isn't stated; report that surface
    tension cannot be predicted instead. Only valid below the CMC, where
    surfactant is present as free monomer at the interface."""
    value = ads.szyszkowski_surface_tension(concentration, gamma0_mN_m, gamma_max_mol_per_m2, K, system_type, temperature_K)
    return {"surface_tension_mN_m": value, "unit": "mN/m"}


@mcp.tool()
def szyszkowski_fit_K(concentrations: list[float], surface_tensions_mN_per_m: list[float], gamma0_mN_m: float, gamma_max_mol_per_m2: float, system_type: str, temperature_K: float = 298.15) -> dict:
    """Fit the Szyszkowski/Langmuir adsorption constant K from a RAW
    PRE-CMC surface-tension-vs-concentration curve -- the fitting step
    szyszkowski_predict_surface_tension's own K parameter otherwise
    assumes is already known. Real, standard two-step workflow:
    gamma_max_mol_per_m2 must come independently from the SAME curve's
    pre-CMC log-slope (gibbs_surface_excess); this tool then does a
    genuine nonlinear least-squares fit for K alone against the
    Szyszkowski equation. concentrations and surface_tensions_mN_per_m
    must be PRE-CMC data only (points at/above the CMC violate the
    model's own assumption and will bias the fit) and the same length,
    at least 3 points. system_type: same required, no-default Gibbs-
    prefactor selector as gibbs_surface_excess -- do not guess. Returns
    the fitted K, the fit's r_squared (a real diagnostic: well below 1
    means this curve does not actually follow Szyszkowski/Langmuir
    behavior), and the model-predicted surface tension at each input
    concentration."""
    result = ads.szyszkowski_fit_K(concentrations, surface_tensions_mN_per_m, gamma0_mN_m, gamma_max_mol_per_m2, system_type, temperature_K)
    return {
        "K": result.K,
        "r_squared": result.r_squared,
        "n_points": result.n_points,
        "fitted_surface_tensions_mN_per_m": result.fitted_surface_tensions_mN_per_m,
    }


@mcp.tool()
def frumkin_predict_surface_tension(concentration: float, gamma0_mN_m: float, gamma_max_mol_per_m2: float, K: float, a: float, system_type: str, temperature_K: float = 298.15) -> dict:
    """Predict surface tension (mN/m) via the Frumkin isotherm, which
    generalizes szyszkowski_predict_surface_tension with a dimensionless
    lateral-interaction parameter a between adsorbed molecules: gamma =
    gamma0 + n_factor*R*T*Gamma_max*[ln(1-theta) + a*theta^2], theta
    solved from K*concentration = theta/(1-theta)*exp(-2*a*theta).
    a > 0 = net ATTRACTIVE interaction (steeper isotherm, more favorable
    packing); a < 0 = net REPULSIVE interaction; a = 0 reduces EXACTLY
    to szyszkowski_predict_surface_tension -- use Szyszkowski (a=0) when
    there is no reason to believe lateral interactions matter, and only
    supply a non-zero a when it is fit from real surface-tension-vs-
    concentration data (never guess a value). For |a| >= 2 the Frumkin
    isotherm is genuinely multi-valued (a real 2D-condensation/phase-
    transition effect) -- a fit landing there should be flagged as near
    a condensation transition, not treated as an ordinary single-valued
    isotherm. K, gamma_max_mol_per_m2, system_type, temperature_K: same
    meaning and unit conventions as szyszkowski_predict_surface_tension;
    system_type is required, no default, do not guess it."""
    value = ads.frumkin_surface_tension(concentration, gamma0_mN_m, gamma_max_mol_per_m2, K, a, system_type, temperature_K)
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
    method (extended with a real, source-verified quaternary ammonium
    number -- see below): HLB = 7 + sum(hydrophilic group numbers) -
    sum(lipophilic group numbers). group_counts keys must be from the
    verified group list: hydrophilic {SO4Na, COOK, COONa,
    N_tertiary_amine, N_quat_trimethyl, N_quat_dimethyl_dialkyl,
    ester_sorbitan_ring, ester_free, COOH, OH_free, O_ether,
    OH_sorbitan_ring}; lipophilic {CH2, CH3, CH, vinyl_CH}.
    N_quat_trimethyl = -N+(CH3)3, a single-tail quaternary ammonium head
    (e.g. CTAB, LTAB, DTAB); N_quat_dimethyl_dialkyl = >N+(CH3)2, a
    double-tail quaternary ammonium head (e.g. DDAB, DODAB) -- both
    assume a Cl-/Br- counterion (not verified for others, e.g. tosylate).
    Raises an error for any other group name rather than guessing a
    value -- most notably, amide and sulfonate groups still have NO
    verified number in this table; do not estimate one."""
    value = hlb_mod.hlb_davies(group_counts)
    return {"hlb": value, "scale": "Davies (7 = neutral reference point)"}


@mcp.tool()
def hlb_davies_guo_ecl(n_carbons_alkyl: int, n_eo: int, extra_group_counts: dict[str, int] | None = None) -> dict:
    """HLB for a NONIONIC polyethoxylated surfactant (e.g. a fatty
    alcohol ethoxylate, R-(OCH2CH2)n-OH) via Guo, Rong & Ying's 2006
    refinement of Davies' method -- a real, sourced ALTERNATIVE to
    hlb_from_groups specifically for this surfactant class: uses
    "effective" (not actual) EO and alkyl chain lengths, which the
    source paper shows fits real HLB data substantially better (224
    surfactants tested, average absolute error < 1.5) than plugging
    actual counts into plain Davies -- because a polyoxyethylene
    chain's hydrophilic contribution is genuinely sub-linear in EO
    count, which plain Davies (linear) cannot capture. n_carbons_alkyl:
    total carbons in the hydrophobic tail. n_eo: actual EO units
    (<=50, the source paper's own valid range). extra_group_counts:
    optional additional Davies group contributions (e.g. an ester
    headgroup) not covered by the alkyl tail or EO chain. IMPORTANT,
    disclosed honestly: the primary source is paywalled -- only the
    effective-chain-length transform equations were independently
    confirmed (via a citing patent); this function assumes Davies'
    original per-unit weights (not independently confirmed from the
    primary source) apply to the effective chain lengths, a reasonable
    but unverified reading of "based on Davies' approach." Treat as a
    real, citable improvement over plain hlb_from_groups for this
    surfactant class, not an exact reproduction of the paper's own
    numbers."""
    value = hlb_mod.hlb_davies_guo_ecl(n_carbons_alkyl, n_eo, extra_group_counts)
    return {"hlb": value, "scale": "Davies (7 = neutral reference point), Guo/Rong/Ying effective-chain-length correction"}


@mcp.tool()
def hld_optimal_salinity(k: float, eacn: float, cc: float, alpha: float = 0.01, delta_T_K: float = 0.0) -> dict:
    """Optimum (balanced, Winsor III, HLD=0) salinity for an IONIC
    surfactant/oil/brine microemulsion system -- Hydrophilic-Lipophilic
    Difference (HLD) framework, S* = exp(k*EACN + alpha*deltaT - Cc).
    Use this for cationic (quaternary ammonium) or other ionic
    surfactants where Davies' HLB group-contribution method (hlb_from_groups)
    has no verified number -- HLD is a modern, physically-measured
    alternative built from real Winsor III phase-equilibrium data, not
    an additive group table. k, eacn, and cc are ALL system-specific
    and must come from real data -- do not guess. For cationic
    quaternary ammonium surfactants specifically, real measured k=0.7
    and Cc values for several common surfactants (LTAB, MTAB, CTAB,
    CMIC, BDHC, DDAB) are available -- see hld_cationic_quat_reference.
    alpha default 0.01 is the real, sourced value for ionic surfactants
    (Schirone et al. 2021); override if a different class's value is
    known. Returns salinity in equivalent grams NaCl per 100 mL."""
    value = hld_mod.optimal_salinity_ionic(k, eacn, cc, alpha, delta_T_K)
    return {"optimal_salinity_pct": value, "unit": "g NaCl equivalent / 100 mL"}


@mcp.tool()
def hld_cationic_quat_reference(surfactant: str) -> dict:
    """Real, measured HLD-framework reference values for a named
    cationic quaternary ammonium surfactant -- characteristic curvature
    (Cc), the shared k constant for this head-group class, and (where
    available) an equivalent Davies-scale HLB. surfactant must be one
    of: LTAB (dodecyltrimethylammonium bromide), MTAB
    (tetradecyltrimethylammonium bromide), CTAB
    (hexadecyltrimethylammonium bromide), CMIC (1-hexadecyl-3-
    methylimidazolium chloride), BDHC (benzyldimethylhexadecylammonium
    chloride), or DDAB (didodecyldimethylammonium bromide, double-tailed).
    Source: Schirone, Tartaro, Gentile & Palazzo, JCIS Open 4 (2021)
    100033, Table 1 -- real HLD-titration measurements, not guessed or
    estimated values. Use these Cc values with hld_optimal_salinity or
    cc_mixing_rule rather than guessing a Cc for these specific
    surfactants."""
    key = surfactant.upper()
    if key not in hld_mod.CATIONIC_QUAT_CC:
        raise ValueError(
            f"surfactant must be one of {sorted(hld_mod.CATIONIC_QUAT_CC)} -- no guessed value "
            "available for anything else"
        )
    cc, cc_uncertainty = hld_mod.CATIONIC_QUAT_CC[key]
    result = {
        "surfactant": key,
        "Cc": cc,
        "Cc_uncertainty": cc_uncertainty,
        "k": hld_mod.CATIONIC_QUAT_K,
        "k_uncertainty": hld_mod.CATIONIC_QUAT_K_UNCERTAINTY,
        "estimated_hlb_davies_scale": hld_mod.hlb_from_cc_cationic_quat(cc),
    }
    if key in hld_mod.CATIONIC_QUAT_HLB_DAVIES:
        result["real_reported_hlb_davies_scale"] = hld_mod.CATIONIC_QUAT_HLB_DAVIES[key]
    return result


@mcp.tool()
def hld_fit_k_and_cc_from_salinity_scan(eacn_series: list[float], optimal_salinity_series_pct: list[float], alpha: float = 0.01, delta_T_K: float = 0.0) -> dict:
    """Fit k and Cc from a RAW multi-oil salinity-scan series -- the real
    experimental method that produces the k/Cc inputs hld_optimal_salinity
    otherwise takes as already-known. A salinity scan measures the
    optimum (Winsor III, HLD=0) salinity S* for several different oils
    spanning a range of EACN, at fixed surfactant/temperature; ln(S*) is
    linear in EACN (ln(S*) = k*EACN + alpha*deltaT - Cc), so fitting that
    line gives k (slope) and Cc directly. eacn_series and
    optimal_salinity_series_pct (in g NaCl equivalent / 100 mL, matching
    hld_optimal_salinity's units) must be the same length, at least 2
    points (3+ recommended). Returns the fitted k, Cc, and the fit's
    r_squared -- HLD-NAC assumes k is constant across EACN for a given
    surfactant class; a poor r_squared signals that assumption doesn't
    hold for this system, not a bug in the fit."""
    result = hld_mod.fit_k_and_cc_from_salinity_scan(eacn_series, optimal_salinity_series_pct, alpha, delta_T_K)
    return {
        "k": result.k,
        "cc": result.cc,
        "r_squared": result.r_squared,
        "n_points": result.n_points,
    }


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
def nagarajan_equilibrium_area_ionic(cmc_M: float, tail_length_A: float, headgroup_prefactor_A: float, temperature_K: float = 298.15, dielectric_constant: float = 80.0) -> dict:
    """Equilibrium area per molecule a_e (square Angstrom) for an IONIC
    surfactant micelle, accounting for chain-length-dependent
    electrostatic headgroup repulsion -- Nagarajan, Langmuir 18 (2002)
    31-38. A real, sourced ALTERNATIVE to plugging a constant/assumed
    head_area_A2 straight into critical_packing_parameter: this paper's
    central finding is that a_e (and hence the packing parameter itself)
    is NOT independent of tail length for ionic surfactants, contrary to
    the common assumption that only the headgroup sets it -- a_e
    decreases with increasing chain length via the Debye screening
    length, which itself depends on the cmc. Use the returned a_e as the
    head_area_A2 input to critical_packing_parameter for an ionic
    surfactant when comparing across a homologous series matters.
    cmc_M: no added electrolyte (uses the paper's own approximation that
    ionic strength ~= cmc near the CMC for a 1:1 ionic surfactant).
    headgroup_prefactor_A: sqrt(2*pi*e^2*d/(epsilon*sigma)) in Angstrom
    -- a real, headgroup/system-specific physical quantity (d = charge-
    separation distance, sigma = bare interfacial free energy per unit
    area); the source paper's own illustrative value for sodium-alkyl-
    sulfate-like headgroups is 82.0 -- do NOT treat 82 as universal for
    other headgroups without a source."""
    a_e = cpp_mod.nagarajan_equilibrium_area_ionic(cmc_M, tail_length_A, headgroup_prefactor_A, temperature_K, dielectric_constant)
    kappa_inv = cpp_mod.nagarajan_debye_huckel_kappa_inverse(cmc_M, temperature_K, dielectric_constant)
    return {"a_e_A2": a_e, "debye_huckel_kappa_inverse_A": kappa_inv, "unit": "square Angstrom (a_e), Angstrom (kappa^-1)"}


@mcp.tool()
def aggregation_number(tail_volume_A3: float, core_radius_A: float) -> dict:
    """Estimated spherical-micelle aggregation number from geometric
    packing: N_agg = core_volume / tail_volume. tail_volume_A3 typically
    comes from tanford_chain_geometry; core_radius_A is usually
    approximated by the critical chain length from that same tool, but
    pass a measured value (e.g. from SANS/SAXS) if you have one. This is
    a geometric estimate (instant, no real measurement), NOT a substitute
    for a directly measured aggregation number -- if you have raw
    fluorescence-quenching or static-light-scattering data, use
    aggregation_number_from_quenching or aggregation_number_from_sls
    instead; those give a real, measured value from real instrumentation
    a smaller lab is actually likely to have (a standard
    spectrofluorometer, or a single-angle light-scattering/DLS
    instrument), not a packing-geometry approximation."""
    value = cpp_mod.aggregation_number_spherical(tail_volume_A3, core_radius_A)
    return {"aggregation_number": value, "unit": "dimensionless (molecules per micelle)",
            "note": "geometric estimate assuming a spherical micelle core"}


@mcp.tool()
def aggregation_number_from_quenching(
    quencher_concentrations: list[float],
    intensities: list[float],
    intensity_without_quencher: float,
    total_surfactant_concentration: float,
    cmc: float,
) -> dict:
    """Micelle aggregation number from a RAW steady-state fluorescence
    quenching (SSFQ) data series -- the Turro-Yekta method (1978):
    ln(I0/I) = N*[Q] / (Ct - cmc), a straight line through the origin.
    The SIMPLEST real (measured, not geometric) aggregation-number
    method available: needs only a standard steady-state
    spectrofluorometer (not the time-resolved/TCSPC instrumentation
    TRFQ needs, not SANS, not an ultracentrifuge) -- a fluorescent probe
    (commonly pyrene) solubilized in the micelles, quenched by a
    compound that partitions almost entirely into the micellar phase
    (e.g. cetylpyridinium chloride, benzophenone), measured at a series
    of quencher concentrations. quencher_concentrations/intensities are
    paired raw data; intensity_without_quencher is I0 (measured
    separately, [Q]=0); total_surfactant_concentration and cmc must use
    the SAME concentration unit as quencher_concentrations, and
    total_surfactant_concentration must exceed cmc."""
    result = curve.aggregation_number_from_quenching_curve(
        quencher_concentrations, intensities, intensity_without_quencher, total_surfactant_concentration, cmc
    )
    return {
        "aggregation_number": result.aggregation_number,
        "unit": "dimensionless (molecules per micelle)",
        "method": result.method,
    }


@mcp.tool()
def aggregation_number_from_sls(
    micellized_concentrations_g_per_mL: list[float],
    rayleigh_ratios_cm_inv: list[float],
    dn_dc_mL_per_g: float,
    wavelength_nm: float,
    monomer_molar_mass_g_per_mol: float,
    refractive_index_solvent: float = 1.333,
) -> dict:
    """Micelle aggregation number from a RAW static light scattering
    (SLS) concentration series -- the single-angle Debye plot method,
    real and standard for small (Rg < 12 nm) scatterers like surfactant
    micelles (no multi-angle instrument needed, unlike larger
    macromolecules): Kc/DeltaR = 1/Mw + 2*A2*c. A real measured value,
    achievable on the same class of light-scattering/DLS instrument
    already common for micelle sizing -- more instrumentation than
    aggregation_number_from_quenching but still far more accessible than
    SANS or an ultracentrifuge. micellized_concentrations_g_per_mL MUST
    be the concentration of MICELLIZED surfactant only (total
    concentration minus the cmc's mass-equivalent), not total surfactant
    concentration -- free monomer contributes no excess scattering and
    would bias the fit. rayleigh_ratios_cm_inv are real, solvent-
    subtracted excess Rayleigh ratios (not raw photon counts).
    dn_dc_mL_per_g (refractive index increment) is system-specific and
    must be a real measured/cited value -- never guess it. wavelength_nm
    is the VACUUM laser wavelength."""
    result = curve.aggregation_number_from_sls_debye_plot(
        micellized_concentrations_g_per_mL, rayleigh_ratios_cm_inv, dn_dc_mL_per_g,
        wavelength_nm, monomer_molar_mass_g_per_mol, refractive_index_solvent,
    )
    return {
        "aggregation_number": result.aggregation_number,
        "micelle_molar_mass_g_per_mol": result.micelle_molar_mass_g_per_mol,
        "second_virial_coefficient_cm3_mol_per_g2": result.second_virial_coefficient_cm3_mol_per_g2,
        "unit": "dimensionless (molecules per micelle)",
        "method": result.method,
    }


@mcp.tool()
def aggregation_number_from_dls(diffusion_coefficient_cm2_per_s: float, viscosity_mPas: float, partial_specific_volume_cm3_per_g: float, monomer_molar_mass_g_per_mol: float, temperature_K: float = 298.15) -> dict:
    """Micelle aggregation number from a RAW DLS-measured diffusion
    coefficient -- a real, sourced ALTERNATIVE to
    aggregation_number_from_sls/aggregation_number_from_quenching, using
    the same instrument class (DLS) already common for micelle sizing.
    Converts the Stokes-Einstein hydrodynamic radius to a hydrodynamic
    volume, then to a molar mass via the micelle's partial specific
    volume, then to an aggregation number via the monomer molar mass.
    diffusion_coefficient_cm2_per_s and viscosity_mPas: same units/
    convention as hydrodynamic_radius. partial_specific_volume_cm3_per_g:
    real, measured/cited value for this surfactant -- system-specific,
    do not guess. Real, disclosed limitation: the hydrodynamic radius
    includes a thin hydration shell, so this gives an aggregation number
    that is a real, usable estimate but typically biased slightly HIGH
    relative to a hydration-corrected value -- same estimate character
    as the geometric aggregation_number tool, not a hydration-corrected
    exact value."""
    result = curve.aggregation_number_from_dls(diffusion_coefficient_cm2_per_s, viscosity_mPas, partial_specific_volume_cm3_per_g, monomer_molar_mass_g_per_mol, temperature_K)
    return {
        "aggregation_number": result.aggregation_number,
        "hydrodynamic_radius_nm": result.hydrodynamic_radius_nm,
        "micelle_molar_mass_g_per_mol": result.micelle_molar_mass_g_per_mol,
        "unit": "dimensionless (molecules per micelle)",
        "method": result.method,
    }


@mcp.tool()
def aggregation_number_from_svedberg(sedimentation_coefficient_S: float, diffusion_coefficient_cm2_per_s: float, partial_specific_volume_cm3_per_g: float, solvent_density_g_per_cm3: float, monomer_molar_mass_g_per_mol: float, temperature_K: float = 298.15) -> dict:
    """Micelle aggregation number from RAW analytical-ultracentrifugation
    data via the classic Svedberg equation -- a real, sourced
    ALTERNATIVE independent of the light-scattering/DLS routes (a
    different instrument, different physical principle: M =
    R*T*s/[D*(1-v_bar*rho)]). sedimentation_coefficient_S must be in
    SVEDBERGS (1 S = 1e-13 s), the standard reporting unit -- do not
    pass raw seconds. diffusion_coefficient_cm2_per_s: same convention
    as elsewhere in this server. partial_specific_volume_cm3_per_g and
    solvent_density_g_per_cm3: real, measured/cited values -- system-
    specific, do not guess (water ~0.997 g/cm^3 at 25C for the solvent
    density is a real physical constant, not a guess). Raises if the
    buoyancy factor (1-v_bar*rho) is non-positive -- the particle would
    not physically sediment under those inputs."""
    result = curve.aggregation_number_from_svedberg_equation(sedimentation_coefficient_S, diffusion_coefficient_cm2_per_s, partial_specific_volume_cm3_per_g, solvent_density_g_per_cm3, monomer_molar_mass_g_per_mol, temperature_K)
    return {
        "aggregation_number": result.aggregation_number,
        "micelle_molar_mass_g_per_mol": result.micelle_molar_mass_g_per_mol,
        "unit": "dimensionless (molecules per micelle)",
        "method": result.method,
    }


@mcp.tool()
def cmc_from_surface_tension_curve(concentrations_mM: list[float], surface_tensions_mN_per_m: list[float]) -> dict:
    """Extract CMC from a RAW tensiometry curve (surface tension vs.
    concentration) via two-segment linear-regression break-point
    detection -- the missing first step other tools assume already
    happened (they all take an already-known CMC as input). Real,
    standard method (Shah/Das/Bhattarai 2025, Heliyon, PMC11835642), not
    invented: fits every possible split into a declining premicellar
    line and a flat post-CMC segment, picks the split minimizing total
    residual sum of squares. Also returns the premicellar slope in the
    exact form gibbs_surface_excess expects, so this chains directly
    into that tool without unit conversion."""
    result = curve.cmc_from_surface_tension_curve(concentrations_mM, surface_tensions_mN_per_m)
    return {
        "cmc_mM": result.cmc_mM,
        "premicellar_slope_mN_per_m_per_lnC": result.premicellar_slope_mN_per_m_per_lnC,
        "gamma_at_cmc_mN_per_m": result.gamma_at_cmc_mN_per_m,
        "r_squared_premicellar": result.r_squared_premicellar,
        "n_premicellar_points": result.n_premicellar_points,
        "n_postmicellar_points": result.n_postmicellar_points,
        "method": result.method,
    }


@mcp.tool()
def cmc_from_conductivity(concentrations_mM: list[float], conductivities: list[float]) -> dict:
    """Extract CMC from a RAW specific-conductivity-vs-concentration curve --
    the standard conductometric break-point method for IONIC surfactants
    only (nonionic surfactants have negligible conductivity response and
    cannot be characterized this way). Unlike cmc_from_surface_tension_curve
    (fit on log10 concentration), this fits on LINEAR concentration: two
    real straight-line segments (premicellar and postmicellar), break point
    at their intersection = CMC. Returns both slopes AND the counterion
    binding degree (beta) computed from them the same way
    thermodynamics.counterion_binding_degree() does, so this one call
    answers "raw conductivity curve -> CMC and beta" end to end."""
    result = curve.cmc_from_conductivity_curve(concentrations_mM, conductivities)
    return {
        "cmc_mM": result.cmc_mM,
        "slope_below_cmc": result.slope_below_cmc,
        "slope_above_cmc": result.slope_above_cmc,
        "counterion_binding_degree": result.counterion_binding_degree,
        "r_squared_below_cmc": result.r_squared_below_cmc,
        "r_squared_above_cmc": result.r_squared_above_cmc,
        "n_points_below_cmc": result.n_points_below_cmc,
        "n_points_above_cmc": result.n_points_above_cmc,
        "method": result.method,
    }


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
def grahame_equation_surface_potential(surface_charge_density_C_per_m2: float, ionic_strength_M: float, z: int = 1, temperature_K: float = 298.15) -> dict:
    """Surface (Stern-layer) potential (mV) from surface charge density
    via the Grahame equation -- the standard, textbook closed-form
    Gouy-Chapman diffuse-double-layer relation (Grahame, Chem. Rev. 41
    (1947) 441-501): sigma = sqrt(8*eps_r*eps0*k_B*T*n0)*sinh(z*e*psi/(2*k_B*T)),
    solved for psi. z: charge number of the symmetric supporting
    electrolyte. ionic_strength_M: same convention as
    debye_screening_length (ionic strength, not raw salt concentration).
    Real use: convert a surfactant's surface excess Gamma (from
    gibbs_surface_excess) to a charge density via sigma=z*F*Gamma
    (F=Faraday constant) to estimate the electrostatic potential of an
    adsorbed ionic monolayer -- a real Gouy-Chapman electrostatics
    calculation, more rigorous than assuming a fixed Gibbs prefactor."""
    value = elec.grahame_equation_surface_potential(surface_charge_density_C_per_m2, ionic_strength_M, z, temperature_K)
    return {"surface_potential_mV": value, "unit": "mV"}


@mcp.tool()
def zeta_potential(electrophoretic_mobility_um_cm_per_Vs: float, viscosity_mPas: float, regime: str, kappa_a: float | None = None) -> dict:
    """Zeta potential (mV) from electrophoretic mobility via the Henry
    equation. regime must be explicitly stated as 'huckel' (small
    particles / low ionic strength, kappa*a << 1), 'smoluchowski'
    (larger colloids, ionic strength >= ~10 mM, kappa*a >> 1),
    'ohshima', 'swan_furst', or 'qin' (three closed-form approximations
    valid and accurate across the FULL kappa*a range, including the
    intermediate zone where neither limit is a good approximation --
    typical real micelles, R~2-4 nm at 1-100 mM ionic strength, land at
    kappa*a ~ 0.2-4.2, squarely in that intermediate zone; all three
    require kappa_a. Accuracy vs. the numerically-exact Henry function:
    ohshima <3% error typically, qin <1.5% (most accurate specifically
    in the kappa*a ~ 1-20 range real micelles land in), swan_furst
    <0.1% overall). There is no default regime, because silently
    picking the wrong one is exactly the kind of error this tool exists
    to prevent. viscosity_mPas: solvent viscosity in mPa.s (=cP; water
    is ~0.89 at 25C -- do not pass Pa.s here, a 1000x unit error).
    kappa_a: the dimensionless product kappa*a (Debye parameter times
    particle radius, same units for both) -- required only for
    regime='ohshima', 'swan_furst', or 'qin'."""
    value = elec.zeta_potential_henry(electrophoretic_mobility_um_cm_per_Vs, viscosity_mPas, regime, kappa_a=kappa_a)
    return {"zeta_potential_mV": value, "unit": "mV", "regime_used": regime.lower()}


@mcp.tool()
def predict_mobility_relaxation_corrected(zeta_mV: float, kappa_a: float, viscosity_mPas: float, temperature_K: float = 298.15, z: int = 1, m: float = 0.15) -> dict:
    """Predict electrophoretic mobility ((um*cm)/(V*s)) from a known/
    assumed zeta potential via the relaxation-corrected model -- the
    forward direction of zeta_potential_relaxation_corrected (see that
    tool's docstring for full sourcing and parameter meaning). Useful
    for comparing a simulated or literature zeta potential against what
    mobility it would actually produce, accounting for the relaxation
    effect Henry's plain formula misses at high zeta."""
    value = elec.electrophoretic_mobility_relaxation_corrected(zeta_mV, kappa_a, viscosity_mPas, temperature_K, z, m)
    return {"electrophoretic_mobility_um_cm_per_Vs": value, "unit": "(um*cm)/(V*s)"}


@mcp.tool()
def zeta_potential_relaxation_corrected(electrophoretic_mobility_um_cm_per_Vs: float, kappa_a: float, viscosity_mPas: float, temperature_K: float = 298.15, z: int = 1, m: float = 0.15) -> dict:
    """Zeta potential (mV) from electrophoretic mobility -- a real,
    sourced ALTERNATIVE to zeta_potential (plain Henry equation) for
    HIGH zeta potential (>50 mV, common for ionic surfactant micelles):
    Henry's formula is LINEAR in zeta and misses the relaxation effect
    entirely (counterion relaxation around the moving particle makes
    mobility grow sub-linearly with zeta, and it can even pass through
    a maximum beyond ~100 mV). Uses O'Brien's simplification of the
    Dukhin-Semenikhin relaxation equation (J. Colloid Interface Sci. 309
    (2007) 194-224, Eq. 24), requires kappa_a > ~20 (the source's own
    stated valid range -- raises otherwise, use zeta_potential/'ohshima'
    etc. for lower kappa_a). z: counterion/co-ion charge number
    (symmetrical z-z electrolyte assumed). m: dimensionless ionic
    mobility parameter (source's Eq. 18) -- default 0.15 is the source's
    own stated typical value for aqueous solutions; pass the real value
    for the actual counterion if known. Solved numerically (grid+
    bisection) and restricted to zeta in [0, 150] mV by default (below
    the reported mobility maximum) -- raises if the given mobility has
    no solution in that range rather than silently returning a wrong
    branch."""
    value = elec.zeta_potential_relaxation_corrected(electrophoretic_mobility_um_cm_per_Vs, kappa_a, viscosity_mPas, temperature_K, z, m)
    return {"zeta_potential_mV": value, "unit": "mV", "model": "relaxation-corrected (O'Brien/Dukhin-Semenikhin)"}


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
def hydrodynamic_radius_perrin_corrected(diffusion_coefficient_cm2_per_s: float, viscosity_mPas: float, axial_ratio: float, temperature_K: float = 298.15) -> dict:
    """TRUE equivalent-sphere hydrodynamic radius (nm) for a non-
    spherical (rodlike/wormlike or disklike) micelle from a DLS-measured
    diffusion coefficient -- a real, sourced ALTERNATIVE to
    hydrodynamic_radius, which assumes a sphere and so returns an
    inflated apparent radius for elongated/flattened aggregates. Uses
    the Perrin friction factor (Perrin 1934) to correct for shape:
    R_true = R_h,apparent / F_P(axial_ratio). axial_ratio = a/b (semi-
    axis along the rotation axis / equatorial semi-axis) must be known
    or assumed independently (e.g. from the critical packing parameter's
    predicted morphology, or SANS/SAXS) -- axial_ratio > 1 is prolate
    (rodlike), < 1 is oblate (disklike), == 1 gives exactly the same
    result as hydrodynamic_radius (no correction needed for a sphere).
    Do not guess axial_ratio."""
    r_true = dyn.hydrodynamic_radius_perrin_corrected(diffusion_coefficient_cm2_per_s, viscosity_mPas, axial_ratio, temperature_K)
    f_p = dyn.perrin_friction_factor(axial_ratio)
    return {"hydrodynamic_radius_true_nm": r_true, "perrin_friction_factor": f_p, "unit": "nm"}


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
def mass_action_free_energy_of_micellization(cmc_M: float, aggregation_number: float, temperature_K: float) -> dict:
    """Standard Gibbs free energy of micellization (kJ/mol) via the
    MASS-ACTION-LAW model -- a real, sourced ALTERNATIVE to
    gibbs_free_energy_micellization's pseudo-phase-separation
    convention: deltaG_mic = R*T*[(1+1/n)*ln(X_cmc) - (1/n)*ln(n)].
    The pseudo-phase result is the n -> infinity LIMIT of this formula
    (they converge for large aggregation number, roughly n>~50-100);
    for smaller aggregation numbers the correction terms here are real
    and non-negligible. NONIONIC SURFACTANTS ONLY -- this is the simple
    single-equilibrium (n monomers <-> 1 micelle) mass-action result
    with no counterion-association step; do not apply to an ionic
    surfactant (use gibbs_free_energy_micellization with the appropriate
    counterion_factor for that case instead). aggregation_number must
    be a real, known value (e.g. from aggregation_number_from_quenching,
    aggregation_number_from_sls, or aggregation_number) -- do not
    guess it."""
    x_cmc = thermo.cmc_to_mole_fraction(cmc_M)
    value = thermo.mass_action_free_energy_of_micellization(x_cmc, aggregation_number, temperature_K)
    return {"deltaG_mic_kJ_per_mol": value, "unit": "kJ/mol", "cmc_mole_fraction_used": x_cmc, "model": "mass-action-law"}


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
def vant_hoff_multi_point_fit(cmc_series_M: list[float], temperatures_K: list[float], counterion_factor: float = 1.0, reference_temperature_K: float | None = None) -> dict:
    """Full multi-point (3+ temperatures, 5+ recommended) Gibbs-Helmholtz
    fit for deltaH, deltaS, AND deltaCp simultaneously from a raw
    CMC-vs-temperature series -- the real fix for vant_hoff_enthalpy's
    strictly-2-point limitation (which cannot separate a real deltaCp
    from noise), per Kantonen, Henriksen & Gilson 2018 (BBA Gen. Subj.
    1862, 692-704). cmc_series_M (molarity) and temperatures_K must be
    the same length and in matching order. counterion_factor: 1.0 for
    nonionic, (2-beta) for ionic (same convention as
    gibbs_free_energy_micellization) -- not auto-detected, do not guess.
    reference_temperature_K defaults to the mean of temperatures_K for a
    well-conditioned fit; deltaH/deltaS are reported AT this reference
    temperature (a real physical consequence of nonzero deltaCp) while
    deltaCp itself is reference-independent. Returns r_squared as a real
    diagnostic: a poor fit signals deltaCp is not actually constant over
    the measured range, or the data is noisy relative to the temperature
    span used."""
    x_values = [thermo.cmc_to_mole_fraction(c) for c in cmc_series_M]
    result = thermo.vant_hoff_multi_point_fit(x_values, temperatures_K, counterion_factor, reference_temperature_K)
    return {
        "deltaH_mic_kJ_per_mol": result.delta_H_kJ_per_mol,
        "deltaS_mic_J_per_mol_K": result.delta_S_J_per_mol_K,
        "deltaCp_mic_J_per_mol_K": result.delta_Cp_J_per_mol_K,
        "reference_temperature_K": result.reference_temperature_K,
        "r_squared": result.r_squared,
        "n_points": result.n_points,
        "deltaG_fitted_kJ_per_mol": result.delta_G_fitted_kJ_per_mol,
    }


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
def owens_wendt_solid_surface_energy(contact_angles_deg: list[float], liquid_gamma_dispersive_mN_m: list[float], liquid_gamma_polar_mN_m: list[float]) -> dict:
    """Decompose a SOLID's surface energy into dispersive and polar
    components via the OWRK method (Owens & Wendt 1969; Rabel 1971;
    Kaelble 1970) -- a real, sourced ALTERNATIVE to
    wetting_work_of_adhesion's single scalar number, needed only when
    decomposing surface energy into dispersive/polar components matters
    (not a more-correct answer to the same question). Combines Young-
    Dupre's work of adhesion with the OWRK mixing rule W = 2*[sqrt(gamma_S^d*gamma_L^d)
    + sqrt(gamma_S^p*gamma_L^p)], solved by ordinary least squares
    across the given test liquids. contact_angles_deg,
    liquid_gamma_dispersive_mN_m, and liquid_gamma_polar_mN_m must be
    the same length (>=2 test liquids, the standard minimum) and in
    matching order -- real, literature-tabulated per-liquid values
    (e.g. water: d~21.8, p~51.0 mN/m; diiodomethane: d~50.8, p~0 mN/m)
    must be supplied by the caller, not guessed. Returns gamma_S^d,
    gamma_S^p, gamma_S_total, and r_squared (a real fit diagnostic when
    3+ liquids are given)."""
    result = wetting.owens_wendt_solid_surface_energy(contact_angles_deg, liquid_gamma_dispersive_mN_m, liquid_gamma_polar_mN_m)
    return {
        "gamma_s_dispersive_mN_m": result.gamma_s_dispersive_mN_m,
        "gamma_s_polar_mN_m": result.gamma_s_polar_mN_m,
        "gamma_s_total_mN_m": result.gamma_s_total_mN_m,
        "r_squared": result.r_squared,
        "n_liquids": result.n_liquids,
        "unit": "mN/m (= mJ/m^2)",
    }


@mcp.tool()
def van_oss_chaudhury_good_solid_surface_energy(contact_angles_deg: list[float], liquid_gamma_lw_mN_m: list[float], liquid_gamma_acid_mN_m: list[float], liquid_gamma_base_mN_m: list[float]) -> dict:
    """Decompose a SOLID's surface energy into Lifshitz-van der Waals
    (LW) and Lewis acid/base components -- van Oss, Chaudhury & Good's
    three-liquid method (Chem. Rev. 88 (1988) 927-941), a real, sourced
    ALTERNATIVE to wetting_work_of_adhesion and
    owens_wendt_solid_surface_energy, needed when the acid-base
    (electron-donor/electron-acceptor) character of the surface matters,
    not just polarity in general. Combines Young-Dupre's work of
    adhesion with the vOCG mixing rule W = 2*[sqrt(gamma_S^LW*gamma_L^LW)
    + sqrt(gamma_S^acid*gamma_L^base) + sqrt(gamma_S^base*gamma_L^acid)]
    -- note the CROSS terms (solid acid pairs with liquid base, and vice
    versa). contact_angles_deg, liquid_gamma_lw_mN_m,
    liquid_gamma_acid_mN_m, and liquid_gamma_base_mN_m must all be
    length EXACTLY 3 (the standard minimum: one purely apolar liquid,
    e.g. diiodomethane, plus two polar liquids with independent acid/
    base character, e.g. water and glycerol) -- real, literature-
    tabulated per-liquid values must be supplied, not guessed. Returns
    gamma_S^LW, gamma_S^acid, gamma_S^base, and gamma_S_total."""
    result = wetting.van_oss_chaudhury_good_solid_surface_energy(contact_angles_deg, liquid_gamma_lw_mN_m, liquid_gamma_acid_mN_m, liquid_gamma_base_mN_m)
    return {
        "gamma_s_lw_mN_m": result.gamma_s_lw_mN_m,
        "gamma_s_acid_mN_m": result.gamma_s_acid_mN_m,
        "gamma_s_base_mN_m": result.gamma_s_base_mN_m,
        "gamma_s_total_mN_m": result.gamma_s_total_mN_m,
        "unit": "mN/m (= mJ/m^2)",
    }


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
    same concentration unit. See micelle_water_partition_coefficient for
    a real complementary metric (intensive, not extensive)."""
    value = solub.molar_solubilization_ratio(total_solubilized_M, intrinsic_water_solubility_M, surfactant_concentration_M, cmc_M)
    return {"msr": value, "unit": "dimensionless (mol solubilizate per mol micellized surfactant)"}


@mcp.tool()
def micelle_water_partition_coefficient(total_solubilized_M: float, intrinsic_water_solubility_M: float, surfactant_concentration_M: float, cmc_M: float) -> dict:
    """Micelle-water partition coefficient (Km, dimensionless),
    MOLE-FRACTION-BASIS convention: Km = X_solubilizate,micelle /
    X_solubilizate,water -- a real, sourced ALTERNATIVE/complementary
    metric to molar_solubilization_ratio, capturing how strongly the
    solubilizate favors the micellar pseudophase over bulk water
    (intensive), whereas MSR is an extensive uptake-per-mole quantity. A
    high Km with a low MSR (or vice versa) is real information MSR alone
    misses. IMPORTANT, disclosed honestly: this is ONE of (at least) two
    real conventions in the literature -- some papers instead use a
    molar-CONCENTRATION-RATIO basis (no mole-fraction normalization),
    which gives a numerically DIFFERENT value for the same raw data.
    Verify which convention a paper you're comparing against actually
    uses before treating numbers as directly comparable. Same four
    inputs as molar_solubilization_ratio, same units, same meaning."""
    value = solub.micelle_water_partition_coefficient(total_solubilized_M, intrinsic_water_solubility_M, surfactant_concentration_M, cmc_M)
    return {"km": value, "unit": "dimensionless (mole-fraction-basis convention)"}


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
