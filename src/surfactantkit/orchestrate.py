"""The actual point of this library: given ONLY a surfactant's real SMILES
and a raw pre-and-post-CMC surface-tension-vs-concentration curve -- the
same two things a real researcher has in hand before doing any analysis,
and nothing else (no stated ionic character, no named isotherm model, no
pre-fit parameters) -- autonomously derive every property this codebase
can legitimately derive, the way a competent researcher would work through
it by hand: classify the surfactant, find the CMC, decide whether the
pre-CMC adsorption is plain Langmuir or genuinely non-ideal (Frumkin), fit
whichever is actually justified, and compute the standard geometric and
thermodynamic quantities that follow.

Every other function in this library is deliberately atomic and refuses to
guess a single missing input (system_type, isotherm model, counterion
binding) rather than silently assuming one -- that discipline is a real
strength for a human calling one function with full context, but it means
NONE of those functions alone can answer "what can we learn from this
surfactant" without a human first making several judgment calls. This
module is that human's judgment, made explicit and inspectable: every
place a real ambiguity exists (zwitterionic/pH-dependent charge, unknown
electrolyte condition, no independent counterion-binding measurement) is
carried through to the final result as an explicit, named gap rather than
resolved by a hidden default -- the same "do not guess" discipline this
library already applies function-by-function, just applied once across
the whole pipeline instead of separately at every step.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .classify import classify_surfactant_charge_type, ChargeClassificationResult
from .curve_analysis import cmc_from_surface_tension_curve, CmcFromCurveResult
from .adsorption import (
    gibbs_gamma_max, gibbs_a_min, select_isotherm_model, IsothermModelSelectionResult,
)
from .thermodynamics import cmc_to_mole_fraction, gibbs_free_energy_micellization


@dataclass
class OrchestratorResult:
    smiles: str
    classification: ChargeClassificationResult
    system_type_used: str | None  # what was actually passed to Gibbs-prefactor-dependent steps, or None if undetermined
    cmc: CmcFromCurveResult
    gamma_max_mol_per_m2: float | None
    a_min_nm2: float | None
    isotherm: IsothermModelSelectionResult | None
    delta_g_mic_kJ_per_mol: float | None
    counterion_factor_used: float | None
    gaps: list[str] = field(default_factory=list)  # named, honest list of everything NOT determined and why


def derive_all_properties_from_smiles_and_curve(
    smiles: str,
    concentrations_mM: list[float],
    surface_tensions_mN_per_m: list[float],
    temperature_K: float = 298.15,
    electrolyte_condition: str | None = None,
    counterion_dissociation_alpha: float | None = None,
    isotherm_a_bounds: tuple[float, float] = (-4.0, 4.0),
) -> OrchestratorResult:
    """Full (SMILES, raw curve) -> all-derivable-properties pipeline.

    electrolyte_condition: optional override, one of 'excess_electrolyte'
    or 'no_added_salt' -- ONLY pass this if it is independently known
    about the real solution being characterized (e.g. stated in the
    experiment's own description). A bare SMILES carries no information
    about what else is dissolved in the water, so classify.py's own
    default ('ionic_no_added_salt', i.e. Gibbs prefactor n=2) is used
    for any ionic surfactant when this is left as None -- surfaced as an
    explicit gap in the result rather than silently assumed to be
    correct, since a real excess-electrolyte system would need n=1
    instead and get a genuinely different Gamma_max.

    counterion_dissociation_alpha: optional, ONLY pass this if
    independently measured (e.g. via the conductivity slope-ratio method,
    counterion_binding_degree, or a literature/EPR value like Bales et
    al.'s alpha=0.272 for SDS -- see gibbs_free_energy_micellization's
    own docstring). Without it, deltaG_mic for an IONIC surfactant is
    left undetermined (reported as a named gap) rather than guessed --
    this is real information a surface-tension curve alone cannot supply,
    the same reasoning gibbs_free_energy_micellization's own docstring
    already applies.

    isotherm_a_bounds: passed through to select_isotherm_model.
    """
    gaps: list[str] = []

    classification = classify_surfactant_charge_type(smiles)

    if classification.charge_type == "unparseable":
        system_type: str | None = None
        gaps.append(f"SMILES could not be parsed ('{smiles}') -- charge type, Gamma_max, A_min, "
                     "isotherm model, and deltaG_mic are all undetermined until a valid SMILES is supplied.")
    elif classification.charge_type in ("zwitterionic", "ambiguous_pH_dependent"):
        system_type = None
        gaps.append(f"charge_type='{classification.charge_type}' -- {classification.caveats[0] if classification.caveats else 'ambiguous from structure alone'} "
                     "Gamma_max, A_min, and the isotherm fit all require a known Gibbs prefactor, which cannot "
                     "be assigned here without resolving this ambiguity (e.g. independently confirming solution pH, "
                     "or the specific zwitterionic-system literature precedent for its prefactor convention).")
    else:
        system_type = classification.system_type_for_gibbs  # 'nonionic' or 'ionic_no_added_salt'
        if electrolyte_condition == "excess_electrolyte" and system_type == "ionic_no_added_salt":
            system_type = "ionic_excess_electrolyte"
        elif electrolyte_condition == "no_added_salt":
            pass  # already the default for ionic; no-op for nonionic
        elif electrolyte_condition is None and system_type == "ionic_no_added_salt":
            gaps.append("electrolyte_condition not supplied for an ionic surfactant -- defaulted to "
                        "'ionic_no_added_salt' (Gibbs prefactor n=2), the correct choice ONLY if the real "
                        "solution truly has no significant excess inert electrolyte; if it does, re-run with "
                        "electrolyte_condition='excess_electrolyte' (n=1) for a materially different Gamma_max.")
        if classification.caveats:
            gaps.extend(f"classification note ({classification.charge_type}): {c}" for c in classification.caveats)

    cmc = cmc_from_surface_tension_curve(concentrations_mM, surface_tensions_mN_per_m)

    gamma_max = a_min = None
    isotherm: IsothermModelSelectionResult | None = None
    if system_type is not None:
        gamma_max = gibbs_gamma_max(cmc.premicellar_slope_mN_per_m_per_lnC, system_type, temperature_K)
        a_min = gibbs_a_min(gamma_max)

        premicellar_pairs = sorted(
            (c, g) for c, g in zip(concentrations_mM, surface_tensions_mN_per_m)
            if cmc.premicellar_x_min_mM <= c < cmc.cmc_mM
        )
        if len(premicellar_pairs) >= 4:
            pre_c = [p[0] for p in premicellar_pairs]
            pre_g = [p[1] for p in premicellar_pairs]
            isotherm = select_isotherm_model(pre_c, pre_g, cmc.gamma0_mN_m, gamma_max, system_type, temperature_K, a_bounds=isotherm_a_bounds)
            gaps.extend(f"isotherm fit note: {w}" for w in isotherm.warnings)
        else:
            gaps.append(f"only {len(premicellar_pairs)} pre-CMC data point(s) available (need >=4) -- "
                        "Langmuir-vs-Frumkin isotherm model selection skipped; Gamma_max/A_min above come "
                        "from the curve's own premicellar slope directly and do not depend on this.")
    else:
        gaps.append("Gamma_max, A_min, and isotherm model selection all skipped: system_type undetermined (see above).")

    delta_g_mic = None
    counterion_factor_used = None
    if system_type == "nonionic":
        counterion_factor_used = 1.0
        delta_g_mic = gibbs_free_energy_micellization(cmc_to_mole_fraction(cmc.cmc_mM / 1000.0), temperature_K, counterion_factor_used)
    elif system_type in ("ionic_no_added_salt", "ionic_excess_electrolyte"):
        if counterion_dissociation_alpha is not None:
            counterion_factor_used = 1.0 + counterion_dissociation_alpha
            delta_g_mic = gibbs_free_energy_micellization(cmc_to_mole_fraction(cmc.cmc_mM / 1000.0), temperature_K, counterion_factor_used)
        else:
            gaps.append("deltaG_mic skipped for this ionic surfactant: no independent counterion-binding "
                        "measurement supplied (counterion_dissociation_alpha). This is real information a "
                        "surface-tension curve and SMILES alone cannot provide -- it requires an independent "
                        "measurement (e.g. conductivity slope-ratio method, see counterion_binding_degree) or "
                        "a literature value for this specific surfactant. Do not assume alpha; report the gap.")
    elif system_type is not None:
        gaps.append(f"deltaG_mic skipped: unrecognized system_type '{system_type}'.")

    return OrchestratorResult(
        smiles=smiles, classification=classification, system_type_used=system_type, cmc=cmc,
        gamma_max_mol_per_m2=gamma_max, a_min_nm2=a_min, isotherm=isotherm,
        delta_g_mic_kJ_per_mol=delta_g_mic, counterion_factor_used=counterion_factor_used, gaps=gaps,
    )
