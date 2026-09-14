"""Tests for surfactantkit.orchestrate.derive_all_properties_from_smiles_and_curve
-- the actual deliverable this library was missing: given ONLY a real
surfactant SMILES and a raw surface-tension-vs-concentration curve (no
stated ionic character, no named isotherm model, no pre-fit parameters),
autonomously derive CMC, Gamma_max, A_min, the Langmuir-vs-Frumkin choice,
and deltaG_mic where the inputs actually support it.

The SDS/DTAB/CTAB cases use the real 8-point robotic-tensiometer datasets
from benchmark/paper3_pilot/build_pilot_50_v2.py (Dankloff et al. 2025,
PendantProp), with gold CMC/Gamma_max/A_min values already independently
verified in benchmark/paper3_pilot/pilot50_v2_questions.json by actually
calling cmc_from_surface_tension_curve + gibbs_gamma_max + gibbs_a_min --
this is a real regression guard for the min_points_per_segment fix in
curve_analysis.py (found via this exact orchestrator, see that module's
own docstring), not just a smoke test.
"""

import pytest

from surfactantkit.orchestrate import derive_all_properties_from_smiles_and_curve
from surfactantkit.adsorption import szyszkowski_surface_tension, frumkin_surface_tension

SDS_SMILES = "CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]"
DTAB_SMILES = "CCCCCCCCCCCC[N+](C)(C)C.[Br-]"
CTAB_SMILES = "CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]"


def test_sds_matches_independently_verified_gold_values():
    concs = [0.40577109375, 0.8115421875, 1.623084375, 3.24616875, 6.4923375, 12.984675, 25.96935, 51.9387]
    gammas = [69.85933451975876, 67.9370105943862, 62.892649575815454, 54.57301518905874,
              42.87619733622359, 37.91402789831731, 37.44295519808575, 36.79102755846683]
    r = derive_all_properties_from_smiles_and_curve(SDS_SMILES, concs, gammas, temperature_K=296.15)

    assert r.classification.charge_type == "anionic"
    assert r.system_type_used == "ionic_no_added_salt"
    assert r.cmc.cmc_mM == pytest.approx(9.181551744003432, rel=1e-6)
    assert r.gamma_max_mol_per_m2 == pytest.approx(3.4266231784578047e-06, rel=1e-6)


def test_dtab_matches_independently_verified_gold_values():
    concs = [0.25578125, 0.5115625, 1.023125, 2.04625, 4.0925, 16.37, 32.74]
    gammas = [70.26907673, 68.20699668, 65.39702756, 57.50386199, 49.23642286, 32.34750118, 36.22382546]
    r = derive_all_properties_from_smiles_and_curve(DTAB_SMILES, concs, gammas, temperature_K=296.35)

    assert r.classification.charge_type == "cationic"
    assert r.cmc.cmc_mM == pytest.approx(8.185, rel=1e-3)
    assert r.gamma_max_mol_per_m2 == pytest.approx(2.3655549359065027e-06, rel=1e-6)


def test_ctab_matches_independently_verified_gold_values():
    concs = [0.016484375, 0.03296875, 0.0659375, 0.131875, 0.26375, 0.5275, 1.055, 2.11]
    gammas = [71.60374462, 71.55862868, 70.38023925, 64.77683517, 57.3369305, 46.13427154, 36.41066498, 32.37490721]
    r = derive_all_properties_from_smiles_and_curve(CTAB_SMILES, concs, gammas, temperature_K=296.85)

    assert r.classification.charge_type == "cationic"
    assert r.cmc.cmc_mM == pytest.approx(0.7459976541518075, rel=1e-6)
    assert r.gamma_max_mol_per_m2 == pytest.approx(2.72426244417475e-06, rel=1e-6)


def test_ionic_surfactant_without_counterion_data_skips_delta_g_mic_with_named_gap():
    concs = [0.40577109375, 0.8115421875, 1.623084375, 3.24616875, 6.4923375, 12.984675, 25.96935, 51.9387]
    gammas = [69.85933451975876, 67.9370105943862, 62.892649575815454, 54.57301518905874,
              42.87619733622359, 37.91402789831731, 37.44295519808575, 36.79102755846683]
    r = derive_all_properties_from_smiles_and_curve(SDS_SMILES, concs, gammas, temperature_K=296.15)

    assert r.delta_g_mic_kJ_per_mol is None
    assert any("counterion" in g for g in r.gaps)


def test_ionic_surfactant_with_supplied_alpha_computes_delta_g_mic():
    concs = [0.40577109375, 0.8115421875, 1.623084375, 3.24616875, 6.4923375, 12.984675, 25.96935, 51.9387]
    gammas = [69.85933451975876, 67.9370105943862, 62.892649575815454, 54.57301518905874,
              42.87619733622359, 37.91402789831731, 37.44295519808575, 36.79102755846683]
    r = derive_all_properties_from_smiles_and_curve(
        SDS_SMILES, concs, gammas, temperature_K=296.15, counterion_dissociation_alpha=0.272,
    )
    assert r.delta_g_mic_kJ_per_mol is not None
    assert r.counterion_factor_used == pytest.approx(1.272)
    assert r.delta_g_mic_kJ_per_mol < 0  # micellization is spontaneous


def test_nonionic_surfactant_computes_delta_g_mic_with_default_counterion_factor():
    smiles = "CCCCCCCCCCCCOCCOCCOCCOCCOCCOCCOCCO"  # C12E8-like nonionic
    gamma0, gmax, K, T = 72.0, 3.0e-6, 8.0, 298.15
    concs = [0.001, 0.003, 0.006, 0.01, 0.02, 0.04, 0.08, 0.15, 0.3, 0.6, 1.2, 2.4]
    gammas = [szyszkowski_surface_tension(c, gamma0, gmax, K, "nonionic", T) for c in concs[:9]] + [31.0, 31.0, 31.0]

    r = derive_all_properties_from_smiles_and_curve(smiles, concs, gammas, temperature_K=T)
    assert r.classification.charge_type == "nonionic"
    assert r.system_type_used == "nonionic"
    assert r.counterion_factor_used == pytest.approx(1.0)
    assert r.delta_g_mic_kJ_per_mol is not None
    assert r.delta_g_mic_kJ_per_mol < 0


def test_zwitterionic_surfactant_refuses_to_guess_gibbs_prefactor():
    smiles = "CCCCCCCCCCCCC(=O)NCCC[N+](C)(C)CC(=O)[O-]"  # CAPB-like zwitterion
    concs = [0.40577109375, 0.8115421875, 1.623084375, 3.24616875, 6.4923375, 12.984675, 25.96935, 51.9387]
    gammas = [69.85933451975876, 67.9370105943862, 62.892649575815454, 54.57301518905874,
              42.87619733622359, 37.91402789831731, 37.44295519808575, 36.79102755846683]

    r = derive_all_properties_from_smiles_and_curve(smiles, concs, gammas, temperature_K=296.15)
    assert r.classification.charge_type == "zwitterionic"
    assert r.system_type_used is None
    assert r.gamma_max_mol_per_m2 is None
    assert r.a_min_nm2 is None
    assert r.isotherm is None
    assert r.delta_g_mic_kJ_per_mol is None
    # CMC itself is charge-independent, so it should still be computed
    assert r.cmc.cmc_mM > 0


def test_unparseable_smiles_still_returns_cmc_with_reported_gap():
    concs = [0.40577109375, 0.8115421875, 1.623084375, 3.24616875, 6.4923375, 12.984675, 25.96935, 51.9387]
    gammas = [69.85933451975876, 67.9370105943862, 62.892649575815454, 54.57301518905874,
              42.87619733622359, 37.91402789831731, 37.44295519808575, 36.79102755846683]

    r = derive_all_properties_from_smiles_and_curve("not a smiles $$$", concs, gammas, temperature_K=296.15)
    assert r.classification.charge_type == "unparseable"
    assert r.system_type_used is None
    assert r.cmc.cmc_mM > 0
    assert any("could not be parsed" in g for g in r.gaps)


def test_isotherm_selection_engages_with_enough_premicellar_points():
    smiles = "CCCCCCCCCCCCOCCOCCOCCOCCOCCOCCOCCO"
    gamma0, gmax, K, a, T = 72.0, 3.2e-6, 12.0, 0.0, 298.15
    concs = [0.005, 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.28, 0.5, 1.0, 2.0]
    gammas = [frumkin_surface_tension(c, gamma0, gmax, K, a, "nonionic", T) for c in concs[:9]] + [28.0, 28.0, 28.0]

    r = derive_all_properties_from_smiles_and_curve(smiles, concs, gammas, temperature_K=T)
    assert r.isotherm is not None
    assert r.isotherm.selected_model in ("langmuir", "frumkin")
