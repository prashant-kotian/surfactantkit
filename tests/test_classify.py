"""Tests for surfactantkit.classify.classify_surfactant_charge_type -- the
SMILES-based charge-type classifier that closes a real gap: every function
elsewhere in this library that needs a surfactant's ionic character
(gibbs_gamma_max, szyszkowski_fit_K, ...) required it as an explicit
human-supplied string, with no way to determine it from the surfactant's
own structure. Validated against real, well-known surfactant SMILES
spanning every charge class, including the deliberately ambiguous
(free-amine, no stated pH) case this function must NOT silently resolve.
"""

import pytest

from surfactantkit.classify import classify_surfactant_charge_type, classify_surfactant_structural_family


def test_sds_anionic_sulfate():
    r = classify_surfactant_charge_type("CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]")
    assert r.charge_type == "anionic"
    assert r.system_type_for_gibbs == "ionic_no_added_salt"
    assert r.confidence == "high"


def test_sds_acid_form_still_anionic():
    r = classify_surfactant_charge_type("CCCCCCCCCCCCOS(=O)(=O)O")
    assert r.charge_type == "anionic"


def test_sdbs_anionic_sulfonate():
    r = classify_surfactant_charge_type("CCCCCCCCCCCCc1ccc(cc1)S(=O)(=O)[O-].[Na+]")
    assert r.charge_type == "anionic"


def test_sodium_oleate_anionic_carboxylate_moderate_confidence():
    r = classify_surfactant_charge_type("CCCCCCCC/C=C\\CCCCCCCC(=O)[O-].[Na+]")
    assert r.charge_type == "anionic"
    assert r.confidence == "moderate"  # carboxylate is pH-sensitive near its own pKa, disclosed


def test_ctab_cationic_quaternary_ammonium():
    r = classify_surfactant_charge_type("CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]")
    assert r.charge_type == "cationic"
    assert r.system_type_for_gibbs == "ionic_no_added_salt"
    assert r.confidence == "high"


def test_dtab_cationic():
    r = classify_surfactant_charge_type("CCCCCCCCCCCC[N+](C)(C)C.[Br-]")
    assert r.charge_type == "cationic"


def test_capb_zwitterionic():
    r = classify_surfactant_charge_type("CCCCCCCCCCCCC(=O)NCCC[N+](C)(C)CC(=O)[O-]")
    assert r.charge_type == "zwitterionic"
    assert r.system_type_for_gibbs is None  # deliberately not defaulted -- real open question


def test_c12e8_nonionic():
    r = classify_surfactant_charge_type("CCCCCCCCCCCCOCCOCCOCCOCCOCCOCCOCCO")
    assert r.charge_type == "nonionic"
    assert r.system_type_for_gibbs == "nonionic"
    assert r.confidence == "moderate"  # absence-based claim, disclosed as such


def test_free_amine_is_ambiguous_not_guessed():
    """A free (non-quaternary) amine's charge state genuinely depends on
    solution pH, which this function is never given -- must NOT silently
    default to cationic."""
    r = classify_surfactant_charge_type("CCCCCCCCCCCCN")
    assert r.charge_type == "ambiguous_pH_dependent"
    assert r.system_type_for_gibbs is None


def test_unparseable_smiles_reported_not_raised():
    r = classify_surfactant_charge_type("not a smiles $$$")
    assert r.charge_type == "unparseable"
    assert r.system_type_for_gibbs is None
    assert r.confidence == "none"


# --- classify_surfactant_structural_family -------------------------------
# Real, sourced structures, not invented: gemini G6/12-4-12 are the same
# literature-verified compounds already used in SurfQSPR's own
# GEMINI_CANDIDATE_STRUCTURES (Azum et al. 2022). Monorhamnolipid and
# sophorolipid SMILES pulled directly from PubChem (CID 162246, CID
# 11856871) via live lookup before writing this test, not recalled from
# training -- see classify.py's own module comment for full provenance.

GEMINI_G6 = "[Br-].[Br-].CCCCCCCCCCCCCCCC[N+](C)(C)CCCCCC[N+](C)(C)CCCCCCCCCCCCCCCC"
GEMINI_12_4_12 = "[Br-].[Br-].CCCCCCCCCCCC[N+](C)(C)CCCC[N+](C)(C)CCCCCCCCCCCC"
MONORHAMNOLIPID = "CCCCCCCC(CC(=O)O)OC(=O)CC(CCCCCCC)O[C@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)C)O)O)O"
SOPHOROLIPID = "CC(CCCCCC/C=C/CCCCCCCC(=O)O)O[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)COC(=O)C)O)O)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)COC(=O)C)O)O)O"
GLUCOSE_ALONE = "C(C1C(C(C(C(O1)O)O)O)O)O"


def test_gemini_g6_detected_as_dimeric():
    r = classify_surfactant_structural_family(GEMINI_G6)
    assert r.structural_family == "dimeric_gemini_type"
    assert r.n_strong_ionic_headgroups == 2
    assert r.confidence == "high"
    assert r.caveats  # the bolaform-ambiguity disclosure must be present


def test_gemini_12_4_12_detected_as_dimeric():
    r = classify_surfactant_structural_family(GEMINI_12_4_12)
    assert r.structural_family == "dimeric_gemini_type"
    assert r.n_strong_ionic_headgroups == 2


def test_dtab_monomeric_not_dimeric():
    r = classify_surfactant_structural_family("CCCCCCCCCCCC[N+](C)(C)C.[Br-]")
    assert r.structural_family == "monomeric"
    assert r.n_strong_ionic_headgroups == 1


def test_ctab_monomeric_not_dimeric():
    r = classify_surfactant_structural_family("CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]")
    assert r.structural_family == "monomeric"
    assert r.n_strong_ionic_headgroups == 1


def test_sds_monomeric_no_sugar_ring():
    r = classify_surfactant_structural_family("CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]")
    assert r.structural_family == "monomeric"
    assert r.n_sugar_rings == 0


def test_monorhamnolipid_detected_as_glycolipid_biosurfactant():
    r = classify_surfactant_structural_family(MONORHAMNOLIPID)
    assert r.structural_family == "glycolipid_biosurfactant"
    assert r.n_sugar_rings == 1
    assert r.n_long_chain_matches >= 1
    assert r.confidence == "high"


def test_sophorolipid_detected_as_glycolipid_biosurfactant():
    r = classify_surfactant_structural_family(SOPHOROLIPID)
    assert r.structural_family == "glycolipid_biosurfactant"
    assert r.n_sugar_rings == 2  # sophorose is a disaccharide -- two rings, unlike rhamnolipid's one


def test_bare_glucose_not_a_biosurfactant():
    """A sugar with no lipid tail is not a surfactant at all -- the long-
    chain requirement must correctly exclude it, not just check for a
    sugar ring alone."""
    r = classify_surfactant_structural_family(GLUCOSE_ALONE)
    assert r.structural_family == "monomeric"
    assert r.n_sugar_rings == 1
    assert r.n_long_chain_matches == 0


def test_structural_family_unparseable_smiles_reported_not_raised():
    r = classify_surfactant_structural_family("not a smiles $$$")
    assert r.structural_family == "unparseable"
    assert r.confidence == "none"
