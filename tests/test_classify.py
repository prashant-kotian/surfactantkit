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

from surfactantkit.classify import classify_surfactant_charge_type


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
