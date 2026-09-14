"""Structure-based surfactant charge-type classification from SMILES.

This closes a real, previously-missing capability: every function in this
project that needs to know whether a surfactant is ionic or nonionic
(gibbs_gamma_max, szyszkowski_fit_K, hld_ionic/hld_nonionic, ...) has always
required that as an explicit, human-supplied string parameter -- there was
no way to determine it from the surfactant's own structure. A researcher
handing over a real dataset AND the surfactant's real SMILES should not
also have to separately tell the toolkit "this is ionic" -- the toolkit
should be able to work that out itself, the same way a real chemist reading
a structure would.

Method: RDKit SMARTS substructure matching for the functional groups that
actually determine surfactant charge behavior in water -- sulfate esters,
sulfonates, carboxylates/carboxylic acids, phosphate esters (anionic,
strong-to-moderate acids, treated as ionized at normal surfactant testing
pH regardless of how the SMILES happens to write the protonation state,
since that's how these groups actually behave in real aqueous
formulations), quaternary ammonium (permanently cationic, pH-independent),
and free primary/secondary/tertiary amines (only conditionally cationic,
protonated below their pKa -- this ambiguity is disclosed explicitly rather
than resolved by guessing, matching this project's standing discipline of
refusing to silently assume unstated information).

A molecule matching both an anionic and a permanently-cationic (quaternary
ammonium) group is classified zwitterionic. A molecule matching neither is
classified nonionic -- but this is a claim by absence, so the low-confidence
case (a completely unfamiliar backbone with no recognized group at all) is
flagged, not silently trusted.
"""
from __future__ import annotations

from dataclasses import dataclass, field

try:
    from rdkit import Chem
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "classify.py requires rdkit (pip install rdkit) -- this is the one "
        "module in SurfactantKit with an external dependency, since correct "
        "functional-group detection from arbitrary real SMILES genuinely "
        "needs a real cheminformatics substructure-matching engine, not "
        "hand-rolled string parsing (which would be fragile against atom "
        "ordering/branching in real SMILES)."
    ) from e


# --- SMARTS patterns, one per real functional group that determines charge
# behavior in water. Each is a (name, smarts, charge_contribution, strength)
# tuple. "strength": "strong" = permanently/always ionized at normal aqueous
# pH (sulfate, sulfonate, quaternary ammonium); "moderate" = ionized at
# typical surfactant-testing pH (>~pH 6) but not at low pH (carboxylic
# acid/carboxylate, pKa ~4-5); "conditional" = only ionized below the
# group's own pKa, genuinely pH-dependent, disclosed as an open question
# rather than resolved (free amines, pKa ~9-11 so USUALLY protonated at
# neutral pH, but this project does not silently assume that).
_ANIONIC_PATTERNS = [
    ("sulfate_ester", "[#6]-[OX2]-[SX4](=[OX1])(=[OX1])-[OX1,OX2H0-1,OX2H1]", "strong"),
    ("sulfonate", "[#6]-[SX4](=[OX1])(=[OX1])-[OX1,OX2H0-1,OX2H1]", "strong"),
    ("phosphate_ester", "[#6]-[OX2]-[PX4](=[OX1])([OX1,OX2H,OX2H0-1])[OX1,OX2H,OX2H0-1]", "strong"),
    ("carboxylate", "[CX3](=[OX1])[OX1H0-1,OX2H1]", "moderate"),
]
_CATIONIC_STRONG_PATTERNS = [
    ("quaternary_ammonium", "[#6][NX4+]([#6])([#6])[#6]", "strong"),
    ("quaternary_ammonium_ring", "[#6][NX4+;R]", "strong"),
]
_CATIONIC_CONDITIONAL_PATTERNS = [
    ("primary_amine", "[NX3H2;!$(NC=O);!$(N=*)]", "conditional"),
    ("secondary_amine", "[NX3H1;!$(NC=O);!$(N=*);!R]", "conditional"),
    ("tertiary_amine", "[NX3H0;!$(NC=O);!$(N=*);!$([N+])]", "conditional"),
]

_COMPILED_ANIONIC = [(name, Chem.MolFromSmarts(sm), strength) for name, sm, strength in _ANIONIC_PATTERNS]
_COMPILED_CATIONIC_STRONG = [(name, Chem.MolFromSmarts(sm), strength) for name, sm, strength in _CATIONIC_STRONG_PATTERNS]
_COMPILED_CATIONIC_COND = [(name, Chem.MolFromSmarts(sm), strength) for name, sm, strength in _CATIONIC_CONDITIONAL_PATTERNS]


@dataclass
class ChargeClassificationResult:
    smiles: str
    charge_type: str  # "anionic" | "cationic" | "zwitterionic" | "nonionic" | "ambiguous_pH_dependent" | "unparseable"
    system_type_for_gibbs: str | None  # "ionic_no_added_salt" | "ionic_excess_electrolyte" | "nonionic" | None if undetermined
    anionic_groups_found: list[str] = field(default_factory=list)
    cationic_strong_groups_found: list[str] = field(default_factory=list)
    cationic_conditional_groups_found: list[str] = field(default_factory=list)
    confidence: str = ""  # "high" | "moderate" | "low" -- honestly stated, not hidden
    caveats: list[str] = field(default_factory=list)


def classify_surfactant_charge_type(smiles: str) -> ChargeClassificationResult:
    """Determine a surfactant's charge type (anionic/cationic/zwitterionic/
    nonionic) from its real SMILES structure via functional-group matching.

    Does NOT guess when the structure is genuinely ambiguous (e.g. a free
    amine with no quaternary nitrogen, whose charge state depends on a pH
    this function is never given) -- returns "ambiguous_pH_dependent" with
    the caveat spelled out, rather than silently picking cationic or
    nonionic, matching this project's standing "do not guess system_type"
    discipline already applied throughout adsorption.py/hld.py.

    system_type_for_gibbs is provided as a convenience default mapping onto
    the exact string this project's other functions (gibbs_gamma_max,
    szyszkowski_fit_K, ...) already expect -- but note it defaults ionic
    surfactants to "ionic_no_added_salt" (Gibbs prefactor n=2), since a bare
    SMILES carries no information about whether excess inert electrolyte is
    present in the actual solution being characterized; the caller must
    override this to "ionic_excess_electrolyte" if that is independently
    known to be true for their real system.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ChargeClassificationResult(
            smiles=smiles, charge_type="unparseable", system_type_for_gibbs=None,
            confidence="none",
            caveats=[f"RDKit could not parse this SMILES: '{smiles}'. Check for typos/invalid valence."],
        )

    anionic_found = [name for name, patt, _ in _COMPILED_ANIONIC if patt is not None and mol.HasSubstructMatch(patt)]
    cationic_strong_found = [name for name, patt, _ in _COMPILED_CATIONIC_STRONG if patt is not None and mol.HasSubstructMatch(patt)]
    cationic_cond_found = [name for name, patt, _ in _COMPILED_CATIONIC_COND if patt is not None and mol.HasSubstructMatch(patt)]

    is_anionic = len(anionic_found) > 0
    is_cationic_strong = len(cationic_strong_found) > 0
    is_cationic_conditional_only = (not is_cationic_strong) and len(cationic_cond_found) > 0

    caveats: list[str] = []

    if is_anionic and is_cationic_strong:
        return ChargeClassificationResult(
            smiles=smiles, charge_type="zwitterionic", system_type_for_gibbs=None,
            anionic_groups_found=anionic_found, cationic_strong_groups_found=cationic_strong_found,
            cationic_conditional_groups_found=cationic_cond_found, confidence="high",
            caveats=["Zwitterionic surfactants are net-neutral overall; the Gibbs prefactor "
                     "convention for zwitterionics is a separate real question this function "
                     "does not resolve -- do not default to either the ionic or nonionic "
                     "prefactor without checking the specific system's own literature "
                     "precedent (e.g. this project's own sulfobetaine/betaine validation "
                     "entries)."],
        )

    if is_anionic and is_cationic_conditional_only:
        return ChargeClassificationResult(
            smiles=smiles, charge_type="ambiguous_pH_dependent", system_type_for_gibbs=None,
            anionic_groups_found=anionic_found, cationic_conditional_groups_found=cationic_cond_found,
            confidence="low",
            caveats=[f"Contains both a real anionic group ({anionic_found}) and a free amine "
                     f"({cationic_cond_found}) with no stated pH -- this could be zwitterionic "
                     "(amino-acid-type), purely anionic (if the amine is unprotonated at the "
                     "real solution pH), or something else. Cannot be resolved from structure "
                     "alone; report the pH-dependence rather than guessing."],
        )

    if is_anionic:
        strengths = {name: s for name, _, s in _ANIONIC_PATTERNS if name in anionic_found}
        if all(strengths.get(n) == "moderate" for n in anionic_found):
            caveats.append("Only a carboxylate/carboxylic acid group found (pKa ~4-5) -- "
                            "this is ionized (anionic) at typical surfactant-testing pH (>~6) "
                            "but would behave as a nonionic (protonated) acid at low pH. "
                            "Confirm the real solution pH if operating near the pKa.")
        return ChargeClassificationResult(
            smiles=smiles, charge_type="anionic", system_type_for_gibbs="ionic_no_added_salt",
            anionic_groups_found=anionic_found, confidence="high" if "moderate" not in strengths.values() else "moderate",
            caveats=caveats,
        )

    if is_cationic_strong:
        return ChargeClassificationResult(
            smiles=smiles, charge_type="cationic", system_type_for_gibbs="ionic_no_added_salt",
            cationic_strong_groups_found=cationic_strong_found, confidence="high", caveats=caveats,
        )

    if is_cationic_conditional_only:
        return ChargeClassificationResult(
            smiles=smiles, charge_type="ambiguous_pH_dependent", system_type_for_gibbs=None,
            cationic_conditional_groups_found=cationic_cond_found, confidence="low",
            caveats=[f"Only a free amine ({cationic_cond_found}, typical pKa ~9-11) found, no "
                     "quaternary (permanently charged) nitrogen. This is USUALLY protonated "
                     "(cationic) at neutral/acidic pH, but this function does not assume that "
                     "-- report as pH-dependent rather than silently picking cationic, matching "
                     "this project's standing 'do not guess system_type' discipline."],
        )

    # No ionizable/charged group matched at all.
    return ChargeClassificationResult(
        smiles=smiles, charge_type="nonionic", system_type_for_gibbs="nonionic",
        confidence="moderate",
        caveats=["No anionic or cationic functional group was matched by this function's "
                 "known patterns. This is a claim by absence, not a positive structural "
                 "confirmation of a nonionic headgroup (e.g. PEG/glycoside) -- correct for "
                 "all real nonionic surfactants tested so far, but if this SMILES contains an "
                 "unusual/unrecognized ionizable group this function's pattern list does not "
                 "cover, it would be silently missed. Flagged as moderate, not high, "
                 "confidence for exactly this reason."],
    )
