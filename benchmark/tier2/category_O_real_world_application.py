"""Category O: real-world application (26 questions).

Formulation/EOR/solubilization reasoning without a clean numeric answer --
judgment questions about which approach/design choice is appropriate for a
stated real-world goal. Not tool-mapped: this is applied judgment, not a
quantity a calculation tool produces. Gold answers are short key phrases
(same discipline established for M/N), full rationale in source_note.
"""

from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from schema import Question, next_id

i = 0
def qid():
    global i
    i += 1
    return next_id("O", i)

QUESTIONS: list[Question] = []

def add(subcat, text, gold, difficulty="medium", source=""):
    QUESTIONS.append(Question(
        id=qid(), category="O", subcategory=subcat, tools_required=[],
        difficulty=difficulty, question_text=text, given_data={},
        gold_answer=gold, grading_method="category_match", source_note=source,
    ))

# --- enhanced oil recovery / formulation judgment (9) ---
add("eor_judgment", "For an enhanced oil recovery formulation, is the PRIMARY lever for increasing the capillary number (and thus mobilizing trapped oil) generally increasing injection velocity, or reducing interfacial tension?",
    "reducing interfacial tension", source="Capillary number Ca = (viscosity*velocity)/IFT; IFT can be driven down by 3-5+ orders of magnitude with the right surfactant system, while velocity increases achievable in a real reservoir are comparatively minor -- IFT reduction is the standard, dominant practical lever in EOR surfactant design.")
add("eor_judgment", "For an EOR surfactant screening program targeting ultra-low interfacial tension, would a surfactant that performs well on pure model oil-brine systems but has NOT been tested with actual crude oil and reservoir brine be considered field-ready?",
    "no", source="Real crude oil and reservoir brine contain complex components (asphaltenes, multivalent ions, etc.) that can substantially change surfactant performance versus simplified model systems -- consistent with this project's own review paper's emphasis on studies using actual crude oil/seawater/elevated temperature as more representative than simplified model systems.")
add("eor_judgment", "If two surfactant formulations achieve similar ultra-low interfacial tension, but one requires a much higher surfactant concentration to do so, which is generally preferred for a cost-sensitive EOR field application?",
    "the lower-concentration formulation", source="Surfactant cost scales directly with concentration used at field scale; achieving the same IFT-reduction performance at lower dosage is a major practical/economic advantage in EOR, independent of the raw physical performance being similar.")
add("eor_judgment", "For an EOR application in a HIGH-salinity/high-temperature reservoir, would a simple anionic surfactant (e.g. plain SDS) generally be a robust first choice, given known salt/electrolyte sensitivity of ionic surfactant CMC and stability?",
    "no", source="High salinity can cause ionic surfactant precipitation/phase separation issues (e.g. Krafft-point and salinity-tolerance concerns), which is why more salt-tolerant formulations (e.g. specific betaines, co-surfactant systems, or nonionic/zwitterionic components) are generally preferred for harsh-condition reservoirs rather than a simple anionic surfactant alone.")
add("eor_judgment", "Does adding a co-solvent (e.g. a glycol ether) to a surfactant EOR formulation, as described in real literature validated in this project (Wazir et al.), act purely as an inert diluent, or does it meaningfully change the formulation's effective HLB/microemulsion behavior?",
    "changes effective HLB/microemulsion behavior", source="Directly consistent with this project's own review paper's discussion of Wazir et al.'s co-solvent work -- the co-solvent shifts the whole formulation's effective hydrophilic-lipophilic balance and changes microemulsion behavior outright, not acting as a passive diluent.")
add("eor_judgment", "For a wettability-alteration application on an oil-wet rock surface, would a surfactant needs to primarily REDUCE the contact angle (toward more water-wet) or INCREASE it, to improve water-flood oil recovery?",
    "reduce", source="Oil-wet rock resists water displacing oil; reducing the contact angle (shifting toward more water-wet behavior) improves water's ability to displace oil from the rock surface, the standard goal of wettability-alteration surfactant treatments.")
add("eor_judgment", "Between a surfactant with strongly negative spreading coefficient S and one with S close to zero on a given surface, which shows better spontaneous spreading/wetting behavior?",
    "S close to zero", source="S closer to zero (approaching the complete-wetting limit) indicates more favorable spontaneous spreading; strongly negative S indicates poor wetting -- directly following from the spreading_coefficient formula's sign convention already implemented in this library.")
add("eor_judgment", "Is achieving the theoretically lowest possible interfacial tension always the correct design goal for a practical EOR formulation, or must formulation stability/robustness under real reservoir conditions also be weighed?",
    "stability must also be weighed", source="An ultra-low-IFT formulation that phase-separates or degrades under real reservoir salinity/temperature is not practically useful even if its IFT number looks best in idealized lab conditions -- practical EOR design requires balancing performance against real-condition robustness, not optimizing one number in isolation.")
add("eor_judgment", "For a clay-stabilization application (preventing clay swelling/migration in a reservoir), would a CATIONIC or ANIONIC surfactant generally be more effective, given that most clay surfaces carry a net negative charge?",
    "cationic", source="Cationic surfactants electrostatically adsorb onto negatively-charged clay surfaces, which is the standard mechanism behind their use in clay stabilization applications, consistent with the poly(oxyethylene)-amidoamine gemini clay/bentonite application discussed in this project's own review paper.")

# --- solubilization / drug delivery judgment (7) ---
add("solubilization_judgment", "For solubilizing a hydrophobic drug or PAH compound, does increasing surfactant concentration WELL ABOVE the CMC generally increase or have negligible further effect on total solubilization capacity?",
    "increase", source="Above the CMC, additional surfactant forms more micelles, which can solubilize additional hydrophobic solute -- solubilization capacity scales with the amount of surfactant present above the CMC, which is exactly the physical basis of the Molar Solubilization Ratio (MSR) formula.")
add("solubilization_judgment", "Below the CMC, does adding more surfactant meaningfully increase a hydrophobic compound's solubility in water?",
    "no", source="Below the CMC there are no micelles to solubilize into; solubility below the CMC is governed by the compound's intrinsic (micelle-independent) water solubility, essentially unaffected by surfactant concentration in that regime -- the physical basis for this library's molar_solubilization_ratio guard rail requiring surfactant concentration above the CMC.")
add("solubilization_judgment", "Would a gemini surfactant (typically lower CMC, denser packing) generally be expected to be MORE or LESS efficient at solubilizing a hydrophobic compound per mole of surfactant, compared to its monomeric analogue, all else equal?",
    "more efficient", difficulty="hard",
    source="Denser interfacial/micellar packing and larger effective hydrophobic domain per aggregate generally correlates with enhanced solubilization capacity for gemini surfactants relative to monomeric analogues, though the exact magnitude is system-dependent.")
add("solubilization_judgment", "For a drug-delivery application requiring biocompatibility, would a synthetic quaternary ammonium cationic surfactant or a bile-salt-derived (e.g. sodium cholate) surfactant generally be preferred, all else being comparable?",
    "bile-salt-derived", source="Bile salts are naturally occurring, generally recognized as more biocompatible than synthetic quaternary ammonium surfactants, which can show higher cytotoxicity/membrane-disruptive effects -- a standard consideration in pharmaceutical surfactant selection.")
add("solubilization_judgment", "If a solubilization study reports a solute concentration in the aqueous phase that is LESS than the solute's known intrinsic water solubility limit, even in the presence of surfactant, does that necessarily indicate the surfactant is solubilizing effectively above baseline?",
    "no", difficulty="hard",
    source="If the reported concentration doesn't even reach the pure-water saturation limit, the system may be operating in a sub-saturation partitioning regime rather than demonstrating surfactant-enhanced solubilization beyond water's own baseline capacity -- directly the real methodological subtlety this project's own literature validation work uncovered when trying to test molar_solubilization_ratio against real rhamnolipid-phenol data.")
add("solubilization_judgment", "For encapsulating a very hydrophobic compound (very low intrinsic water solubility), would surfactant micellization generally provide a LARGER or SMALLER relative (fold) increase in apparent solubility, compared to a moderately hydrophobic compound?",
    "larger", source="The MSR-driven enhancement is added on top of a very small baseline (intrinsic water solubility), so the same absolute micellar solubilization contribution represents a much larger relative/fold increase for a very hydrophobic compound than for a more water-soluble one.")
add("solubilization_judgment", "Would using a mixed surfactant system with a demonstrated synergistic (negative beta) interaction generally be expected to improve or worsen solubilization capacity, compared to either pure surfactant alone at the same total concentration?",
    "improve", difficulty="hard",
    source="Synergistic mixing generally produces a lower mixed CMC and often denser/more effective micellar packing, both of which tend to favor greater solubilization capacity per unit surfactant compared to either single component -- though this isn't guaranteed in every system and would ideally be confirmed experimentally.")

# --- general formulation reasoning (10) ---
add("formulation_judgment", "For a household detergent formulation prioritizing low cost and biodegradability, would a synthetic gemini surfactant or a simple, well-established anionic surfactant (e.g. LAS, linear alkylbenzene sulfonate) typically be the default industrial choice?",
    "simple anionic surfactant", source="Cost and established biodegradability profiles (plus decades of regulatory/manufacturing infrastructure) make simple, well-characterized anionic surfactants the practical default for high-volume, cost-sensitive detergent applications, despite gemini surfactants' superior CMC/interfacial performance in lab settings.")
add("formulation_judgment", "If a formulator needs to STABILIZE a foam (e.g. in a firefighting or personal-care foam application), would they generally prioritize a surfactant with high or low surface elasticity/viscoelastic monolayer behavior?",
    "high", source="Viscoelastic, mechanically robust adsorbed surfactant monolayers resist film thinning and drainage, which is the standard mechanism behind foam stabilization -- surfactants forming rigid/viscoelastic films at the interface are preferred for foam-stability-critical applications.")
add("formulation_judgment", "For an emulsion formulation requiring a stable oil-in-water emulsion, would a surfactant with a Griffin HLB value around 4 (low, lipophilic) or around 12-16 (higher, hydrophilic) generally be preferred?",
    "12-16", source="Griffin's HLB scale guideline: low HLB (~3-6) favors water-in-oil emulsions; higher HLB (~8-18, commonly cited 12-16 as a practical o/w range) favors oil-in-water emulsions -- standard HLB-based emulsifier selection heuristic.")
add("formulation_judgment", "For a water-in-oil emulsion instead, would a LOW or HIGH Griffin HLB surfactant generally be preferred?",
    "low", source="The inverse of the o/w guideline -- low-HLB (lipophilic-dominant) surfactants preferentially stabilize water-in-oil emulsions.")
add("formulation_judgment", "If a single surfactant cannot achieve a target HLB value exactly, would blending two surfactants with HLB values bracketing the target (one above, one below) be a standard, accepted formulation strategy?",
    "yes", source="Blending surfactants of different individual HLB to hit an intermediate target HLB (weighted average) is a standard, widely-used formulation technique, precisely because pure single-surfactant HLB values are fixed and rarely match an arbitrary target exactly.")
add("formulation_judgment", "Would a formulator generally expect a nonionic surfactant's CMC and micellization behavior to be strongly or weakly sensitive to solution ionic strength/added salt, compared to an ionic surfactant?",
    "weakly", source="Nonionic surfactants lack charged headgroups, so they aren't subject to the electrostatic-screening mechanism (Debye screening, counterion effects) that makes ionic surfactant CMC strongly salt-sensitive -- a standard, well-established distinction.")
add("formulation_judgment", "For a formulation intended to work reliably across a WIDE pH range (e.g. both acidic and basic conditions), would a simple amphoteric surfactant (whose charge state changes with pH) or a permanently-charged surfactant (e.g. quaternary ammonium cationic) generally show more consistent behavior?",
    "permanently-charged surfactant", difficulty="hard",
    source="An amphoteric surfactant's net charge (and therefore CMC, interfacial behavior, and compatibility with other charged components) changes across the pH range, making its behavior inherently more pH-dependent/variable, whereas a permanently cationic (or anionic) surfactant's charge state doesn't change with pH -- more consistent, if less versatile in other respects.")
add("formulation_judgment", "If a formulation needs to remain stable and effective across a wide temperature range (e.g. a product stored/used from cold to hot climates), is CMC's known temperature dependence (per this project's own real literature findings) a factor worth considering in surfactant selection?",
    "yes", source="Directly consistent with this project's own literature-validated finding that real surfactant CMC values shift measurably with temperature (e.g. [C12mpy][Br] CMC rising from 9.74 to 11.47 mM across a 20K range) -- formulators do need to consider this if wide-temperature-range performance matters.")
add("formulation_judgment", "For a formulation combining a surfactant with a polymer (e.g. a Pluronic block copolymer) for enhanced viscosity/stability, would the resulting mixed-system behavior generally be assumed identical to a simple two-surfactant mixture, or does polymer-surfactant coassembly typically introduce additional/different physics?",
    "additional/different physics", difficulty="hard",
    source="Polymer-surfactant systems (e.g. gemini-Pluronic coassembly, discussed in this project's own review paper) involve polymer-specific interactions (e.g. surfactant binding along polymer chains, polymer-induced micelle bridging) not captured by simple two-surfactant regular-solution theory (Rubingh/Clint) alone -- a genuinely distinct physical regime, not just 'another mixed surfactant pair.'")
add("formulation_judgment", "Would a formulator generally expect a stimuli-responsive emulsion gel system (e.g. one incorporating a pH- or temperature-responsive amidoamine surfactant) to have MORE or FEWER formulation variables to control than a simple fixed-composition surfactant system?",
    "more", source="Stimuli-responsiveness inherently adds design variables (the trigger condition, the response magnitude/reversibility, compatibility of the responsive mechanism with the rest of the formulation) beyond a simple fixed-composition system, consistent with this project's own review paper's discussion of stimuli-responsive emulsion gel systems as a more complex application category.")

if __name__ == "__main__":
    print(f"Category O: {len(QUESTIONS)} questions")
    from collections import Counter
    print(Counter(q.subcategory for q in QUESTIONS))
    long_ones = [q.id for q in QUESTIONS if isinstance(q.gold_answer, str) and len(q.gold_answer.split()) > 8]
    print("gold answers over 8 words (should be empty):", long_ones)
