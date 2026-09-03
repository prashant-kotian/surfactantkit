"""Category N: gemini-specific structural reasoning (30 questions).

Ties directly to the PhD's own domain (amidoamine-derived gemini cationic
surfactants). Not tool-mapped -- structural/architectural reasoning about
WHY a design choice matters, not a calculation. Gold answers are short key
phrases (per the lesson learned building Category M: long sentence-style gold
answers break the category_match substring grader), full rationale lives in
source_note.
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
    return next_id("N", i)

QUESTIONS: list[Question] = []

def add(subcat, text, gold, difficulty="medium", source=""):
    QUESTIONS.append(Question(
        id=qid(), category="N", subcategory=subcat, tools_required=[],
        difficulty=difficulty, question_text=text, given_data={},
        gold_answer=gold, grading_method="category_match", source_note=source,
    ))

# --- spacer length/type effects (10) ---
add("spacer_length", "Between a gemini surfactant with a C2 (ethylene) spacer and one with a C6 (hexamethylene) spacer, same tail length, which generally has the LOWER CMC?",
    "C6", source="CMC typically decreases as spacer length increases from very short toward an intermediate optimum (often around C6-C10), before rising again for very long/flexible spacers -- a well-established non-monotonic trend in gemini surfactant literature.")
add("spacer_length", "For very long, flexible alkyl spacers (e.g. C12 or longer), gemini surfactant behavior starts to resemble what, as the spacer becomes long enough to fold back into the micelle core?",
    "two independent monomeric surfactants", difficulty="hard",
    source="At very long spacer lengths, the spacer itself can fold into the hydrophobic core rather than remaining extended between headgroups, causing the gemini molecule's effective behavior to converge toward that of two separate monomeric surfactant molecules rather than a tightly coupled dimer.")
add("spacer_length", "A rigid spacer (e.g. containing an aromatic or stilbene unit) compared to a flexible alkyl spacer of similar length generally does what to the range of accessible headgroup spacings?",
    "restricts it", source="Rigidity prevents the spacer from bending/folding to adjust headgroup separation, locking in a narrower range of accessible geometries compared to a flexible chain.")
add("spacer_length", "A polar/hydrophilic spacer (e.g. containing ether oxygens) tends to sit closer to which region of a micelle, compared to a purely hydrophobic alkyl spacer?",
    "aqueous interface", source="A polar spacer is thermodynamically favored to remain near the water-facing interface rather than folding into the hydrophobic core, unlike a purely alkyl spacer, which may fold into the core especially at longer lengths.")
add("spacer_length", "Does increasing spacer length generally increase or decrease the effective headgroup area a0 relevant to the critical packing parameter, for short-to-intermediate spacer lengths?",
    "increase", source="A longer spacer generally holds the two headgroups farther apart, increasing the effective area they jointly occupy at the interface, which lowers CPP and favors more curved (spherical/cylindrical) aggregate morphologies.")
add("spacer_length", "What is the general effect of spacer length on a gemini surfactant's surface activity (efficiency at lowering surface tension), independent of its effect on CMC?",
    "also depends on spacer length", difficulty="hard",
    source="Spacer length affects packing efficiency at the interface (surface activity/pi_CMC) somewhat independently of its effect on aggregate formation (CMC) -- the two properties don't always track together, which is part of why morphology and interfacial performance are described in this project's own review as carrying more design information than CMC alone.")
add("spacer_length", "Two gemini surfactants have identical tail length and identical spacer LENGTH, but one has a flexible alkyl spacer and the other a rigid aromatic spacer. Are their CMCs guaranteed to be similar?",
    "no", source="Spacer chemistry (rigidity, polarity), not just spacer carbon count, materially affects packing and CMC -- two spacers of the same nominal length are not interchangeable if their chemical nature differs.")
add("spacer_length", "In gemini surfactant nomenclature, does a SHORTER spacer number (e.g. 12-2-12) or a LONGER spacer number (e.g. 12-12-12) generally correspond to tighter headgroup packing at the interface?",
    "shorter", source="Shorter spacers hold the two headgroups closer together, generally giving tighter/denser headgroup packing at the interface, all else equal.")
add("spacer_length", "For amidoamine-derived gemini surfactants specifically, what is the most common spacer length range used in the baseline synthetic literature this project's review paper surveys?",
    "C3 or C6", source="Directly consistent with this project's own review paper's baseline platform, which surveys amidoamine gemini surfactants primarily using C3 and C6 alkyl spacers.")
add("spacer_length", "Why does spacer length affect micelle AGGREGATION NUMBER, not just CMC?",
    "changes effective headgroup area and packing geometry", difficulty="hard",
    source="Since aggregation number depends on the packing geometry (via CPP, which itself depends on effective headgroup area), and spacer length directly modulates effective headgroup area, spacer length indirectly but genuinely affects the predicted/observed aggregation number, not only the CMC.")

# --- monomeric vs gemini vs bolaform comparisons (10) ---
add("architecture_comparison", "Compared to a monomeric surfactant with the same single tail length, does a gemini surfactant (two such tails, linked) generally have a HIGHER or LOWER CMC?",
    "lower", source="Well-established general trend -- gemini surfactants typically show CMCs one to two orders of magnitude lower than their monomeric analogues.")
add("architecture_comparison", "Compared to a gemini surfactant, does a bolaform surfactant (one long chain, headgroup at each end) generally form MICELLES of similar aggregation number, or does it tend toward different (e.g. more unusual/smaller) aggregate types?",
    "different aggregate types", difficulty="hard",
    source="Bolaform surfactants, having only one hydrophobic chain that must span between two headgroups (rather than two independent chains able to pack into a core), often show distinct self-assembly behavior -- e.g. they can be less prone to forming conventional micelles and more prone to monolayer or membrane-spanning arrangements -- compared to gemini surfactants.")
add("architecture_comparison", "For the SAME total hydrophobic tail volume, does a gemini surfactant (2 tails, half the length each) or a monomeric surfactant (1 tail, full length) generally pack into a smaller critical packing parameter regime, favoring more curved (spherical) structures?",
    "gemini", difficulty="hard",
    source="Splitting the same total tail volume into two shorter tails linked near two separate headgroups tends to increase the effective headgroup-area contribution relative to tail volume in typical gemini architectures, generally favoring lower CPP / more curved structures compared to one long single tail -- though this depends on spacer length and headgroup chemistry too.")
add("architecture_comparison", "A trimeric surfactant (three linked head-tail units) compared to a gemini (two-unit) surfactant of the same individual tail length generally shows what direction of change in CMC?",
    "further decrease", source="Extending the cooperative multi-tail effect from two to three linked units generally continues the trend of CMC lowering seen going from monomeric to gemini, though the literature on trimeric/higher-order surfactants is considerably sparser than for gemini systems.")
add("architecture_comparison", "Why is direct chain-length-resolved comparison between monomeric and gemini surfactant series considered valuable for a mixed-system review, according to this project's own review paper's stated design principle?",
    "isolates the architecture effect", source="Keeping tail chain length, spacer, and other variables matched while varying only monomeric-vs-gemini architecture isolates the specific effect of architecture itself, rather than confounding it with unrelated chain-length differences -- directly the design principle used in this project's own review (e.g. citing Ghumare et al.'s matched C8-C18 monomeric/gemini series).")
add("architecture_comparison", "Does gemini architecture ALONE (lower CMC) guarantee favorable (synergistic) mixing when combined with a second surfactant, according to this project's own review paper's stated finding?",
    "no", source="Directly stated as a key finding of this project's own review: gemini architecture lowers CMC but does not guarantee favorable mixing outcomes -- cationic-cationic gemini pairs specifically were found to swing between synergistic and antagonistic depending on spacer structure and packing fit.")
add("architecture_comparison", "What single structural feature most directly distinguishes a gemini surfactant from simply mixing two separate monomeric surfactant molecules in solution?",
    "covalent linkage between the two units", source="The defining feature of a gemini surfactant is the covalent spacer bond holding the two head-tail halves together as one molecule, as opposed to two chemically independent monomeric molecules that merely happen to co-exist in a mixture.")
add("architecture_comparison", "Does a gemini surfactant's covalent linkage make its two 'halves' free to diffuse and adsorb independently at the interface, the way two separate monomeric molecules in a mixed system would?",
    "no", source="The covalent spacer constrains the two halves to move and adsorb together as one unit, fundamentally different from an unlinked mixed system of two independent monomeric molecules, even if the two gemini halves are chemically identical to two monomeric molecules.")
add("architecture_comparison", "Amidoamine monomeric and gemini surfactants sharing the same fatty-acid-derived synthetic platform allow what kind of structure-property comparison that unrelated surfactant families cannot?",
    "genuine homologous comparison", source="Directly consistent with this project's own review paper's stated rationale -- a shared synthetic route lets chain length, architecture, and spacer be varied independently without changing the underlying chemistry, giving cleaner structure-property conclusions than comparing surfactants pulled from unrelated synthetic routes.")
add("architecture_comparison", "Why might chain-length-resolved comparisons between monomeric and gemini amidoamine surfactants still be described as a 'gap' in the current literature, according to this project's own review paper?",
    "few matched chain-length studies exist", difficulty="hard",
    source="Directly stated as a critical gap in this project's own review: most published amidoamine studies don't provide matched, systematically chain-length-varied comparisons across both monomeric and gemini architectures, limiting how confidently structure-property trends can be generalized.")

# --- amidoamine-specific structural reasoning (10) ---
add("amidoamine_structure", "What functional group links the fatty-acid-derived hydrophobic tail to the cationic headgroup in an amidoamine surfactant?",
    "amide", source="The defining structural feature of the amidoamine class -- an amide bond formed between a fatty acid and an aminoalkylamine, as opposed to a direct alkyl-ammonium bond (no linking functional group) found in simpler quaternary ammonium surfactants like CTAB.")
add("amidoamine_structure", "How does the amide linkage in amidoamine surfactants affect intermolecular hydrogen bonding compared to a simple alkyl-ammonium surfactant with no amide group?",
    "amide enables additional hydrogen bonding", source="The amide C=O and N-H groups can participate in hydrogen bonding between adjacent surfactant molecules, an interaction not available to a simple alkyl-quaternary-ammonium surfactant lacking any amide group -- directly tied to this project's own review paper's discussion of the amide linkage 'doing real chemical work.'")
add("amidoamine_structure", "What common synthetic precursor is used to convert a long-chain fatty acid into an amidoamine intermediate, in the baseline synthetic route surveyed by this project's own review paper?",
    "dimethylaminopropylamine", source="Directly stated in this project's own review paper as the standard aminoalkylamine reagent (or a related aminoalkylamine) used to convert fatty acids to amidoamine intermediates, prior to quaternization to the final cationic surfactant.")
add("amidoamine_structure", "After forming the amidoamine intermediate from a fatty acid and an aminoalkylamine, what reaction step converts it into the final CATIONIC surfactant?",
    "quaternization", source="The tertiary amine in the amidoamine intermediate is quaternized (typically via alkylation, e.g. with an alkyl halide) to install the permanent positive charge that makes the final molecule cationic.")
add("amidoamine_structure", "Does the shared amidoamine synthetic platform allow chain length, architecture (monomeric vs gemini), and spacer structure to be varied together only, or independently of each other?",
    "independently", source="Directly stated in this project's own review paper as the key advantage of this synthetic route -- chain length, architecture, and spacer can each be varied independently without touching the underlying amide-forming/quaternization chemistry.")
add("amidoamine_structure", "Amidobetaine surfactants, related to amidoamine surfactants, carry what type of net charge behavior, distinct from a simple amidoamine-derived cationic surfactant?",
    "zwitterionic", difficulty="hard",
    source="Amidobetaines add a carboxylate (or similar anionic) group onto the amidoamine-derived cationic headgroup, making the overall molecule zwitterionic (permanent dual charge) rather than purely cationic -- relevant to the amidobetaine oilfield studies discussed in this project's own review.")
add("amidoamine_structure", "Why might an amidoamine surfactant show different pH-dependent behavior than a simple quaternary-ammonium (permanently cationic) surfactant like CTAB?",
    "the amide/amine chemistry can be pH sensitive", difficulty="hard",
    source="Depending on the exact amidoamine structure (e.g. whether a tertiary amine remains unquaternized, or the amide itself is base/acid sensitive under extreme pH), amidoamine-derived surfactants can show pH-dependent behavior that a permanently quaternized simple alkylammonium surfactant like CTAB does not.")
add("amidoamine_structure", "What is the main practical/synthetic reason the amide linkage is chosen over a direct (non-amide) alkyl-to-headgroup bond in amidoamine surfactant synthesis?",
    "synthetic convenience from fatty acid feedstocks", source="Directly stated in this project's own review paper -- the amide-forming route offers synthetic convenience (fatty acids are readily available feedstocks reacted with an aminoalkylamine), not only the resulting hydrogen-bonding/property benefits.")
add("amidoamine_structure", "For an 'as-synthesized' internal mixture of amidoamine/amidobetaine surfactants (not purified to one exact structure), what commonly co-occurs alongside the intended amide product from the same synthesis batch?",
    "residual side products", difficulty="hard",
    source="Synthesis-derived internal mixtures, as discussed in this project's own review paper (e.g. the as-synthesized oleic amidopropyl betaine internal mixture), commonly contain residual starting material, unreacted amide, or related side-products alongside the intended surfactant, which can itself act as part of the real formulation's effective behavior.")
add("amidoamine_structure", "Why does an amidoamine gemini surfactant's amide linkage specifically (as opposed to its spacer) matter for interpreting its mixed-system behavior with other surfactants?",
    "amide hydrogen bonding effects", difficulty="hard",
    source="The amide group's hydrogen-bonding capacity is a separate structural variable from spacer length/rigidity -- both can independently affect how an amidoamine gemini surfactant interacts with a second surfactant in a mixed system, so attributing an observed mixing effect to 'spacer' alone without considering the amide linkage risks an incomplete structure-property conclusion.")

if __name__ == "__main__":
    print(f"Category N: {len(QUESTIONS)} questions")
    from collections import Counter
    print(Counter(q.subcategory for q in QUESTIONS))
    long_ones = [q.id for q in QUESTIONS if isinstance(q.gold_answer, str) and len(q.gold_answer.split()) > 8]
    print("gold answers over 8 words (should be empty):", long_ones)
