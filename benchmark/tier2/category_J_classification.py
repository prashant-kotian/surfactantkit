"""Category J: classification & nomenclature (40 questions).

Not tool-mapped by design (per taxonomy.md) -- no SurfactantKit calculation
applies to "what class of surfactant is this." Hand-curated, grounded in
standard, well-established surfactant-science definitions (Rosen's
"Surfactants and Interfacial Phenomena" and equivalent standard texts), not
obscure or contested classifications. Every fact here is the kind stated
identically across essentially all colloid-science textbooks -- deliberately
avoiding borderline/contested classification calls, since a benchmark
question needs one unambiguous right answer.

grading_method="category_match": the model's answer just needs to contain the
correct classification keyword (checked by the same category_match grading
already used for Tier 1 no-solution traps), not match verbatim phrasing.
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
    return next_id("J", i)

QUESTIONS: list[Question] = []

def add(subcat, text, gold, difficulty="easy", source=""):
    QUESTIONS.append(Question(
        id=qid(), category="J", subcategory=subcat, tools_required=[],
        difficulty=difficulty, question_text=text, given_data={},
        gold_answer=gold, grading_method="category_match", source_note=source,
    ))

# --- ionic-class identification from headgroup description (12) ---
add("ionic_class", "A surfactant headgroup is -SO4^- Na+. What ionic class does this surfactant belong to?",
    "anionic", source="Sulfate headgroup carries a negative charge -> anionic surfactant (standard definition, e.g. SDS).")
add("ionic_class", "A surfactant headgroup is -N(CH3)3^+ Br-. What ionic class does this surfactant belong to?",
    "cationic", source="Quaternary ammonium headgroup carries a positive charge -> cationic surfactant (standard definition, e.g. CTAB).")
add("ionic_class", "A surfactant headgroup carries both a permanent positive charge (quaternary ammonium) and a permanent negative charge (sulfonate) on the same molecule, at all pH. What ionic class is this?",
    "zwitterionic", source="Permanent, pH-independent dual charge on one molecule is the defining feature of zwitterionic (not just 'amphoteric') surfactants.")
add("ionic_class", "A surfactant headgroup is a carboxylic acid group that is protonated (neutral) at low pH and deprotonated (negatively charged) at high pH. What class of surfactant is this, with respect to pH dependence?",
    "amphoteric", difficulty="medium",
    source="Amphoteric surfactants change net charge with pH (positive/neutral/negative depending on pH), distinct from zwitterionic, which carries both charges simultaneously and permanently.")
add("ionic_class", "A surfactant has a polyoxyethylene (PEG) headgroup with no ionizable or charged group at any pH. What ionic class is this?",
    "nonionic", source="No charge at any pH is the defining feature of nonionic surfactants (e.g. Triton X-100, Tween 80).")
add("ionic_class", "Sodium dodecyl sulfate (SDS) belongs to which ionic class of surfactant?",
    "anionic", source="Standard textbook classification -- sulfate headgroup, sodium counterion.")
add("ionic_class", "Cetyltrimethylammonium bromide (CTAB) belongs to which ionic class of surfactant?",
    "cationic", source="Standard textbook classification -- quaternary ammonium headgroup.")
add("ionic_class", "Triton X-100 belongs to which ionic class of surfactant?",
    "nonionic", source="Standard textbook classification -- polyoxyethylene headgroup, no charge.")
add("ionic_class", "Cocamidopropyl betaine belongs to which ionic class of surfactant?",
    "zwitterionic", source="Carries a permanent quaternary ammonium cation and a carboxylate anion simultaneously -- the defining 'betaine' structure.")
add("ionic_class", "Sodium cholate (a bile salt) belongs to which ionic class of surfactant?",
    "anionic", source="Carboxylate headgroup, ionized at physiological/most working pH.")
add("ionic_class", "Docusate sodium (AOT) belongs to which ionic class of surfactant?",
    "anionic", source="Sulfosuccinate headgroup, standard anionic classification.")
add("ionic_class", "An amino-acid-derived surfactant with a free carboxylic acid group that is anionic at high pH and cationic (protonated amine) at low pH is classified as what type, with respect to pH-dependent charge?",
    "amphoteric", difficulty="medium", source="Charge sign changes with pH, matching the amphoteric definition.")

# --- architectural classification: monomeric vs gemini vs bolaform (10) ---
add("architecture_class", "A surfactant molecule has ONE hydrophobic tail and ONE polar headgroup. What architectural class is this?",
    "monomeric", source="Single tail, single headgroup is the baseline 'conventional'/monomeric surfactant architecture.")
add("architecture_class", "A surfactant molecule has TWO hydrophobic tails and TWO polar headgroups, covalently linked near the headgroups by a spacer. What architectural class is this?",
    "gemini", source="Two head-tail units joined near the headgroups by a spacer is the defining gemini (dimeric) surfactant architecture.")
add("architecture_class", "A surfactant molecule has a single long hydrophobic chain with a polar/ionic headgroup at BOTH ends (no headgroup in the middle). What architectural class is this?",
    "bolaform", difficulty="medium",
    source="Two headgroups at opposite ends of one hydrophobic chain is the defining bolaform (bola-amphiphile) architecture -- distinct from gemini, which has two full head-tail units, not one chain with two end-groups.")
add("architecture_class", "12-6-12 is standard shorthand naming a surfactant with two C12 tails, two headgroups, and a hexamethylene (C6) spacer linking them. What architectural class does this name describe?",
    "gemini", source="The 'm-s-m' gemini naming convention (tail carbons - spacer carbons - tail carbons).")
add("architecture_class", "In gemini surfactant nomenclature '12-3-12', what does the middle number (3) represent?",
    "spacer length", difficulty="medium", source="Standard gemini 'm-s-m' naming: m=tail carbon number, s=spacer carbon number.")
add("architecture_class", "Compared to its monomeric analogue with the same tail length, a gemini surfactant's CMC is typically (higher / lower / the same order of magnitude but usually lower)?",
    "lower", source="Well-established general trend: gemini surfactants typically have CMCs one to two orders of magnitude lower than their monomeric analogues, due to two tails cooperating at the interface/in the micelle core.")
add("architecture_class", "A surfactant with three or more covalently linked head-tail units (rather than two) is generally termed what, relative to a gemini (two-unit) surfactant?",
    "trimeric", difficulty="medium", source="Extending the gemini (dimeric, two-unit) concept to three linked units gives 'trimeric' surfactants -- a recognized, if less common, architectural class in the literature.")
add("architecture_class", "What is the primary structural role of the 'spacer' group in a gemini surfactant?",
    "links the two headgroups", source="The spacer covalently connects the two head-tail halves near the headgroups, and its length/rigidity/chemistry controls headgroup spacing and packing.")
add("architecture_class", "A surfactant synthesized with a mixture of chain lengths and minor structural variants from a single reaction batch (not purified to one exact structure) is generally described in the literature as what kind of mixture?",
    "internal mixture", difficulty="medium",
    source="Synthesis-derived internal/technical mixtures are a recognized category distinct from a deliberately-blended two-component mixed system -- relevant to amidoamine synthesis routes discussed in this project's own review paper.")
add("architecture_class", "Amidoamine-derived cationic surfactants are typically synthesized from a fatty acid via what class of chemical linkage connecting the hydrophobic tail to the cationic headgroup?",
    "amide", source="The defining feature of amidoamine surfactants -- amide bond links the fatty-acid-derived tail to the aminoamine-derived headgroup, as opposed to a direct alkyl-ammonium (no linker) or ester linkage.")

# --- micelle/aggregate morphology terms (10) ---
add("morphology_term", "A surfactant aggregate that is roughly spherical, with the hydrophobic tails forming an inner core and headgroups facing the aqueous exterior, is called what?",
    "spherical micelle", source="Standard, most common micelle morphology for CPP <= 1/3.")
add("morphology_term", "A surfactant aggregate that has grown elongated, roughly cylindrical, is called what (two acceptable synonymous terms)?",
    "wormlike micelle", source="Also called 'threadlike' or 'rodlike' micelle depending on aspect ratio; 'wormlike' specifically implies long, flexible, entangled cylinders capable of imparting viscoelasticity.")
add("morphology_term", "A closed bilayer shell of surfactant enclosing an aqueous interior is called what?",
    "vesicle", source="Standard term for a bilayer-enclosed aqueous compartment, distinct from a solid-core micelle.")
add("morphology_term", "An extended, essentially infinite (relative to molecular scale) planar surfactant bilayer sheet is called what?",
    "bilayer", source="The lamellar/bilayer phase, the CPP~1 limiting morphology.")
add("morphology_term", "A surfactant aggregate in a nonpolar/oil continuous phase, with polar headgroups facing inward around a small aqueous or polar core, is called what (as opposed to a normal micelle in water)?",
    "reverse micelle", source="Also 'inverted micelle' -- headgroups face inward, tails face the nonpolar continuous phase, opposite orientation to a normal aqueous micelle.")
add("morphology_term", "Long, entangled wormlike micelles that impart significant viscosity/viscoelasticity to a solution, similar to polymer solutions, are often described using what polymer-physics-borrowed term?",
    "viscoelastic", difficulty="medium", source="Wormlike micelle solutions are frequently described as 'viscoelastic' or forming a 'micellar network,' borrowing directly from polymer rheology language.")
add("morphology_term", "What is the term for the surfactant concentration at which the FIRST aggregates (typically spherical micelles) begin to form in solution?",
    "critical micelle concentration", source="CMC -- the standard, most fundamental defined threshold concentration in surfactant science.")
add("morphology_term", "In DLS/SANS/SAXS characterization, the parameter describing how elongated/non-spherical a particle or aggregate is (ratio-based shape descriptor) is called what?",
    "asphericity", difficulty="medium", source="Or equivalently 'aspect ratio' / 'eccentricity' depending on the exact measure -- standard shape-characterization terminology.")
add("morphology_term", "The transition of surfactant aggregates from spherical micelles to wormlike/cylindrical micelles as CPP increases past 1/3 is generally called what kind of transition?",
    "morphology transition", source="Or 'shape transition' -- standard terminology in the self-assembly literature for CPP-driven aggregate shape changes.")
add("morphology_term", "What is the standard term for the number of surfactant monomers making up a single micelle?",
    "aggregation number", source="Standard defined quantity (often denoted N_agg), directly computed by SurfactantKit's aggregation_number tool from geometric packing.")

# --- gemini-specific spacer/architecture terms, tied to the PhD's own domain (8) ---
add("gemini_spacer", "A gemini surfactant spacer made of a flexible alkyl chain (e.g. -(CH2)n-) is generally classified as what type of spacer, as opposed to a rigid aromatic spacer?",
    "flexible", source="Standard spacer classification: flexible (alkyl) vs rigid (e.g. stilbene, aromatic-containing) vs polar/hydrophilic (e.g. ethoxylated) spacers, each affecting headgroup spacing and packing differently.")
add("gemini_spacer", "As a gemini surfactant's alkyl spacer length increases from very short (e.g. C2-C3) toward intermediate length (e.g. C6), the CMC typically does what?",
    "decreases", difficulty="medium",
    source="Well-established non-monotonic trend: CMC generally decreases as spacer length increases from very short spacers, reaches a minimum around intermediate spacer length (~C6-C10 depending on tail length), then can increase again for very long/flexible spacers that start behaving more like two separate monomeric surfactants.")
add("gemini_spacer", "A gemini surfactant spacer containing ether oxygen linkages (e.g. -(CH2CH2O)n-) rather than plain alkyl carbons is generally classified as what type of spacer?",
    "polar/hydrophilic", source="Ethoxylated or otherwise heteroatom-containing spacers are classified as polar/hydrophilic spacers, distinct from purely hydrophobic alkyl spacers -- they can sit closer to the aqueous interface rather than folding into the micelle core.")
add("gemini_spacer", "A rigid, non-flexible gemini spacer (e.g. containing a stilbene or other aromatic/unsaturated rigid unit) primarily restricts what molecular behavior compared to a flexible alkyl spacer?",
    "conformational flexibility", difficulty="medium",
    source="Rigid spacers restrict the spacer's ability to fold/bend, which in turn constrains headgroup spacing and micelle packing geometry differently than a flexible spacer can.")
add("gemini_spacer", "What is the standard term for a gemini surfactant in which the two tail-headgroup halves are chemically identical (same tail length, same headgroup)?",
    "symmetric gemini", source="Standard terminology distinguishing symmetric (identical halves) from asymmetric/heterogemini surfactants (different tail lengths and/or different headgroups on each half).")
add("gemini_spacer", "What is the standard term for a gemini surfactant in which the two halves have different tail lengths and/or different headgroup chemistries?",
    "asymmetric gemini", source="Also called 'heterogemini' -- standard terminology for non-identical gemini halves.")
add("gemini_spacer", "Amidoamine-derived gemini surfactants typically link the two amidoamine head-tail halves via what type of spacer group in the most common synthetic route?",
    "alkyl", difficulty="medium",
    source="The common synthetic route links two amidoamine units via a flexible alkyl diamine or dihalide spacer (e.g. C3 or C6), consistent with this project's own review paper's baseline synthetic platform.")
add("gemini_spacer", "For a fixed tail length, does a gemini surfactant or its monomeric analogue generally show a MORE favorable (more negative) standard Gibbs free energy of micellization, deltaG_mic?",
    "gemini", difficulty="medium",
    source="Consistent with gemini surfactants' lower CMC -- a lower CMC corresponds to a more negative (more favorable) deltaG_mic via the standard deltaG_mic = RT ln(X_cmc) relation.")

if __name__ == "__main__":
    print(f"Category J: {len(QUESTIONS)} questions")
    from collections import Counter
    print(Counter(q.subcategory for q in QUESTIONS))
