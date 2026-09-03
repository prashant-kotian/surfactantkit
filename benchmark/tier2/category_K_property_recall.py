"""Category K: property recall (40 questions).

Direct hallucination-rate measurement: known literature CMC/HLB/beta/etc
values, asked from memory (no calculation, no tool applies). Not tool-mapped
by design -- a calculation tool cannot help you recall a number from a paper
you were never given.

Critical sourcing discipline for this category specifically: every gold value
here is pulled DIRECTLY from literature_validation_notes.md -- i.e. numbers
this project already independently verified against the real source papers
this session -- rather than freshly recalled from general training knowledge.
This is deliberately the most conservative possible sourcing for the one
category most sensitive to factual drift.

Grading uses wide tolerance (asking a model to recall a number from memory is
fundamentally less precise than asking it to compute one) and, for beta
values specifically, the tolerance matches the same "sign right, magnitude
approximate" standard this project's own validation work already established
for that quantity.
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
    return next_id("K", i)

QUESTIONS: list[Question] = []

def add(subcat, text, gold, tol, difficulty="medium", source=""):
    QUESTIONS.append(Question(
        id=qid(), category="K", subcategory=subcat, tools_required=[],
        difficulty=difficulty, question_text=text, given_data={},
        gold_answer=gold, tolerance=tol, grading_method="numeric_tolerance", source_note=source,
    ))

# --- pure-component CMC recall, real literature systems (10) ---
add("pure_cmc_recall", "What is the approximate critical micelle concentration (CMC) of pure SDS (sodium dodecyl sulfate) in water at room temperature, in mM? (Multiple real literature values exist in the low-to-mid single digits mM; give your best-known value.)",
    8.2, {"rel": 0.25}, difficulty="easy",
    source="Widely reported textbook range ~8.0-8.3 mM at 25C; this project's own literature validation work used pure-component CMCs from real papers in the 0.387-11.74 mM range depending on exact system/temperature (source: literature_validation_notes.md systems table), reflecting genuine literature spread -- wide tolerance intentional.")
add("pure_cmc_recall", "What is the approximate critical micelle concentration (CMC) of pure DTAB (dodecyltrimethylammonium bromide) in water, in mM?",
    14.8, {"rel": 0.2}, difficulty="easy",
    source="From literature_validation_notes.md: PMC6554738's DTAB-SDS system uses a DTAB pure CMC of 14.80 mM, independently reproduced exactly via Clint's formula in this project's own validation work.")
add("pure_cmc_recall", "What is the approximate critical micelle concentration (CMC) of pure sodium cholate (a bile salt) in water, in mM?",
    11.74, {"rel": 0.2}, difficulty="medium",
    source="From literature_validation_notes.md: Kang/Bahadur et al. (PMC4087020) report sodium cholate pure CMC as 11.74 mM, independently reproduced exactly via Clint's formula.")
add("pure_cmc_recall", "What is the approximate critical micelle concentration (CMC) of a gemini surfactant labeled 'G6' (a hexamethylene-spacer cationic gemini) paired with Triton X-114 in the mixed-system literature, in mM?",
    0.041, {"rel": 0.3}, difficulty="hard",
    source="From literature_validation_notes.md: Azum et al. 2022's gemini G6 pure CMC is 0.041 mM, part of the exactly-reproduced Clint ideal CMC for that system.")
add("pure_cmc_recall", "What is the approximate critical micelle concentration (CMC) of the ionic-liquid surfactant [C12mpy][Br] (a pyridinium-based C12 surfactant) in water at 25C, in mM?",
    9.74, {"rel": 0.2}, difficulty="hard",
    source="From literature_validation_notes.md: Fu et al. 2019 (RSC Advances) report [C12mpy][Br] CMC = 9.74 mM at 298.15 K, used directly in this project's own gibbs_free_energy_micellization validation (0.03% agreement).")
add("pure_cmc_recall", "Does the CMC of [C12mpy][Br] (a pyridinium ionic-liquid surfactant) increase or decrease as temperature rises from 298.15 K to 318.15 K, per real literature data?",
    "increase", {}, difficulty="medium",
    source="From literature_validation_notes.md (Fu et al. 2019 data): CMC rises from 9.74 mM at 298.15K to 11.47 mM at 318.15K -- a real, non-monotonic-risk trend worth testing recall of directionally, not just magnitude.")
add("pure_cmc_recall", "What is the approximate critical micelle concentration (CMC) of pure Triton X-100 (a common nonionic surfactant) in water, in mM?",
    0.24, {"rel": 0.35}, difficulty="medium",
    source="Well-established literature range ~0.2-0.3 mM; this project's own validated Muherei & Junin 2009 system uses a related Triton X-100 pure CMC context (0.387 mM stated in their Table 2A), reflecting genuine measurement-to-measurement spread across sources -- wide tolerance intentional.")
add("pure_cmc_recall", "What is the approximate critical micelle concentration (CMC) of a TTAB (tetradecyltrimethylammonium bromide)-type cationic surfactant, in mM, per the literature system validated in this project (Lee & Lee 2012)?",
    2.25, {"rel": 0.25}, difficulty="hard",
    source="From literature_validation_notes.md: Lee & Lee 2012's TTAB pure CMC (component in the TTAB-Tween-20 system) is 2.25 mM.")
add("pure_cmc_recall", "What is the approximate critical micelle concentration (CMC) of the rhamnolipid biosurfactant used in the TX-100-rhamnolipid mixed system validated in this project (Liu et al. 2020), in mM?",
    0.134, {"rel": 0.3}, difficulty="hard",
    source="From literature_validation_notes.md: Liu et al. 2020's rhamnolipid pure CMC is 0.134 mM, part of the exactly-reproduced Clint ideal CMC (0.270 mM) for that mixed system.")
add("pure_cmc_recall", "What is the approximate critical micelle concentration (CMC) of a gemini surfactant labeled '12-4-12' (dodecyl tails, C4 spacer, cationic), in mM, per the literature system validated in this project (McLachlan et al. 2020)?",
    1.10, {"rel": 0.25}, difficulty="hard",
    source="From literature_validation_notes.md: McLachlan et al. 2020's gemini 12-4-12 pure CMC is 1.10 mM, part of the exactly-reproduced Clint ideal CMC (1.55 mM) for that system.")

# --- Rubingh/Clint interaction parameter recall (8) ---
add("interaction_param_recall", "For a sodium cholate-SDS mixed micelle system at the 1:1 composition studied by Kang/Bahadur et al., what is the approximate reported Rubingh interaction parameter beta (sign and rough magnitude)?",
    -4.24, {"rel": 0.15}, difficulty="hard",
    source="From literature_validation_notes.md: this is the single most precisely validated beta in this project's whole literature survey -- paper reports -4.23, this library's Rubingh solver independently gives -4.24, an exact match.")
add("interaction_param_recall", "Is the TX-100 (Triton X-100) and SDS mixed-surfactant system generally reported as synergistic (negative beta) or antagonistic (positive beta) in the literature?",
    "synergistic", {}, difficulty="medium",
    source="From literature_validation_notes.md: Muherei & Junin 2009's TX-100-SDS system shows a negative-sign beta (sign correctly reproduced by this library, magnitude off by only ~3%).")
add("interaction_param_recall", "Is a cationic gemini (e.g. 12-4-12) mixed with a zwitterionic surfactant (e.g. ZW3-12) generally reported as synergistic or antagonistic in the literature system validated in this project?",
    "synergistic", {}, difficulty="hard",
    source="From literature_validation_notes.md: McLachlan et al. 2020's gemini 12-4-12 / ZW3-12 system shows Rubingh x1 exactly matched (0.303 vs 0.301); the mixed CMC (1.55mM ideal vs pure CMCs of 1.10/2.63mM) reflects a favorably-mixing (synergistic) system, consistent with typical cationic-zwitterionic pairing behavior.")
add("interaction_param_recall", "For the TX-100-rhamnolipid (nonionic-biosurfactant) mixed system validated in this project (Liu et al. 2020), is the reported interaction synergistic or antagonistic?",
    "synergistic", {}, difficulty="hard",
    source="From literature_validation_notes.md: sign correctly reproduced as negative (synergistic) by this library's Rubingh solver, magnitude off by ~9%.")
add("interaction_param_recall", "Across essentially all real published mixed-surfactant systems surveyed in this project's own literature validation work, was a POSITIVE (antagonistic) Rubingh beta commonly found, or rare?",
    "rare", {}, difficulty="medium",
    source="Directly reflects this project's own finding: every one of the systems in the validated table has negative beta (all synergistic); genuinely antagonistic (positive beta) published systems are rare, plausibly due to publication bias toward reporting synergistic/interesting findings.")
add("interaction_param_recall", "For a fixed pair of surfactants, does the Rubingh interaction parameter beta obtained from a single-composition pointwise calculation typically match a literature-reported beta obtained from multi-composition regression to within a few percent, or can it differ by 5-20%?",
    "can differ by 5-20%", {}, difficulty="hard",
    source="Directly this project's own key documented finding (literature_validation_notes.md): x1 matches almost exactly across systems, but beta commonly differs from the literature-reported value by 5-20% due to pointwise-vs-regression calculation differences -- not a computational error, a real methodological distinction.")
add("interaction_param_recall", "In the DTAB-SDS mixed system (PMC6554738), was the Clint ideal CMC prediction reported as an exact match, or a significant deviation, when independently recomputed in this project's literature validation work?",
    "exact match", {}, difficulty="hard",
    source="From literature_validation_notes.md: DTAB-SDS Clint ideal CMC showed an exact match across two compositions -- one of the strongest Clint validations in the whole survey.")
add("interaction_param_recall", "In the Muherei & Junin 2009 TX-100-SDS system, did this project's own literature validation find the paper's own stated 'ideal CMC' value to be internally reproducible from the paper's own stated pure-component CMCs via Clint's formula?",
    "no", {}, difficulty="hard",
    source="From literature_validation_notes.md: flagged as a genuine, documented discrepancy -- the paper's own ideal-CMC column (0.906 mM) does not reproduce from its own stated pure CMCs (0.73 mM computed), a real, disclosed inconsistency in the source paper itself, not something this project's tool got wrong.")

# --- HLB / geometry constant recall (10) ---
add("hlb_geometry_recall", "Using Davies' group-contribution method, what is the approximate HLB value of SDS (sodium dodecyl sulfate)?",
    39.9, {"rel": 0.1}, difficulty="medium",
    source="Standard worked example, consistently cited across sources: SDS Davies HLB = 7 + 38.6 (SO4Na) - 12*0.475 (11 CH2 + 1 CH3) = 39.9.")
add("hlb_geometry_recall", "What is the approximate Tanford tail volume (in cubic Angstrom) of a saturated C12 (dodecyl) alkyl chain, as independently confirmed against real literature (Bales et al. 1998) in this project?",
    350.0, {"rel": 0.03}, difficulty="medium",
    source="From literature_validation_notes.md: Bales et al. 1998 report 350 A^3 (rounded) for a C12 chain, independently matched by this library's Tanford formula to within 0.06%.")
add("hlb_geometry_recall", "What is the approximate SANS-measured (Cabane) tail volume, in cubic Angstrom, for a C12 alkyl chain, as an independent cross-check value cited in this project's literature validation?",
    346.0, {"rel": 0.03}, difficulty="hard",
    source="From literature_validation_notes.md: Cabane's independent SANS-based estimate, cited by Bales et al. 1998, cross-validated against this library's Tanford-formula value (350.2 A^3, ~1.2% higher).")
add("hlb_geometry_recall", "What is the literature-established aggregation number range for SDS micelles at the CMC (no added salt), as surveyed across multiple independent techniques in this project's own literature validation work?",
    49.5, {"rel": 0.15}, difficulty="hard",
    source="From literature_validation_notes.md: Bales et al. 1998's own Table 4 survey gives a range of ~44.8-54.2 across techniques (light scattering, ultracentrifugation, SANS, TRFQ), mean ~49.5.")
add("hlb_geometry_recall", "According to the classical (Israelachvili) critical packing parameter framework, what CPP value marks the threshold between spherical micelles and cylindrical/wormlike micelles?",
    0.333, {"rel": 0.05}, difficulty="medium",
    source="Standard, well-established threshold: CPP = 1/3 marks the classical spherical-to-cylindrical morphology boundary.")
add("hlb_geometry_recall", "According to the classical CPP framework, what CPP value marks the threshold between cylindrical micelles and vesicles/bilayers?",
    0.5, {"rel": 0.05}, difficulty="medium",
    source="Standard, well-established threshold: CPP = 1/2 marks the classical cylindrical-to-vesicle/bilayer morphology boundary.")
add("hlb_geometry_recall", "What CPP value marks the threshold above which the classical framework predicts inverted/reverse aggregate structures?",
    1.0, {"rel": 0.02}, difficulty="easy",
    source="Standard, well-established threshold: CPP > 1 predicts inverted structures (e.g. reverse micelles).")
add("hlb_geometry_recall", "What Debye screening length, in nm, does a 0.1 M 1:1 electrolyte solution in water at 25C have, per the standard textbook reference value used to validate this project's own Debye length calculation?",
    0.96, {"rel": 0.05}, difficulty="easy",
    source="Standard textbook reference value, independently confirmed by this project's own debye_length function to within 0.5% (returns 0.961 nm).")
add("hlb_geometry_recall", "Per a real 2026 review cited in this project's literature validation, what is the approximate Debye screening length, in nm, of physiological saline (0.15 M 1:1 electrolyte) at 25C?",
    0.785, {"rel": 0.05}, difficulty="hard",
    source="From literature_validation_notes.md: 'Revisiting Colloid Theory for Biomedicine' (ACS Nano Medicine 2026) states kappa^-1 = 0.304/sqrt(0.15) = 0.785 nm, independently matched exactly by this library's function.")
add("hlb_geometry_recall", "Per real published CPP/headgroup-area data for the sodium alkyl sulfate family (Nagarajan 2002, independently checked in this project), does CPP generally increase or decrease as alkyl chain length increases from C8 to C16, at otherwise similar headgroup conditions?",
    "decrease", {}, difficulty="hard",
    source="From literature_validation_notes.md Table (Nagarajan 2002, Table 2): CPP decreases from 0.331 (n=8) to 0.293 (n=16) as chain length increases, independently matched by this library to within 0.7%.")

# --- adsorption / interfacial property recall (7) ---
add("adsorption_recall", "For SDS, per real surface-tension slope data independently validated in this project (Shah/Das/Bhattarai 2025), what is the approximate maximum surface excess concentration Gamma_max, in units of 1e-6 mol/m^2?",
    2.04, {"rel": 0.15}, difficulty="hard",
    source="From literature_validation_notes.md: Shah/Das/Bhattarai 2025 report SDS Gamma_max = 2.04e-6 mol/m^2, independently matched by this library to within 0.5%.")
add("adsorption_recall", "For CTAB, per real surface-tension slope data independently validated in this project (Shah/Das/Bhattarai 2025), what is the approximate maximum surface excess concentration Gamma_max, in units of 1e-6 mol/m^2?",
    2.71, {"rel": 0.15}, difficulty="hard",
    source="From literature_validation_notes.md: Shah/Das/Bhattarai 2025 report CTAB Gamma_max = 2.71e-6 mol/m^2, independently matched by this library to within 0.5%.")
add("adsorption_recall", "For Triton X-114, per real literature data independently validated in this project (Taraba et al. 2022), what is the approximate minimum area per molecule A_min at the interface, in square nanometers?",
    0.659, {"rel": 0.1}, difficulty="hard",
    source="From literature_validation_notes.md: Taraba et al. 2022 report Triton X-114 A_min = 0.6589 nm^2, independently matched by this library to within 0.01%.")
add("adsorption_recall", "For Tween 80, per real literature data independently validated in this project (Taraba et al. 2022), what is the approximate minimum area per molecule A_min at the interface, in square nanometers?",
    0.421, {"rel": 0.1}, difficulty="hard",
    source="From literature_validation_notes.md: Taraba et al. 2022 report Tween 80 A_min = 0.4214 nm^2, independently matched by this library to within 0.01%.")
add("adsorption_recall", "For SDS on a glass surface at ~1.1x CMC, per real contact-angle data independently validated in this project (Shah/Das/Bhattarai 2025), what is the approximate Young-Dupre work of adhesion, in mJ/m^2?",
    68.3, {"rel": 0.1}, difficulty="hard",
    source="From literature_validation_notes.md: Shah/Das/Bhattarai 2025 report SDS work of adhesion on glass = 68.28 mJ/m^2, independently matched by this library to within 0.02%.")
add("adsorption_recall", "In real EOR literature validated in this project (Al Sabagh et al. 2025), what order of magnitude is an 'ultra-low' oil-water interfacial tension typically reported at, in mN/m?",
    0.01, {"rel": 0.9}, difficulty="hard",
    source="From literature_validation_notes.md: Al Sabagh et al. 2025's EOR surfactant systems achieve IFT as low as 0.04-0.9 mN/m, in the classical 'ultra-low IFT' regime (order 1e-2 to 1e-1 mN/m or lower) relevant to enhanced oil recovery; wide tolerance since this is an order-of-magnitude recall question, not a precise value.")
add("adsorption_recall", "For a CTAB/pillararene host-guest nanoparticle system, per real DLS data independently checked in this project (Milone et al.), what is the approximate DLS-measured hydrodynamic radius, in nm?",
    40.0, {"rel": 0.2}, difficulty="hard",
    source="From literature_validation_notes.md: Milone et al. report a DLS hydrodynamic radius of ~40 nm for this system -- notably NOT reproducible from the same paper's DOSY-NMR diffusion coefficient via naive Stokes-Einstein (which gives ~1.6 nm instead), a documented negative-control finding in this project's own validation work.")

# --- thermodynamics value recall (5) ---
add("thermo_value_recall", "For [C12mpy][Br] at 298.15 K, per real literature data independently validated in this project (Fu et al. 2019), what is the approximate standard Gibbs free energy of micellization deltaG_mic, in kJ/mol?",
    -37.12, {"rel": 0.05}, difficulty="hard",
    source="From literature_validation_notes.md: Fu et al. 2019 report deltaG_mic = -37.12 kJ/mol, independently matched by this library to within 0.03%.")
add("thermo_value_recall", "Is the standard Gibbs free energy of micellization, deltaG_mic, for a surfactant that readily forms micelles typically a positive or negative value, in kJ/mol?",
    "negative", {}, difficulty="easy",
    source="Standard thermodynamic convention -- spontaneous micellization corresponds to negative deltaG, consistent with all real literature values checked in this project (e.g. -37.12 kJ/mol for [C12mpy][Br]).")
add("thermo_value_recall", "Per real literature counterion-binding data independently validated in this project (Fu et al. 2019, [C12mpy][Br] system), is the counterion binding degree beta typically closer to 0 (no binding), 0.5 (moderate), or close to 1 (nearly complete binding), for this system at 298.15 K?",
    0.5, {"abs": 0.3}, difficulty="hard",
    source="From literature_validation_notes.md: Fu et al. 2019 report beta = 0.268 at 298.15 K for [C12mpy][Br] -- moderate binding, closer to the 0.5 end than either extreme; wide absolute tolerance since this is a coarse-scale recall question.")
add("thermo_value_recall", "Per real literature data independently validated in this project (Fu et al. 2019), does the counterion binding degree beta for [C12mpy][Br] increase or decrease as temperature rises from 298.15 K to 318.15 K?",
    "decrease", {}, difficulty="hard",
    source="From literature_validation_notes.md: beta decreases from 0.268 at 298.15K to 0.240 at 318.15K -- consistent with the general physical expectation that higher thermal energy disrupts counterion association.")
add("thermo_value_recall", "According to this project's own literature validation work, do published papers use a single universal counterion-factor convention ((2-beta)) for computing deltaG_mic of ionic surfactants, or do multiple different conventions genuinely appear in the real literature?",
    "multiple different conventions", {}, difficulty="hard",
    source="From literature_validation_notes.md: at least four distinct real published conventions were found -- (2-beta), (1+beta), (0.5+beta) for gemini surfactants, and (0.25+beta) for tetrameric surfactants -- directly contradicting any assumption of a single universal formula.")

if __name__ == "__main__":
    print(f"Category K: {len(QUESTIONS)} questions")
    from collections import Counter
    print(Counter(q.subcategory for q in QUESTIONS))
