"""Category M: conceptual theory (40 questions).

"Why does X happen" / mechanism / interpretation questions -- Clint, Rubingh,
Gibbs, CPP theory, not calculation. Not tool-mapped: a calculation tool cannot
help explain WHY a formula behaves a certain way, only compute its output.
Grounded in the same theory already implemented and literature-validated in
SurfactantKit itself (mixed_micelle.py, adsorption.py, cpp.py docstrings and
literature_validation_notes.md), so every fact here traces to something
already verified elsewhere in this project, not fresh recall.
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
    return next_id("M", i)

QUESTIONS: list[Question] = []

def add(subcat, text, gold, difficulty="medium", source=""):
    QUESTIONS.append(Question(
        id=qid(), category="M", subcategory=subcat, tools_required=[],
        difficulty=difficulty, question_text=text, given_data={},
        gold_answer=gold, grading_method="category_match", source_note=source,
    ))

# --- Rubingh / Clint interpretation (12) ---
add("rubingh_interpretation", "In Rubingh's regular-solution theory for mixed micelles, what does a NEGATIVE interaction parameter beta physically indicate about the two surfactants' interaction in the mixed micelle, compared to ideal (Clint) mixing?",
    "synergistic", source="beta < 0 means the mixed system's CMC is lower than ideal mixing predicts -- net attractive interaction between the two surfactants in the mixed micelle, i.e. synergism.")
add("rubingh_interpretation", "In Rubingh's regular-solution theory, what does a POSITIVE interaction parameter beta physically indicate?",
    "antagonistic", source="beta > 0 means the mixture behaves worse than ideal mixing -- net repulsive/unfavorable interaction, i.e. antagonism (rare in practice, per this project's own literature survey -- most published systems show negative beta, since antagonistic systems are less commonly published).")
add("rubingh_interpretation", "What does a Rubingh interaction parameter beta close to ZERO indicate about the mixed system?",
    "near-ideal mixing", source="beta near 0 means the mixture behaves close to Clint's ideal-mixing prediction -- negligible net interaction beyond simple dilution.")
add("rubingh_interpretation", "Why does Clint's ideal mixing equation alone fail to predict synergistic CMC lowering in real ionic-nonionic surfactant mixtures?",
    "assumes no interaction between components", source="Clint's equation assumes the two surfactants behave as if independent (no interaction term) -- it cannot capture the electrostatic/steric interactions that Rubingh's beta parameter explicitly models.")
add("rubingh_interpretation", "Physically, why do cationic-anionic surfactant mixtures typically show STRONGLY synergistic (very negative beta) mixed-CMC behavior, more so than cationic-nonionic mixtures?",
    "electrostatic attraction between oppositely charged headgroups", source="Opposite headgroup charges attract each other directly, providing a strong favorable interaction beyond what steric/dilution effects alone (as in ionic-nonionic mixing) can provide.")
add("rubingh_interpretation", "In Rubingh's theory, what does the micellar mole fraction x1 represent, as distinct from the bulk mole fraction alpha1?",
    "micelle composition", source="alpha1 is the overall/bulk composition of surfactant 1 in solution; x1 is specifically the composition WITHIN the micelle, which generally differs from alpha1 due to preferential partitioning of one component into micelles.")
add("rubingh_interpretation", "Rosen's mixed-monolayer theory uses the same regular-solution mathematics as Rubingh's mixed-micelle theory, but applied to what different type of experimental data?",
    "surface tension data", source="Rosen's beta^sigma describes interaction in the adsorbed surface monolayer, using surface-tension-derived concentrations, not CMC data -- same math, different physical system, as documented in this project's own mixed_micelle.py.")
add("rubingh_interpretation", "Why might a mixed surfactant system show a NEGATIVE beta in the mixed MICELLE (Rubingh) but a different-magnitude (though often still negative) beta^sigma in the adsorbed MONOLAYER (Rosen), for the same two surfactants?",
    "different packing environments", difficulty="hard",
    source="The bulk micelle core and the 2D air-water interface impose different packing/curvature constraints on the same two surfactants, so the interaction parameter -- which reflects how favorably the two headgroups/tails pack together -- need not be identical in magnitude between the two environments, even though both are captured by the same regular-solution mathematical form.")
add("rubingh_interpretation", "The Corrin-Harkins relation describes CMC as a log-linear function of what variable, for ionic surfactants?",
    "counterion (added salt) concentration", source="log(CMC) decreases linearly with log(counterion concentration) for ionic surfactants -- electrostatic screening by added salt lowers headgroup repulsion, lowering CMC.")
add("rubingh_interpretation", "Physically, why does adding salt (e.g. NaCl) to an IONIC surfactant solution generally lower its CMC, while having comparatively little effect on a NONIONIC surfactant's CMC?",
    "screens electrostatic repulsion", source="Added electrolyte ions screen the repulsion between charged (ionic) headgroups via the Debye screening mechanism, easing micellization; nonionic headgroups have no charge to screen, so the mechanism doesn't apply.")
add("rubingh_interpretation", "Rubingh's regular-solution theory assumes the two surfactants mix non-ideally but that the excess entropy of mixing is what, relative to the excess enthalpy?",
    "zero (regular solution assumption)", difficulty="hard",
    source="'Regular solution' theory specifically assumes ideal entropy of mixing but non-ideal (non-zero) enthalpy of mixing -- the interaction parameter beta captures only the enthalpic non-ideality, which is a stated theoretical assumption/limitation, not a universal physical fact.")
add("rubingh_interpretation", "If a Rubingh calculation returns no valid root for x1 in (0,1) given a measured mixed CMC, what does that most likely indicate about the input data (assuming no other numerical error)?",
    "not physically consistent with the theory", difficulty="hard",
    source="A missing valid root is a genuine, real, occasionally-occurring outcome of the theory being applied outside its physically consistent parameter range for that specific data point -- not necessarily a bug -- as directly documented and tested in this project's own rubingh_solve no-solution trap questions.")

# --- Gibbs adsorption / Szyszkowski interpretation (10) ---
add("gibbs_interpretation", "In the Gibbs adsorption isotherm, why does the surfactant-related prefactor n=2 (rather than n=1) apply to a simple 1:1 ionic surfactant with no added salt?",
    "counterion also adsorbs", source="For an ionic surfactant with no swamping electrolyte, both the surfactant ion and its dissociated counterion contribute to the interfacial excess, requiring the factor of 2 in the Gibbs equation; adding excess inert electrolyte (or using a nonionic surfactant) removes this requirement, giving n=1.")
add("gibbs_interpretation", "Physically, what does a STEEPER (more negative) slope of surface tension vs ln(concentration), below the CMC, indicate about the surfactant's adsorption behavior?",
    "higher/more efficient adsorption at the interface (higher Gamma_max)", source="A steeper pre-CMC slope directly corresponds to a larger maximum surface excess Gamma_max via the Gibbs equation -- the surfactant packs more densely/adsorbs more strongly at the interface.")
add("gibbs_interpretation", "Why is the Gibbs adsorption isotherm's pre-CMC slope measured specifically BELOW the CMC, rather than at or above it?",
    "goes into micelles instead", source="Once micelles form, added surfactant partitions predominantly into micelles rather than continuing to adsorb at the interface, so surface tension becomes essentially constant above the CMC -- the pre-CMC region is where the slope actually reflects ongoing interfacial adsorption.")
add("gibbs_interpretation", "What does the minimum area per molecule, A_min, computed from Gamma_max, physically represent?",
    "area per molecule at saturation", source="A_min = 1/(N_A * Gamma_max) is the reciprocal-scaled inverse of the maximum packing density -- smaller A_min means denser interfacial packing.")
add("gibbs_interpretation", "In the Szyszkowski equation, what does the constant K physically represent?",
    "adsorption affinity constant", source="Larger K means the surfactant reaches significant surface-tension reduction at lower bulk concentration -- a measure of adsorption affinity, analogous to a Langmuir equilibrium constant.")
add("gibbs_interpretation", "Why is the Szyszkowski equation only considered valid/applicable BELOW the CMC?",
    "models monomer adsorption", source="The equation relates surface tension to bulk MONOMER concentration; above the CMC, added surfactant goes into micelles rather than increasing free monomer concentration, so the equation's underlying assumption breaks down.")
add("gibbs_interpretation", "Physically, why does a gemini surfactant typically show a LOWER equilibrium surface tension at the CMC (pi_CMC higher / more surface-tension reduction) than its monomeric analogue?",
    "denser interfacial packing", source="Consistent with gemini surfactants' generally lower CMC and stronger interfacial packing -- two covalently linked tails occupy the interface more efficiently per adsorbed molecule than two independent monomeric molecules would.")
add("gibbs_interpretation", "What is the key physical/experimental difference between the Gibbs adsorption slope method and the Szyszkowski equation, given the same surface-tension-vs-concentration dataset?",
    "slope versus full curve fit", difficulty="hard",
    source="Gibbs adsorption is a local-slope method (one number, Gamma_max); Szyszkowski is a full functional fit across the whole isotherm -- related (Szyszkowski's own slope recovers the same Gamma_max, as directly cross-validated in this project's own test suite) but operationally different tools.")
add("gibbs_interpretation", "Why does an ionic surfactant in the presence of excess inert electrolyte (e.g. high background NaCl) behave like n=1 in the Gibbs equation, similar to a nonionic surfactant, even though it's still charged?",
    "excess electrolyte dominates", source="With a large excess of inert electrolyte already present, the surfactant's own counterion's contribution to the interfacial excess becomes negligible relative to the swamping electrolyte, so the simpler n=1 form applies -- a standard, well-documented exception to the n=2 rule for ionic surfactants.")
add("gibbs_interpretation", "What does it mean, physically, for a Szyszkowski-fitted K value to be described as 'not a universal constant'?",
    "system-specific fitted parameter", source="Unlike e.g. Avogadro's number, K is an empirical fitted parameter for one specific system's adsorption isotherm -- reusing a K value from a different surfactant or condition would be a hallucination-style error, as explicitly documented in this project's own szyszkowski_surface_tension docstring.")

# --- CPP / geometry interpretation (10) ---
add("cpp_interpretation", "Physically, what does the critical packing parameter (CPP) represent, as a ratio?",
    "volume divided by area times length", source="CPP = v/(a0*lc) -- a geometric packing ratio predicting which aggregate shape can be assembled from cone-shaped or cylinder-shaped molecular geometry.")
add("cpp_interpretation", "Why does a LARGE headgroup area (bulky headgroup) relative to tail volume tend to favor SPHERICAL micelles rather than bilayers?",
    "cone-shaped molecule", source="Geometrically, cone-shaped molecules (large head, narrow tail) pack efficiently into high-curvature spherical aggregates; cylinder-shaped molecules (head and tail similar cross-section) pack into flat bilayers instead.")
add("cpp_interpretation", "Why does REDUCING headgroup repulsion (e.g. by adding salt to screen an ionic headgroup's charge) tend to favor a transition from spherical micelles toward wormlike micelles or bilayers?",
    "shrinks effective headgroup area", source="Screening electrostatic repulsion between charged headgroups reduces the effective headgroup area a0 in the CPP formula, increasing CPP toward the cylindrical/bilayer regime -- the same underlying mechanism connecting Debye screening/added salt to observed morphology transitions.")
add("cpp_interpretation", "Why do gemini surfactants often show a lower CPP-driven tendency to form wormlike micelles / higher-order structures than a comparable single-tail surfactant with the same total tail volume, WHEN the spacer keeps the two headgroups far apart?",
    "increases effective headgroup area", difficulty="hard",
    source="A spacer that keeps the two headgroups apart (rather than merged) increases the effective a0 term relative to what a single fused double-tail headgroup would give, which can lower CPP relative to naive expectation -- this is why spacer length/rigidity is such an important, non-obvious design variable in gemini surfactant self-assembly.")
add("cpp_interpretation", "What aggregate morphology does the classical CPP framework predict for CPP values greater than 1?",
    "inverted/reverse structures", source="CPP > 1 corresponds to an inverted-cone molecular shape (headgroup smaller than tail cross-section), predicting inverted micelles or other reverse structures rather than any normal (headgroup-out) aggregate.")
add("cpp_interpretation", "Tanford's tail-volume formula (v = 27.4 + 26.9*nc cubic Angstrom) applies specifically to what kind of hydrocarbon chain?",
    "saturated, unbranched alkyl chains", source="The formula is derived for straight-chain, fully saturated hydrocarbons; branched chains or chains with unsaturation (double bonds) would need a different volume estimate, as explicitly noted in this project's own cpp.py docstring.")
add("cpp_interpretation", "Why is the geometric aggregation number (computed from tail volume and core radius) described as an ESTIMATE rather than an exact measured value?",
    "idealized close packing assumption", source="The geometric calculation assumes the micelle core is a perfect sphere fully and uniformly packed with tail volume -- real micelles have some water penetration, polydispersity, and non-ideal packing, so techniques like fluorescence quenching or scattering give the actual measured aggregation number, which the geometric estimate should only be expected to approximate in order of magnitude.")
add("cpp_interpretation", "Why does the CPP framework break down or become less reliable for predicting the morphology of MIXED surfactant systems, compared to single-component systems?",
    "depends on both components' interaction", difficulty="hard",
    source="CPP as classically formulated is a single-molecule geometric argument; in a mixed aggregate, the true controlling factors also include the two components' relative composition and interaction (captured separately by Rubingh's beta, not CPP), so CPP alone is a less complete predictor for mixed systems -- an important limitation for the mixed-system-focused review this project's PhD work is built around.")
add("cpp_interpretation", "What physically distinguishes the 'critical chain length' lc in the CPP formula from the tail's fully extended (all-trans) length?",
    "less than the fully extended length", difficulty="hard",
    source="Tanford's lc formula (1.5 + 1.265*nc) accounts for the chain's practical maximum extension within a micelle core, which is somewhat less than a perfectly rigid, fully extended all-trans conformation would give, reflecting real chain flexibility/disorder.")
add("cpp_interpretation", "Why does increasing tail chain length (nc) generally increase CPP, for a fixed headgroup area?",
    "volume grows faster than length", difficulty="hard",
    source="Both v and lc scale roughly linearly with nc, but because a0 is held fixed (a headgroup property, not a tail property), the net CPP = v/(a0*lc) ratio's dependence on nc is governed by how the volume-to-length ratio itself scales, which for typical saturated alkyl chains yields a net CPP increase with tail length at fixed headgroup area -- part of why longer-tailed surfactants of the same headgroup class tend toward less curved (more bilayer/vesicle-like) morphologies.")

# --- thermodynamics interpretation (8, to reach the taxonomy target of 40) ---
add("thermo_interpretation", "Why is the standard Gibbs free energy of micellization, deltaG_mic, typically NEGATIVE for a surfactant that readily forms micelles?",
    "micellization is spontaneous", source="A negative deltaG indicates the micellization process is thermodynamically favorable/spontaneous under standard conditions -- consistent with surfactants readily self-assembling above their CMC.")
add("thermo_interpretation", "Why does the counterion factor in the deltaG_mic formula use (2-beta) for an ionic surfactant rather than simply 1, as used for a nonionic surfactant?",
    "counterion binding must be accounted for", source="For an ionic surfactant, both the surfactant ion and its partially-bound counterion contribute to the free energy of micelle formation; beta (the binding degree) adjusts how much the counterion's presence modifies the simple nonionic form -- using factor 1 for an ionic surfactant would omit this real physical contribution.")
add("thermo_interpretation", "A positive deltaH_mic (endothermic micellization) combined with a negative deltaG_mic implies the entropy term deltaS_mic must be what sign, for micellization to still be spontaneous overall?",
    "positive", source="From deltaG = deltaH - T*deltaS, if deltaH is positive and deltaG is negative, deltaS must be sufficiently positive (T*deltaS > deltaH) -- classic 'entropy-driven' micellization, often attributed to the hydrophobic effect releasing ordered water molecules around the tail.")
add("thermo_interpretation", "Why can the two-point van't Hoff method for computing deltaH_mic give a substantially different answer than a paper's own polynomial-fit-derived deltaH at a specific temperature?",
    "assumes constant enthalpy over the interval", difficulty="hard",
    source="The two-point method assumes deltaH is constant across the whole T1-T2 range; if the real CMC(T) curve has curvature (deltaH genuinely varies with T), a two-point slope only gives an average over that interval, which can differ substantially -- including in sign -- from a local derivative at one specific temperature. This was directly confirmed empirically in this project's own literature validation work (Fu et al. 2019 data), not just a theoretical caveat.")
add("thermo_interpretation", "Physically, why does the counterion binding degree beta generally DECREASE (less counterion binding) as temperature increases, for a typical ionic surfactant?",
    "increased thermal motion disrupts counterion association", difficulty="hard",
    source="Higher temperature increases the counterions' thermal/kinetic energy, making them less likely to remain closely associated with the micelle surface against thermal disruption -- a general trend, though the magnitude varies by system.")
add("thermo_interpretation", "What does the entropy of micellization, deltaS_mic, physically capture, beyond just 'how much entropy increases'?",
    "net ordering and disordering effects of micelle formation", source="deltaS_mic reflects the balance between the surfactant molecules becoming more ordered (confined into a micelle) and the surrounding water becoming LESS ordered (released from structured hydration shells around the hydrophobic tail) -- the net sign depends on which effect dominates.")
add("thermo_interpretation", "Why is deltaG_mic generally considered a more reliable/robust quantity to compare across surfactants than deltaH_mic derived from a narrow-range van't Hoff fit?",
    "needs only one CMC value", difficulty="hard",
    source="deltaG_mic can be computed from a single, well-measured CMC at one temperature; deltaH_mic requires CMC values at multiple temperatures and is sensitive to the fitting method and temperature range used, making it inherently noisier/more method-dependent, as directly demonstrated by this project's own literature validation findings.")
add("thermo_interpretation", "Why do researchers sometimes report deltaG_mic using different counterion-factor conventions -- e.g. (2-beta) versus (1+alpha) versus (0.5+beta) for a gemini surfactant -- for what is nominally the same physical quantity?",
    "no single universally standardized convention exists", difficulty="hard",
    source="Different research groups/papers use different, non-interchangeable formulations for how counterion binding enters the free-energy expression, particularly for less-standard architectures like gemini surfactants -- directly confirmed by this project's own literature survey finding at least four distinct real published conventions ((2-beta), (1+beta), (0.5+beta), (0.25+beta)) across different surfactant classes.")

if __name__ == "__main__":
    print(f"Category M: {len(QUESTIONS)} questions")
    from collections import Counter
    print(Counter(q.subcategory for q in QUESTIONS))
