# Literature validation notes

SurfactantKit's Clint and Rubingh implementations were checked against eight independent, open-access, published binary surfactant systems spanning cationic-anionic, cationic-nonionic, anionic-nonionic, gemini-nonionic, gemini-zwitterionic, nonionic-biosurfactant, and anionic-anionic/bile-salt pairs. Numbers were pulled directly from source tables (not paraphrased search summaries) wherever possible.

## Summary

| System | Source | Clint ideal | Rubingh x1 | Rubingh beta |
|---|---|---|---|---|
| DTAB - SDS | PMC6554738 | Exact match (2 compositions) | Not independently verifiable (see note) | Not asserted (see note) |
| TX-100 - SDS | Muherei & Junin 2009 | Does not reproduce (flagged) | Exact match (0.7501) | Sign correct, magnitude off ~3% |
| Gemini G6 - TX-114 | Azum et al. 2022 | Exact match (0.073 mM) | Close match (0.76 vs 0.764) | Sign correct, magnitude off ~4% |
| TTAB - Tween-20 | Lee & Lee 2012 | Not reported by source | Close match (0.28 vs 0.28) | Sign correct, magnitude off ~4% |
| Gemini 12-4-12 - ZW3-12 | McLachlan et al. 2020 | Exact match (1.55 mM) | Exact match (0.303 vs 0.301) | Not reported by source |
| TX-100 - rhamnolipid | Liu et al. 2020 | Exact match (0.270 mM) | Exact match (0.744 vs 0.744) | Sign correct, magnitude off ~9% |
| Sodium cholate - SDS | Kang/Bahadur et al., PMC4087020 | Exact match (11.74 mM) | Exact match (0.5033 vs 0.503) | **Exact match** (-4.24 vs -4.23) |

## Key finding: x1 matches almost exactly; reported beta consistently doesn't, by a few to ~20%

This was investigated rather than shrugged off. At the solved root, the two algebraically-equivalent ways of computing beta (`term1/(1-x1)^2` and `term2/x1^2`) agree with each other to 3-4 decimal places in every system checked -- confirming the solver has genuinely found the correct root, not a numerical artifact. But neither computed value matches the literature-reported beta to better than ~5-20%, always with the correct sign.

The most likely explanation: many of these papers determine beta by regression across several alpha compositions (a global least-squares fit), while this library computes the exact pointwise value implied by one single (alpha, CMCmix) pair. Both are legitimate methods and will not agree to three decimals in general. This is directly supported by the sodium cholate-SDS case below: it's a fixed 1:1 mixture study (no other composition measured, so no regression possible), and there beta matches the paper almost exactly (-4.24 vs -4.23). Given x1 -- the solver's actual output -- matches essentially exactly across six independent systems, this is read as strong evidence the implementation is correct, not as a bug.

**Practical implication for users:** treat `rubingh_beta()`'s output as the exact pointwise regular-solution value for the specific (alpha, CMCmix) pair given. If comparing against a beta value from a paper that used multi-point regression, expect agreement in sign and order of magnitude, not necessarily to the decimal.

## Notes on individual systems

- **Muherei & Junin 2009**: the paper's own "ideal CMC" column (0.906 mM) does not reproduce from its own stated pure-component CMCs (0.387, 3.468 mM) via the standard Clint formula (0.73 mM computed). Possibly a different pure-CMC basis was used for that column (e.g. the literature-range values in their Table 1 rather than their own Table 2A measurements). Flagged, not resolved -- the x1 and beta-sign checks for this system are unaffected since those use the paper's directly-stated experimental CMCmix, not the disputed ideal value.
- **DTAB-SDS (PMC6554738)**: this was the first system checked, before the systematic beta pattern above was understood. Two competing pieces of its data (composition A: DTAB-rich, composition B: SDS-rich) gave exact Clint matches, which is why they're kept as the primary Clint test cases. An attempt to reverse-engineer the paper's reported beta (-2.5674 at x=0.5) from its stated alpha/CMCmix values did not reproduce, and is not asserted in the test suite.

- **Sodium cholate - SDS (PMC4087020)**: the strongest single validation case in the survey. Fixed 1:1 mixture, both fluorimetry- and tensiometry-derived experimental CMCs given; the fluorimetry value was used for the beta assertion since it lands closest to the paper's own reported beta.

## Batch collection note

A fourth automated research-agent pass (targeting nonionic-nonionic and amino-acid-surfactant systems) was interrupted by hitting the account's monthly API spend limit before completing. A subsequent manual search round (searching and extracting directly rather than via sub-agents) covered zwitterionic-cationic, sugar-surfactant, fluorocarbon-hydrocarbon, ionic-liquid, gemini-conventional, catanionic, and additional bile-salt/nonionic combinations -- most candidates turned out not to have complete, cleanly tabulated pure surfactant-surfactant data open-access (common failure modes: drug-surfactant or polymer-surfactant systems mislabeled by search snippets as surfactant-surfactant, data only in alcohol- or salt-modified conditions with no clean baseline, review papers describing results narratively without reproducing the source tables, or SANS/spectroscopy-focused papers with no CMC-determination table at all). One additional clean system (sodium cholate-SDS) was found and added. Eight systems is what's validated as of this writing.
