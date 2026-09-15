# Where LLMs genuinely refuse, and where a toolkit can convert that into confidence

Purpose: not to lower LLM credibility, but to find exactly where SurfactantKit
can give an LLM a real path to an answer it currently lacks — versus where the
correct behavior for both the LLM and the tool is to keep refusing, because no
tool in the world can supply information that was never measured. Every claim
below is sourced from a real fetch/search performed 2026-09-15.

## 1. Documented taxonomy of LLM refusal

Real, current research treats abstention as a first-class capability, not a
safety afterthought. The clearest practical split: **data uncertainty**
(the question lacks a unique objective answer because information is
genuinely missing or ambiguous) versus **model uncertainty** (the question
has a unique answer but exceeds the model's own capability). A second, more
granular framework names three concrete triggers: **P-Ambiguity** (linguistic
ambiguity), **P-Contradiction** (logically inconsistent facts), and
**P-MissingInfo** (absence of critical information) — this last one is
exactly the shape of every real refusal point this project has found (a
missing counterion-binding degree, a missing electrolyte condition, a group
number that was never tabulated).
[Epistemic Abstention in AI Systems](https://www.emergentmind.com/topics/epistemic-abstention),
[Know Your Limits: A Survey of Abstention in Large Language Models](https://arxiv.org/pdf/2407.18418),
[Ambiguity Collapse by LLMs: A Taxonomy of Epistemic Risks](https://arxiv.org/abs/2603.05801)

**Abstention does not improve with scale.** AbstentionBench evaluated 20
frontier LLMs across 20 datasets spanning unknown answers, underspecification,
false premises, and outdated information, and found abstention "an unsolved
problem" — bigger/newer models are not reliably better at knowing when to
say "I don't know." This matters directly for this project: don't expect the
refusal-calibration gap to close on its own as models improve; it's a
real, durable axis to keep testing.
[AbstentionBench: Reasoning LLMs Fail on Unanswerable Questions](https://arxiv.org/pdf/2506.09038)

**When a required parameter is missing, models default to silently assuming
a value far more often than they refuse.** Real documented behavior:
"LLMs implicitly complete missing information with a default assumption,"
and in financial reasoning specifically, models "predict answers as correct
even though the conclusion is not logically supported by the observed
premises alone" — i.e. they proceed on an unstated assumption rather than
flagging the gap. This is the opposite failure mode from safety-trained
refusal, and it's the one this project's own "don't guess a system
parameter" trap category is built to catch.
[LLMs as Implicit Imputers: Uncertainty Should Scale with Missing Information](https://arxiv.org/html/2605.13188)

## 2. How tool access changes refusal behavior — a directly relevant, important finding

**Merely making an irrelevant tool available can make models refuse MORE,
even when they never call it.** A real, controlled study found six pooled
models' answer rate on closed-domain questions dropped from 98.2% to 63.5%
once a related-but-unnecessary tool was made available — with Gemini 2.5
Flash-Lite falling from 99.4% to 23.4% while calling the tool in only 7.8% of
those trials. The mechanism: models misread a narrow tool's presence as a
**scope restriction on what they're allowed to answer**, not as an optional
extra capability — one model, when a weather tool was present, said it
"can only provide weather information" even for a question it answered
correctly with no tools at all. The fix that recovered up to 45.6 points
of lost answers was a single explicit instruction: tools are optional,
answer directly from your own knowledge when they don't apply.
[When Tools Get in the Way: The Effect of Unnecessary Tool Availability on LLM Answering](https://arxiv.org/html/2609.14157)

**This is directly, urgently relevant to this project's own stated design
principle** ("our tool must enhance the ability of LLMs and not add
unnecessary MCP which AI doesn't need at all"). The risk isn't hypothetical:
attaching SurfactantKit's MCP tools broadly, without a clear signal about
scope, could make a model MORE conservative/refusal-prone on questions it
would have answered fine on its own — the opposite of the intended effect.
The same study's utility-control condition shows a properly-scoped, genuinely
relevant tool gets used 85.8-100% of the time when actually needed — so the
fix isn't "expose fewer tools," it's "make each tool's own scope and
optionality explicit," which SurfactantKit's MCP tool docstrings already do
reasonably well (each one states its own real applicability limits), but
this is worth auditing explicitly against this finding before any re-module.

## 3. This project's own four real refusal points, classified

**(a) = genuinely unanswerable, no tool could ever resolve it — correct
behavior is refusal, and the toolkit's value is enforcing/confirming that,
not letting either the LLM or the tool guess.**
**(b) = the toolkit already has (or could cheaply add) a real path to an
actual answer — the refusal is currently premature, not correct.**

1. **Missing counterion-binding degree for deltaG_mic** — **(a), pure.**
   This is a real physical quantity (conductivity slope-ratio, EPR, or a
   literature value for THAT specific compound) that does not exist anywhere
   for most real compounds unless someone measured it. No SMILES parsing,
   no curve-fit, no amount of computation invents a measurement that was
   never taken. This is the textbook P-MissingInfo case. The toolkit's real
   value here is exactly what it already does: refuse cleanly and name what's
   missing, rather than let a model's documented "silently assume a default"
   tendency (section 1 above) substitute a plausible-looking beta.

2. **Missing Davies HLB group numbers for gemini/glycolipid headgroups** —
   **currently (a), but convertible to (b) with real, disclosed work, not
   guessing.** This is NOT a case where the information doesn't exist — it's
   that Davies' 1957 table was never extended to these structural classes.
   SurfactantKit already has the infrastructure to do this properly:
   `derive_davies_group_number_from_griffin` (hlb.py) uses Davies' own
   cross-calibration methodology to back-solve a real new group number for a
   functional group his original table never assigned one to — already used
   for several other previously-missing groups (sulfonate, sultaine,
   phosphate-diNa, amphodiacetate). Doing the same real derivation for a
   glycoside/gemini-quaternary headgroup would genuinely convert this refusal
   into an answer — this is the single highest-value, most concrete "give the
   LLM more real confidence" opportunity found in this pass, not a new tool
   category, just applying an existing, already-validated method to two more
   real structural classes.

3. **HLD's cationic-quaternary-only HLB branch, misapplied to a non-quat
   cationic surfactant** — **(a), pure.** The empirical k/Cc bridge was
   calibrated specifically for quaternary-ammonium surfactants; using it for
   a different cationic headgroup is a real methodological mismatch, not a
   data gap. No tool should ever produce a number here — the toolkit's value
   is exactly the applicability check that already exists, and a good
   benchmark question is one that tests whether an LLM recognizes this limit
   unprompted, matching how GZ-03/GZ-04's HLB-Davies refusals were designed.

4. **Missing electrolyte condition for the Gibbs prefactor (n=1 vs n=2)** —
   **splits into two different sub-cases, a genuinely interesting nuance.**
   For the Gibbs-prefactor-DEPENDENT path itself (Gamma_max/deltaG_mic via
   the pseudo-phase route) — **(a), pure**: whether excess electrolyte was
   added to the beaker is a real fact about the solution that no amount of
   curve-fitting or SMILES parsing can recover; it must be told or refused.
   But for the geometry/packing-parameter path specifically — **(b), already
   solved and unused**: Nagarajan's method (already in `cascade_
   cmc_surface_tension.py`) doesn't need an externally-stated electrolyte
   condition at all; it derives the relevant screening length directly from
   the CMC itself via Debye-Huckel. GZ-10 refused the corrected-CPP question
   for a reason that only applies to a DIFFERENT downstream quantity than the
   one actually needed — a real, concrete instance of a refusal trigger
   being applied one step too broadly. This is worth designing future
   questions around precisely: the same missing input can legitimately gate
   one sub-question while leaving a sibling sub-question fully answerable.

## 4. Synthesis — where to invest, where to reinforce refusal

**Highest-value real investment**: extend `derive_davies_group_number_from_
griffin`'s already-validated cross-calibration methodology to gemini and
glycolipid headgroups (item 2). This converts a currently-premature refusal
into a genuine answer, using a method this project has already proven works,
not a new capability class.

**Second-highest-value investment, no new code**: build questions that
separate the Nagarajan-geometry path from the Gibbs-prefactor path on the
SAME real dataset (item 4) — this doesn't require building anything, just
asking the right combination of sub-questions, and it directly tests whether
an LLM (or a careless refusal) conflates two genuinely different missing-input
situations that happen to share a surface-level "electrolyte condition" label.

**Reinforce, don't fix**: items 1 and 3 are genuinely correct refusals. The
toolkit's real contribution there isn't new capability, it's the toolkit
holding the line a model's own documented guessing tendency (section 1) would
otherwise cross — this is precisely the "MCP gives LLMs the confidence to say
no correctly" framing the user asked for, and it should stay a refusal in
every future question built around these two cases.

**Structural risk to audit before any MCP re-module**: section 2's finding
that unnecessary tool availability can suppress correct answers by tens of
percentage points, even unused, is a real, direct threat to this project's
own explicit design principle. Before adding any new MCP tool, confirm its
docstring states its own scope/optionality as explicitly as the source study's
fix requires — this is a real, concrete design check, not a vague caution.
