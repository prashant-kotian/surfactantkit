# Current LLM capabilities and limitations (research pass, September 2026)

Purpose: ground the next SurfactantKit/Paper 3 redesign in real, current, sourced
evidence about what LLMs are actually good and bad at — not stale assumptions,
and not this project's intuition alone. Every claim below is sourced from a
real search result or direct fetch performed 2026-09-15; nothing is presented
as fact without a link.

## 1. Anthropic's own stated rationale for MCP

Direct from anthropic.com/news/model-context-protocol (fetched directly): the
core problem MCP addresses is that AI models are "constrained by their
isolation from data — trapped behind information silos and legacy systems,"
where "every new data source requires its own custom implementation, making
truly connected systems difficult to scale." Anthropic frames this as an
**M×N integration problem** (M models × N tools = M×N custom integrations);
MCP turns it into M+N by giving every model and every tool one shared
protocol to implement once. [Introducing the Model Context Protocol](https://www.anthropic.com/news/model-context-protocol)

**Important, and directly relevant to the user's framing**: Anthropic's own
stated rationale is about **data access and integration**, not primarily
about "LLMs can't compute." MCP's original pitch is closer to "give the model
a USB-C port to your actual systems" than "give the model a calculator."
This matters for how SurfactantKit should think of its own MCP value
proposition: the strongest, most defensible case is not "LLMs can't do
surfactant math" (increasingly false, see below) but "LLMs have no access to
a validated, literature-cross-checked computational core, and no way to know
when a computation requires information (a counterion-binding degree, an
electrolyte condition) that only exists in a specific paper or a specific
researcher's own unpublished data."

Sources: [Introducing the Model Context Protocol](https://www.anthropic.com/news/model-context-protocol), [Model Context Protocol (MCP) — A Deep Dive, WWT](https://www.wwt.com/blog/model-context-protocol-mcp-a-deep-dive)

## 2. Real, currently-documented LLM weaknesses relevant to this project

**(a) Confident guessing is a trained-in incentive problem, not just a random
failure.** OpenAI's own September 2025 paper found that "next-token training
objectives and common leaderboards reward confident guessing over calibrated
uncertainty, so models learn to bluff" — benchmarks typically penalize
abstention, and RLHF can amplify this when human raters prefer long, detailed
answers over merely correct ones. Explicit "don't guess" instructions reduce
hallucination rates by up to 15% in some settings, and "sample disagreement
across multiple decodes" (asking the same thing several times and checking
for agreement) is described as the most reliable practical uncertainty signal
in 2026. **This directly explains this session's own real finding**: this
project's whole "refuse rather than guess" trap category (missing
counterion-binding degree, missing Davies group numbers) is testing exactly
this incentive-driven bias, and the observed inconsistency between models
(one model refusing an HLB computation, another force-fitting generic group
numbers to produce a number anyway) is consistent with a genuine, documented,
model-to-model calibration difference, not a fluke.
[LLM Hallucinations in 2026, Lakera](https://www.lakera.ai/blog/guide-to-hallucinations-in-large-language-models)

**(b) Counting repeated tokens/substructures fails for a specific, now
mechanistically-understood reason, not "the model can't count."** A 2026
mechanistic-interpretability paper found that linear probes on a model's own
internal residual stream decode the CORRECT count with near-perfect accuracy
at every layer — even at the exact layer where the model's own output
crystallizes into the WRONG number. The failure is a "format-triggered MLP
block" that overwrites the correctly-encoded internal count with a fixed
wrong answer at ~88-93% network depth. Separately, tokenization mismatch is
a real contributing factor: byte-level BPE tokenizers group characters in a
way that doesn't align with the unit being counted (e.g. counting letters
when the model's own token boundaries don't respect letter boundaries).
**This is a direct, mechanistic explanation for this session's own finding**:
both Claude and GPT independently miscounted the repeated -OCH2CH2- units in
a real SMILES string (said 7, the real answer was 8) — this isn't a fluke or
a chemistry-specific gap, it's a documented, general LLM architecture
limitation for counting repeated substructures in any long string, chemistry
notation included.
[Repeated-Token Counting Reveals a Dissociation, arXiv 2605.09239](https://arxiv.org/html/2605.09239v1),
[Counting Ability of LLMs and Impact of Tokenization, arXiv 2410.19730](https://arxiv.org/html/2410.19730v2)

**(c) Arithmetic/multi-step numeric reasoning is still a genuinely open,
unsolved problem in 2026, not a fully-closed gap.** Despite strong performance
on many reasoning benchmarks, LLMs "consistently struggle with simple
arithmetic tasks, such as the addition of multiple or large numbers," and
strong performance on a benchmark's original numbers does NOT guarantee
robustness to the same problem with different numeric values substituted in
— i.e. models may be pattern-matching a problem shape rather than genuinely
computing. Recent work also shows LLM chain-of-thought is not always a
faithful record of the actual computation happening internally (models can
use semantically-irrelevant filler tokens to boost accuracy on synthetic
reasoning tasks by up to 13 points) — meaning a model's stated derivation
steps are not fully trustworthy as a record of HOW it got a number, even when
the number itself is right.
[LLM Reasoning Is Latent, Not the Chain of Thought, arXiv 2604.15726](https://arxiv.org/abs/2604.15726),
[Testing LLM Arithmetic Reasoning Generalization, arXiv 2606.03606](https://arxiv.org/html/2606.03606v1)

**(d) Reproducibility is fundamentally probabilistic, not just an API
setting.** Even with temperature=0 and a fixed seed, different eval scores
were observed on the same prompt suite across runs, caused by "batch-dependent
floating point in the inference engine plus silent provider routing." The
working 2026 guidance is to design evaluations to "tolerate semantic
equivalence rather than byte-exact match" and bound variance rather than
expect determinism. **Directly relevant**: this project's planned "same
question, multiple fresh sessions" reproducibility axis (proposed earlier
this session, not yet run) has real literature support as a legitimate,
expected-to-show-real-variance test, not a methodology flaw if answers differ
run to run.
[temperature=0 didn't make LLM evals reproducible, dev.to](https://dev.to/marcuswwchen/temperature0-didnt-make-our-llm-evals-reproducible-5ae6),
[A statistical framework for repeatability/reproducibility of LLMs, PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12637745/)

**(e) Chart/graph numeric extraction remains a real, current, primary weak
point for vision-capable models, exactly as this session assumed when
scoping the graph-digitizer work.** A 2026 CVPR benchmark (MeasureBench)
and a CHI 2026 paper both confirm multimodal LLMs "struggle with chart data
extraction, particularly in achieving high numerical accuracy," even though
they perform well on general vision-language tasks — "numeric extraction is
identified as the primary failure mode" specifically (unit recognition is
comparatively more stable). A self-ensembling technique (running multiple
passes and combining) improves accuracy but only by ~23% relative — meaning
even the improved approach is well short of solved. This directly confirms
the design decision, made earlier this session, that our own calibrated
pixel-transform digitizer has a real, current, evidenced accuracy advantage
over asking a vision LLM to eyeball a chart.
[Do Vision-Language Models Measure Up? MeasureBench, CVPR 2026](https://arxiv.org/html/2510.26865v1),
[Self-Ensembling VLMs for Chart Data Extraction, arXiv 2605.27298](https://arxiv.org/abs/2605.27298)

**(f) Qwen/DashScope's fake-tool-call behavior (found directly in this
session's own testing) is NOT independently corroborated by any public
documentation found in this search pass.** DashScope's own docs describe the
Anthropic-compatible endpoint as supporting real "thinking and tool calling,"
with no disclosed reliability caveat. This session's own direct evidence
(raw API response showing `stop_reason: "end_turn"` instead of `"tool_use"`,
confirmed via the actual response body, not assumed) remains the strongest
and only evidence found for this specific behavior — worth stating in the
paper as this project's own original finding, not something backed by
existing published documentation, since none was found.
[DashScope Anthropic-compatible endpoint docs](https://docs.agno.com/models/providers/native/dashscope/overview)

## 3. Real, currently-documented LLM strengths — where NOT to build redundant tooling

**(a) PhD-level scientific reasoning on GPQA-Diamond is now genuinely strong
and still climbing.** As of February 2026: Gemini 3.1 Pro led at 94.3%,
Claude Opus 4.6 at 91.3%, Qwen3.5-plus at 88.4%, GPT-5.3 Codex at 81% — on a
benchmark specifically designed so non-expert PhD holders score ~34% (a real
floor). By September 2026 the top of the leaderboard (GPT-6 Astra) reached
0.960. The benchmark is "approaching saturation at the very top" but "still
clearly differentiates models in the 60-90% range." **This directly matches
and explains this session's own finding**: 87.5-95% unaugmented scores on
genuinely hard, unhinted surfactant-chemistry judgment questions is
consistent with, not an outlier from, the broader current state of frontier
PhD-level science reasoning — this project should not expect raw "does it
know the science" questions to differentiate models much further no matter
how well-designed, because this is now a broadly strong capability across
providers, confirmed by an independent, credible, non-domain-specific
benchmark.
[GPQA-Diamond Benchmark Scores, IntuitionLabs](https://intuitionlabs.ai/articles/gpqa-diamond-ai-benchmark),
[GPQA Leaderboard, llm-stats.com](https://llm-stats.com/benchmarks/gpqa)

**(b) Structured data extraction from documents is now a strong, largely-solved
general capability.** General-purpose LLMs "achieve state-of-the-art
performance on structured information extraction from semi-structured
documents, surpassing classical machine learning approaches" with "superior
generalization to unseen document layouts." Implication for this project: a
tool whose only job is "pull structured values out of a well-formatted table
a researcher pastes in" is not a strong differentiator anymore — the real
value is downstream (what real chemistry judgment happens AFTER extraction),
not the extraction step itself.
[Information Extraction from Electricity Invoices with LLMs, arXiv 2604.25927](https://arxiv.org/pdf/2604.25927)

**(c) Tool-calling itself, once a model is properly integrated with a given
protocol, is now considered a mature, reliable capability for the major
proprietary providers** — sources describe Claude's tool-call accuracy as
strong enough to be "the safe pick when a mis-call writes bad data or fires a
high-stakes action," and note real per-provider differences in tool-call
PATTERNS (not raw reliability) between Anthropic and OpenAI. This is
consistent with — not contradicted by — this session's own finding that
Gemini reliably invoked real MCP tools while Qwen did not: the general
claim in public sources is about maturity at the LEADING providers
specifically, and Qwen/DashScope's Anthropic-compatible bridge is a newer,
less-established integration path than either provider's own native API.
[Best LLM for Coding Agents in 2026, Evolink](https://evolink.ai/blog/best-llm-for-coding-agents-api-cost-reliability)
(Note: most tool-calling comparison sources found in this pass were SEO/
marketing blogs, not primary benchmarks — treat this specific claim as
lower-confidence than the sourced items above, and worth re-verifying
against a real benchmark paper if it becomes load-bearing for a paper claim.)

## 4. Synthesis — what this means for the redesign

The broader 2026 literature **confirms, rather than contradicts,** every one
of this session's own four real findings:
1. Curve-fitting-by-hand disagreement (~20% off a tool's exact segment fit)
   is consistent with documented arithmetic/multi-step-numeric weaknesses
   that persist even in otherwise-strong models (item 2c above).
2. Inconsistent refusal-vs-approximate calibration matches a well-documented,
   training-incentive-driven bias toward confident guessing over abstention
   (item 2a above) — this is not noise, it's a real, expected, reproducible
   axis to keep testing.
3. Qwen's fake-tool-call behavior has no independent public corroboration
   found — flag it honestly in the paper as this project's own original,
   directly-observed finding, not something the literature already knew.
4. The chart-reading weakness this project's graph-digitizer was built to
   address is real, current, and still far from solved even with
   state-of-the-art techniques (item 2e above) — that build decision holds up.

**The single most important redirect for the next question batch**: item 3a
(GPQA-Diamond saturation) is real and should reshape expectations. Pure
"does the model know surfactant chemistry" questions — however well-designed,
however unhinted — are testing a capability that is now broadly strong across
providers and still improving. The durable, defensible differentiator for
both the benchmark AND the toolkit's own value proposition is NOT raw
domain knowledge; it's (a) precise, deterministic multi-step numeric
computation where hand-reasoning genuinely drifts (~20% CMC disagreement is
real and reproducible), (b) calibrated refusal under genuine data gaps
(a documented, training-incentive-driven weakness, not a knowledge gap),
(c) counting/parsing repeated structural motifs in long strings (a
mechanistically-understood architecture limitation, not chemistry-specific),
and (d) chart/image-based numeric extraction (still far from solved). Build
and test toward those four, not toward "harder chemistry questions" in the
abstract — the field's own current literature says that lever is running out
of room.
