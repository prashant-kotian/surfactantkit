# SurfactantKit

Surfactant and interfacial-science theory as a tested, importable Python library: Clint ideal mixing, Rubingh regular-solution theory (and its Rosen extension to mixed monolayers, plus a real BIC-based orchestrator that fits Clint/Rubingh/EOMMM to a raw composition series and picks the winner), the Gibbs and Szyszkowski adsorption relations (including a Langmuir-vs-Frumkin BIC model selector), HLB (Griffin's, Davies', and Guo/Rong/Ying's methods), HLD-NAC (optimal salinity, salinity-scan fitting), the critical packing parameter and aggregation number (Tanford's formulas, plus Nagarajan's chain-length-dependent electrostatic headgroup model and a real CPP-geometry axial-ratio estimator), full micellization thermodynamics (ΔG/ΔH/ΔS, counterion binding, mass-action and multi-point van't Hoff fits), the Corrin-Harkins salt-CMC relation, Debye screening length, zeta potential (Henry equation and the O'Brien-White relaxation correction for high zeta), hydrodynamic radius (Stokes-Einstein and Perrin shape-corrected), wetting/adhesion (Young-Dupré, OWRK and van Oss-Chaudhury-Good surface-energy decomposition, with a real cited standard-probe-liquid lookup table), the capillary number (enhanced oil recovery), micellar solubilization (Molar Solubilization Ratio, the micelle-water partition coefficient, and a QSPR fallback estimator for compounds with no measured solubility), raw-curve extraction (CMC from tensiometry/conductivity, aggregation number from fluorescence quenching/SLS/DLS/analytical ultracentrifugation), SMILES-based structural classification (charge type, including a real pH-conditional resolver, and structural family: gemini/dimeric vs. glycolipid biosurfactant), and a full autonomous orchestrator that derives everything it can directly from (SMILES, raw curve) alone, naming every real gap rather than guessing. Built to stop AI assistants (and researchers) from hallucinating these values, and to be genuinely correct rather than plausible-looking.

## Status

Early alpha, but functionally broad: 14 theory/capability modules, 67 MCP tools, 438 passing tests (current as of 2026-09-15; re-verify with `pytest tests/ --collect-only -q` and `grep -c "^@mcp.tool" src/surfactantkit/mcp_server.py` rather than trusting this number blind if it's been a while). Core theory is implemented and tested, and SurfMCP (`surfactantkit.mcp_server`) exposes all of it as MCP tools for Claude Desktop, Cursor, and other MCP-aware AI assistants. Full build history and per-capability validation detail lives in `ROADMAP.md` and `literature_validation_notes.md`; do not assume either "Status" line above is exhaustive detail — they are current counts, not a substitute for reading those files.

## Why

No open-source library implements Rubingh/Clint mixed-surfactant theory or the Gibbs adsorption isotherm as reusable, tested code — this kind of calculation is still routinely done by hand or in one-off spreadsheets. This library grew out of a real thermodynamic calculator built for a PhD project on amidoamine-derived gemini cationic surfactants, validated against real conductometry and tensiometry data, and generalized here into a standalone tool.

## Install (development)

```bash
pip install -e ".[dev]"
pytest
```

## Quick example

```python
from surfactantkit import clint_ideal_cmc, solve_rubingh_x, rubingh_beta

# DTAB-SDS mixed system, SDS-rich composition
cmc_id = clint_ideal_cmc(alpha1=0.25, cmc1=14.80, cmc2=8.00)  # mM
x1 = solve_rubingh_x(alpha1=0.25, cmc_mix=6.011, cmc1=14.80, cmc2=8.00)
beta = rubingh_beta(x1, alpha1=0.25, cmc_mix=6.011, cmc1=14.80)
```

## SurfMCP: using this from Claude Desktop / Cursor / other AI assistants

Install with the `mcp` extra, then point your MCP client at the `surfactantkit-mcp` entry point:

```bash
pip install -e ".[mcp]"
```

Claude Desktop config (`claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "surfactantkit": {
      "command": "surfactantkit-mcp"
    }
  }
}
```

This exposes 67 tools (current as of 2026-09-15; re-count with `grep -c "^@mcp.tool" src/surfactantkit/mcp_server.py` rather than trusting this number blind later) covering every capability listed at the top of this file — mixed-system theory and its own BIC model-selection orchestrator, Gibbs/Szyszkowski/Frumkin adsorption with BIC isotherm selection, HLB (three methods) and HLD-NAC, critical packing parameter/aggregation number/axial-ratio geometry, electrostatics, dynamics, full micellization thermodynamics, wetting/adhesion including a cited standard-probe-liquid lookup, the EOR capillary number, micellar solubilization including a QSPR fallback, raw-curve/raw-dataset extraction, SMILES-based structural classification, and the full SMILES+curve autonomous orchestrator. Every tool returns units explicitly in its response — the point of this server is to stop unit confusion (mN/m vs dyn/cm, cP vs Pa·s, cm²/s vs m²/s, etc.), so a bare unlabeled number would defeat the purpose. Tools that would otherwise require guessing an unverified value (a Davies group number with no cited source, a Rubingh equation with no valid root, an electrophoretic-mobility regime, a solution pH, an electrolyte condition) raise an explicit error or return a named, disclosed gap instead of a plausible-looking wrong number.

Built against MCP spec 2026-07-28 / Python SDK v2.

## Validation

See `tests/test_mixed_micelle.py` and `literature_validation_notes.md`. Clint and Rubingh are checked against eight independent published binary surfactant systems (numbers pulled directly from source tables, not paraphrased). Clint ideal CMC and the Rubingh solver's micellar mole fraction (x1) match literature values essentially exactly across the board; the Rubingh beta parameter matches sign and order of magnitude but usually not to the decimal — investigated and documented, not swept under the rug (see the notes file for why). The Rubingh solver is also checked with mathematically guaranteed round-trip tests independent of any external source.

HLB and CPP (see `tests/test_hlb_cpp.py`) use constants verified against cited sources before being hardcoded: Davies' (1957) group numbers, and Tanford's chain-volume/length formulas (cross-checked internally — the commonly-cited shorthand `27.4 + 26.9*nc` is algebraically derived from, and matches exactly, the more detailed per-CH2/CH3 formula). The Davies table intentionally omits groups with no confirmed value (e.g. quaternary ammonium, amide) rather than guessing — `hlb_davies()` raises on an unknown group instead of silently returning a wrong number. `tanford_tail_volume` matches a real independent re-derivation of the same Tanford formula (Bales et al., *J. Phys. Chem. B* 1998) to within 0.06% for SDS; `critical_packing_parameter` matches a real sodium-alkyl-sulfate series (Nagarajan, *Langmuir* 2002) to within ~1% across five chain lengths.

**Electrostatics, dynamics, and thermodynamics** (see `tests/test_electrostatics_dynamics.py`, `tests/test_thermodynamics.py`, `tests/test_rosen_and_corrin_harkins.py`, and `literature_validation_notes.md` rounds 2-4 for the real-paper checks below) use CODATA/SI-2019 physical constants and formulas checked against real sources before implementation. `debye_length` reproduces both the classic textbook reference value (~0.96 nm at 0.1 M, 1:1 electrolyte, 25°C) and a 2026 review's stated physiological-saline value (0.785 nm at 0.15 M) exactly. `zeta_potential_henry`'s Smoluchowski case is confirmed algebraically identical to a real paper's own equation (Khademi et al., *Langmuir* 2017) and reproduces its reported zeta potentials exactly in a formula round-trip. `hydrodynamic_radius_stokes_einstein` was tested against a real paper reporting both a DOSY-NMR diffusion coefficient and a DLS hydrodynamic radius for the same system (Milone et al.) — the formula correctly reproduces what Stokes-Einstein predicts from that specific D (a small nanometer-scale value), which does *not* match the paper's much larger DLS-derived radius; documented as an expected physical-chemistry caveat (DOSY reports a population-averaged D dominated by fast free-monomer diffusion, not the large-aggregate D) rather than a formula error — a useful negative control, not a failure. `gibbs_free_energy_micellization` matches a real paper (Fu et al., *RSC Advances* 2019) to within 0.03%; `vant_hoff_enthalpy`'s two-point method was found to genuinely diverge from a paper's own polynomial-derivative ΔH on the same real CMC(T) data — documented honestly as a real, now-confirmed limitation of the two-point assumption (already flagged in the function's own docstring), not glossed over. The geometric aggregation number for a C12 chain lands at 55.5 — squarely inside the commonly-cited literature range (~55-70) for SDS/DTAB-scale surfactants. The Rosen monolayer functions are verified to produce results identical to the underlying Rubingh functions, since Rosen's theory is mathematically the same regular-solution equation applied to surface-tension data instead of CMC data — a distinction in *meaning*, not math, and exactly the kind of subtlety an AI assistant is likely to blur.

**Szyszkowski, adsorption, wetting, and solubilization** (see `tests/test_szyszkowski_wetting_solubilization.py` and `literature_validation_notes.md` round 3) — `gibbs_gamma_max` matches four real ionic surfactant systems (SDS, CTAB, AOT, CPC; Shah/Das/Bhattarai 2025) to within 0.5%, and `gibbs_a_min` matches five real systems across two papers to within 0.4%. `work_of_adhesion` matches six real systems across two papers (including an EOR interfacial-tension-based case) to within 0.34%; `spreading_coefficient` agreement is tight in absolute terms (looks like a larger percentage only because the tested S values are close to zero). The Szyszkowski surface-tension equation is separately cross-checked against the independently-implemented Gibbs adsorption module: at high concentration, the Szyszkowski curve's numerical slope, fed back through `gibbs_gamma_max`, recovers the same Γmax used to construct it to within 0.02% — two separately-coded physics relations agreeing exactly where theory says they must. Young-Dupré/spreading-coefficient boundary conditions (complete wetting, complete non-wetting) are exact by construction. `molar_solubilization_ratio` has round-trip and boundary-condition tests but no confirmed real-paper match yet — two literature candidates were found to use a below-saturation partitioning regime that doesn't map onto this function's above-saturation convention (documented in the notes file rather than forced through with a mismatched substitution); still open. **The micelle-water partition coefficient (Km) IS implemented** (`micelle_water_partition_coefficient`, added 2026-09-10, mole-fraction-basis convention specifically) — a literature check found its exact formula is study-specific with no single dominant convention (some papers instead use a molar-concentration-ratio basis, giving a numerically different value for the same raw data), so that ambiguity is prominently disclosed in the function's own docstring rather than silently picking a convention or omitting the tool entirely.

**Newer capabilities not yet given their own prose paragraph here** (all real, tested, MCP-wired — see `ROADMAP.md` for full build detail and `BOTTLENECK_RESOLUTION_PLAN.md`/`GENUINE_BOTTLENECK_AUDIT.md` in the `benchmark/paper3_groundzero/` folder of the companion PhD-research repo for the most recent additions): the `hld.py` module (optimal salinity, salinity-scan fitting); `mixture_model_selection.py` (BIC-based Clint/Rubingh/EOMMM selection on a raw composition series); `orchestrate.py` (a full autonomous SMILES+raw-curve derivation pipeline that names every real gap instead of guessing); `classify.py` (SMILES-based charge-type and structural-family classification, including a pH-conditional Henderson-Hasselbalch resolver); `adsorption.py`'s `frumkin_fit_K_and_a`/`select_isotherm_model` (joint nonlinear Frumkin fit and BIC-based Langmuir-vs-Frumkin selection); `cpp.py`'s `estimate_axial_ratio_from_cpp_geometry` (real prolate-ellipsoid geometry for Perrin shape correction); `curve_analysis.py`'s `estimate_partial_specific_volume_from_tail_and_headgroup` (Tanford-tail-volume-based v_bar, explicitly partial by design); `solubilization.py`'s `estimate_intrinsic_water_solubility_qspr` (ESOL QSPR fallback, disclosed as a coarse estimate, never a measurement); and `wetting.py`'s `get_vocg_standard_liquid`/`get_owens_wendt_standard_liquid` (real, cited standard-probe-liquid values for water/diiodomethane/glycerol/formamide).

## License

MIT
