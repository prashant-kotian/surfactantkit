# SurfactantKit

Surfactant and interfacial-science theory as a tested, importable Python library: Clint ideal mixing, Rubingh regular-solution theory (and its Rosen extension to mixed monolayers), the Gibbs and Szyszkowski adsorption relations, HLB (Griffin's and Davies' methods), the critical packing parameter and aggregation number (Tanford's formulas), full micellization thermodynamics (ΔG/ΔH/ΔS, counterion binding), the Corrin-Harkins salt-CMC relation, Debye screening length, zeta potential (Henry equation), hydrodynamic radius (Stokes-Einstein), wetting/adhesion (Young-Dupré, spreading coefficient), the capillary number (enhanced oil recovery), and micellar solubilization (Molar Solubilization Ratio). Built to stop AI assistants (and researchers) from hallucinating these values, and to be genuinely correct rather than plausible-looking.

## Status

Early alpha, but functionally broad: 11 theory modules, 25 MCP tools, 106 passing tests. Core theory is implemented and tested, and SurfMCP (`surfactantkit.mcp_server`) exposes all of it as MCP tools for Claude Desktop, Cursor, and other MCP-aware AI assistants.

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

This exposes 25 tools covering Clint/Rubingh/Rosen mixed-system theory, Gibbs and Szyszkowski adsorption, HLB (both methods), critical packing parameter and aggregation number, Corrin-Harkins salt-CMC prediction, Debye screening length, zeta potential, hydrodynamic radius, full micellization thermodynamics (ΔG/ΔH/ΔS, counterion binding), wetting/adhesion, the EOR capillary number, and micellar solubilization. Every tool returns units explicitly in its response — the point of this server is to stop unit confusion (mN/m vs dyn/cm, cP vs Pa·s, cm²/s vs m²/s, etc.), so a bare unlabeled number would defeat the purpose. Tools that would otherwise require guessing an unverified value (a Davies group number with no cited source, a Rubingh equation with no valid root, an electrophoretic-mobility regime) raise an explicit error instead of returning a plausible-looking wrong number.

Built against MCP spec 2026-07-28 / Python SDK v2.

## Validation

See `tests/test_mixed_micelle.py` and `literature_validation_notes.md`. Clint and Rubingh are checked against eight independent published binary surfactant systems (numbers pulled directly from source tables, not paraphrased). Clint ideal CMC and the Rubingh solver's micellar mole fraction (x1) match literature values essentially exactly across the board; the Rubingh beta parameter matches sign and order of magnitude but usually not to the decimal — investigated and documented, not swept under the rug (see the notes file for why). The Rubingh solver is also checked with mathematically guaranteed round-trip tests independent of any external source.

HLB and CPP (see `tests/test_hlb_cpp.py`) use constants verified against cited sources before being hardcoded: Davies' (1957) group numbers, and Tanford's chain-volume/length formulas (cross-checked internally — the commonly-cited shorthand `27.4 + 26.9*nc` is algebraically derived from, and matches exactly, the more detailed per-CH2/CH3 formula). The Davies table intentionally omits groups with no confirmed value (e.g. quaternary ammonium, amide) rather than guessing — `hlb_davies()` raises on an unknown group instead of silently returning a wrong number. `tanford_tail_volume` matches a real independent re-derivation of the same Tanford formula (Bales et al., *J. Phys. Chem. B* 1998) to within 0.06% for SDS; `critical_packing_parameter` matches a real sodium-alkyl-sulfate series (Nagarajan, *Langmuir* 2002) to within ~1% across five chain lengths.

**Electrostatics, dynamics, and thermodynamics** (see `tests/test_electrostatics_dynamics.py`, `tests/test_thermodynamics.py`, `tests/test_rosen_and_corrin_harkins.py`, and `literature_validation_notes.md` rounds 2-4 for the real-paper checks below) use CODATA/SI-2019 physical constants and formulas checked against real sources before implementation. `debye_length` reproduces both the classic textbook reference value (~0.96 nm at 0.1 M, 1:1 electrolyte, 25°C) and a 2026 review's stated physiological-saline value (0.785 nm at 0.15 M) exactly. `zeta_potential_henry`'s Smoluchowski case is confirmed algebraically identical to a real paper's own equation (Khademi et al., *Langmuir* 2017) and reproduces its reported zeta potentials exactly in a formula round-trip. `hydrodynamic_radius_stokes_einstein` was tested against a real paper reporting both a DOSY-NMR diffusion coefficient and a DLS hydrodynamic radius for the same system (Milone et al.) — the formula correctly reproduces what Stokes-Einstein predicts from that specific D (a small nanometer-scale value), which does *not* match the paper's much larger DLS-derived radius; documented as an expected physical-chemistry caveat (DOSY reports a population-averaged D dominated by fast free-monomer diffusion, not the large-aggregate D) rather than a formula error — a useful negative control, not a failure. `gibbs_free_energy_micellization` matches a real paper (Fu et al., *RSC Advances* 2019) to within 0.03%; `vant_hoff_enthalpy`'s two-point method was found to genuinely diverge from a paper's own polynomial-derivative ΔH on the same real CMC(T) data — documented honestly as a real, now-confirmed limitation of the two-point assumption (already flagged in the function's own docstring), not glossed over. The geometric aggregation number for a C12 chain lands at 55.5 — squarely inside the commonly-cited literature range (~55-70) for SDS/DTAB-scale surfactants. The Rosen monolayer functions are verified to produce results identical to the underlying Rubingh functions, since Rosen's theory is mathematically the same regular-solution equation applied to surface-tension data instead of CMC data — a distinction in *meaning*, not math, and exactly the kind of subtlety an AI assistant is likely to blur.

**Szyszkowski, adsorption, wetting, and solubilization** (see `tests/test_szyszkowski_wetting_solubilization.py` and `literature_validation_notes.md` round 3) — `gibbs_gamma_max` matches four real ionic surfactant systems (SDS, CTAB, AOT, CPC; Shah/Das/Bhattarai 2025) to within 0.5%, and `gibbs_a_min` matches five real systems across two papers to within 0.4%. `work_of_adhesion` matches six real systems across two papers (including an EOR interfacial-tension-based case) to within 0.34%; `spreading_coefficient` agreement is tight in absolute terms (looks like a larger percentage only because the tested S values are close to zero). The Szyszkowski surface-tension equation is separately cross-checked against the independently-implemented Gibbs adsorption module: at high concentration, the Szyszkowski curve's numerical slope, fed back through `gibbs_gamma_max`, recovers the same Γmax used to construct it to within 0.02% — two separately-coded physics relations agreeing exactly where theory says they must. Young-Dupré/spreading-coefficient boundary conditions (complete wetting, complete non-wetting) are exact by construction. `molar_solubilization_ratio` has round-trip and boundary-condition tests but no confirmed real-paper match yet — two literature candidates were found to use a below-saturation partitioning regime that doesn't map onto this function's above-saturation convention (documented in the notes file rather than forced through with a mismatched substitution); still open. The micelle-water partition coefficient (Km), despite being a commonly-reported companion to the Molar Solubilization Ratio, is deliberately **not** implemented — a literature check found its exact formula is study-specific with no single dominant convention, and guessing one risks silently reproducing the wrong definition.

## License

MIT
