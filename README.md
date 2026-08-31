# SurfactantKit

Surfactant and interfacial-science theory as a tested, importable Python library: Clint ideal mixing, Rubingh regular-solution theory (and its Rosen extension to mixed monolayers), the Gibbs adsorption isotherm, HLB (Griffin's and Davies' methods), the critical packing parameter and aggregation number (Tanford's formulas), full micellization thermodynamics (ΔG/ΔH/ΔS, counterion binding), the Corrin-Harkins salt-CMC relation, Debye screening length, zeta potential (Henry equation), and hydrodynamic radius (Stokes-Einstein). Built to stop AI assistants (and researchers) from hallucinating these values, and to be genuinely correct rather than plausible-looking.

## Status

Early alpha, but functionally broad: 9 theory modules, 20 MCP tools, 86 passing tests. Core theory (`surfactantkit.mixed_micelle`, `surfactantkit.adsorption`, `surfactantkit.hlb`, `surfactantkit.cpp`, `surfactantkit.electrostatics`, `surfactantkit.dynamics`, `surfactantkit.thermodynamics`) is implemented and tested, and SurfMCP (`surfactantkit.mcp_server`) exposes all of it as MCP tools for Claude Desktop, Cursor, and other MCP-aware AI assistants.

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

This exposes 20 tools covering Clint/Rubingh/Rosen mixed-system theory, Gibbs adsorption, HLB (both methods), critical packing parameter and aggregation number, Corrin-Harkins salt-CMC prediction, Debye screening length, zeta potential, hydrodynamic radius, and full micellization thermodynamics (ΔG/ΔH/ΔS, counterion binding). Every tool returns units explicitly in its response — the point of this server is to stop unit confusion (mN/m vs dyn/cm, cP vs Pa·s, cm²/s vs m²/s, etc.), so a bare unlabeled number would defeat the purpose. Tools that would otherwise require guessing an unverified value (a Davies group number with no cited source, a Rubingh equation with no valid root, an electrophoretic-mobility regime) raise an explicit error instead of returning a plausible-looking wrong number.

Built against MCP spec 2026-07-28 / Python SDK v2.

## Validation

See `tests/test_mixed_micelle.py` and `literature_validation_notes.md`. Clint and Rubingh are checked against eight independent published binary surfactant systems (numbers pulled directly from source tables, not paraphrased). Clint ideal CMC and the Rubingh solver's micellar mole fraction (x1) match literature values essentially exactly across the board; the Rubingh beta parameter matches sign and order of magnitude but usually not to the decimal — investigated and documented, not swept under the rug (see the notes file for why). The Rubingh solver is also checked with mathematically guaranteed round-trip tests independent of any external source.

HLB and CPP (see `tests/test_hlb_cpp.py`) use constants verified against cited sources before being hardcoded: Davies' (1957) group numbers, and Tanford's chain-volume/length formulas (cross-checked internally — the commonly-cited shorthand `27.4 + 26.9*nc` is algebraically derived from, and matches exactly, the more detailed per-CH2/CH3 formula). The Davies table intentionally omits groups with no confirmed value (e.g. quaternary ammonium, amide) rather than guessing — `hlb_davies()` raises on an unknown group instead of silently returning a wrong number.

**Electrostatics, dynamics, and thermodynamics** (see `tests/test_electrostatics_dynamics.py`, `tests/test_thermodynamics.py`, `tests/test_rosen_and_corrin_harkins.py`) use CODATA/SI-2019 physical constants and formulas checked against real sources before implementation. The Debye length calculation reproduces the classic textbook reference value (~0.96 nm for 0.1 M 1:1 electrolyte at 25°C) to within 0.5%. The geometric aggregation number for a C12 chain lands at 55.5 — squarely inside the commonly-cited literature range (~55-70) for SDS/DTAB-scale surfactants. Functions without a clean external reference point (zeta potential, hydrodynamic radius, van't Hoff enthalpy, Corrin-Harkins) are checked with mathematically guaranteed round-trip tests instead. The Rosen monolayer functions are verified to produce results identical to the underlying Rubingh functions, since Rosen's theory is mathematically the same regular-solution equation applied to surface-tension data instead of CMC data — a distinction in *meaning*, not math, and exactly the kind of subtlety an AI assistant is likely to blur.

## License

MIT
