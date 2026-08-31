# SurfactantKit

Mixed-surfactant micellization theory as a tested, importable Python library: Clint ideal mixing, Rubingh regular-solution theory, and the Gibbs adsorption isotherm. Built to stop AI assistants (and researchers) from hallucinating CMC and interaction-parameter values, and to be genuinely correct rather than plausible-looking.

## Status

Early alpha. Core theory (`surfactantkit.mixed_micelle`, `surfactantkit.adsorption`) is implemented and tested. An MCP server exposing these functions as tools for Claude/Cursor/other AI assistants is planned next.

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

## Validation

See `tests/test_mixed_micelle.py`. The Clint ideal CMC formula is checked against real literature values (DTAB-SDS system) to an exact match. The Rubingh solver is checked with mathematically guaranteed round-trip tests, plus a qualitative sign check (synergistic vs. antagonistic) against the same literature system. We deliberately do not assert an exact literature beta value we could not independently reproduce with confidence — see the test module docstring for why.

## License

MIT
