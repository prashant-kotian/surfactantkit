"""SurfactantKit: mixed-surfactant micellization theory as a tested, importable library.

Implements Clint ideal mixing, Rubingh regular-solution theory, and the
Gibbs adsorption isotherm -- the calculations an AI assistant or a
researcher needs to reason correctly about mixed-surfactant systems
instead of hallucinating numbers.
"""

from .mixed_micelle import (
    clint_ideal_cmc,
    solve_rubingh_x,
    rubingh_beta,
    activity_coefficients,
    excess_free_energy,
)
from .adsorption import gibbs_gamma_max, gibbs_a_min

__version__ = "0.1.0"

__all__ = [
    "clint_ideal_cmc",
    "solve_rubingh_x",
    "rubingh_beta",
    "activity_coefficients",
    "excess_free_energy",
    "gibbs_gamma_max",
    "gibbs_a_min",
]
