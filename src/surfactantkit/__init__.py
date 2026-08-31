"""SurfactantKit: surfactant and interfacial-science theory as a tested,
importable library.

Implements Clint ideal mixing, Rubingh/Rosen regular-solution theory,
the Gibbs adsorption isotherm, HLB, the critical packing parameter,
electrostatics (Debye length, zeta potential), dynamics (hydrodynamic
radius), and full micellization thermodynamics -- the calculations an
AI assistant or a researcher needs to reason correctly about surfactant
systems instead of hallucinating numbers.
"""

from .mixed_micelle import (
    clint_ideal_cmc,
    solve_rubingh_x,
    rubingh_beta,
    activity_coefficients,
    excess_free_energy,
    solve_rosen_monolayer_x,
    rosen_beta_sigma,
    corrin_harkins_predict_cmc,
)
from .adsorption import gibbs_gamma_max, gibbs_a_min
from .hlb import hlb_griffin, hlb_davies
from .cpp import (
    tanford_tail_volume,
    tanford_critical_length,
    critical_packing_parameter,
    classify_aggregate_morphology,
    aggregation_number_spherical,
)
from .electrostatics import (
    ionic_strength,
    debye_length,
    henry_function,
    zeta_potential_henry,
)
from .dynamics import hydrodynamic_radius_stokes_einstein
from .thermodynamics import (
    cmc_to_mole_fraction,
    counterion_binding_degree,
    gibbs_free_energy_micellization,
    vant_hoff_enthalpy,
    entropy_micellization,
)

__version__ = "0.2.0"

__all__ = [
    "clint_ideal_cmc",
    "solve_rubingh_x",
    "rubingh_beta",
    "activity_coefficients",
    "excess_free_energy",
    "solve_rosen_monolayer_x",
    "rosen_beta_sigma",
    "corrin_harkins_predict_cmc",
    "gibbs_gamma_max",
    "gibbs_a_min",
    "hlb_griffin",
    "hlb_davies",
    "tanford_tail_volume",
    "tanford_critical_length",
    "critical_packing_parameter",
    "classify_aggregate_morphology",
    "aggregation_number_spherical",
    "ionic_strength",
    "debye_length",
    "henry_function",
    "zeta_potential_henry",
    "hydrodynamic_radius_stokes_einstein",
    "cmc_to_mole_fraction",
    "counterion_binding_degree",
    "gibbs_free_energy_micellization",
    "vant_hoff_enthalpy",
    "entropy_micellization",
]
