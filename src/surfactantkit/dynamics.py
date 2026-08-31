"""Translational diffusion and hydrodynamic radius (dynamic light
scattering)."""

from __future__ import annotations

K_B = 1.380649e-23  # J/K, exact SI 2019


def hydrodynamic_radius_stokes_einstein(diffusion_coefficient_cm2_per_s: float, viscosity_mPas: float, temperature_K: float = 298.15) -> float:
    """Hydrodynamic radius (nm) from the Stokes-Einstein equation:
    R_h = k_B*T / (6*pi*eta*D).

    diffusion_coefficient_cm2_per_s: translational diffusion coefficient
    in the commonly-reported DLS unit, cm^2/s (converted internally to
    m^2/s -- do not pass an m^2/s value here, this is a classic unit trap).
    viscosity_mPas: solvent viscosity in mPa.s (=cP; water is ~0.89 at 25C).
    """
    import math

    if diffusion_coefficient_cm2_per_s <= 0:
        raise ValueError("diffusion_coefficient_cm2_per_s must be positive")
    if viscosity_mPas <= 0:
        raise ValueError("viscosity_mPas must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")

    D_SI = diffusion_coefficient_cm2_per_s * 1e-4  # cm^2/s -> m^2/s
    eta_SI = viscosity_mPas * 1e-3  # mPa.s -> Pa.s
    r_h_m = (K_B * temperature_K) / (6.0 * math.pi * eta_SI * D_SI)
    return r_h_m * 1e9  # m -> nm
