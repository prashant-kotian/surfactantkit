"""Wetting, adhesion, and the capillary number (enhanced oil recovery
relevance): Young-Dupre equation and the dimensionless ratio that
connects interfacial tension reduction to actual oil-displacement
efficiency."""

from __future__ import annotations
import math


def work_of_adhesion(gamma_LV_mN_m: float, contact_angle_deg: float) -> float:
    """Young-Dupre work of adhesion (mN/m, numerically = mJ/m^2):
    W_a = gamma_LV * (1 + cos(theta)).

    gamma_LV_mN_m: liquid-vapor surface tension. contact_angle_deg: the
    measured (equilibrium) contact angle in degrees.
    """
    if gamma_LV_mN_m <= 0:
        raise ValueError("gamma_LV_mN_m must be positive")
    if not (0.0 <= contact_angle_deg <= 180.0):
        raise ValueError("contact_angle_deg must be between 0 and 180")
    theta_rad = math.radians(contact_angle_deg)
    return gamma_LV_mN_m * (1.0 + math.cos(theta_rad))


def spreading_coefficient(gamma_LV_mN_m: float, contact_angle_deg: float) -> float:
    """Spreading coefficient (mN/m) derived from Young's equation:
    S = gamma_LV * (cos(theta) - 1).

    S <= 0 always in this formulation (equality only in the limiting
    case of complete wetting, theta = 0); more negative S means poorer
    spreading. gamma_LV_mN_m and contact_angle_deg as in work_of_adhesion.
    """
    if gamma_LV_mN_m <= 0:
        raise ValueError("gamma_LV_mN_m must be positive")
    if not (0.0 <= contact_angle_deg <= 180.0):
        raise ValueError("contact_angle_deg must be between 0 and 180")
    theta_rad = math.radians(contact_angle_deg)
    return gamma_LV_mN_m * (math.cos(theta_rad) - 1.0)


def capillary_number(viscosity_mPas: float, velocity_m_per_s: float, interfacial_tension_mN_m: float) -> float:
    """Capillary number (dimensionless): Ca = (viscosity * velocity) /
    interfacial_tension -- the ratio of viscous to capillary forces that
    governs oil-displacement efficiency in enhanced oil recovery. Below
    Ca ~ 1e-5, flow is capillary-dominated and residual oil is trapped;
    surfactants raise Ca primarily by driving interfacial tension toward
    ultra-low values.

    viscosity_mPas: displacing-fluid viscosity in mPa.s (=cP).
    velocity_m_per_s: Darcy/interstitial velocity in m/s.
    interfacial_tension_mN_m: oil-water interfacial tension in mN/m.
    """
    if viscosity_mPas <= 0 or velocity_m_per_s <= 0 or interfacial_tension_mN_m <= 0:
        raise ValueError("viscosity, velocity, and interfacial tension must all be positive")
    eta_SI = viscosity_mPas * 1e-3  # mPa.s -> Pa.s
    gamma_SI = interfacial_tension_mN_m * 1e-3  # mN/m -> N/m
    return (eta_SI * velocity_m_per_s) / gamma_SI
