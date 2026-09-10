"""Wetting, adhesion, and the capillary number (enhanced oil recovery
relevance): Young-Dupre equation and the dimensionless ratio that
connects interfacial tension reduction to actual oil-displacement
efficiency."""

from __future__ import annotations
import math
from dataclasses import dataclass


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


@dataclass
class OwensWendtResult:
    gamma_s_dispersive_mN_m: float
    gamma_s_polar_mN_m: float
    gamma_s_total_mN_m: float
    r_squared: float
    n_liquids: int


def owens_wendt_solid_surface_energy(
    contact_angles_deg: list[float],
    liquid_gamma_dispersive_mN_m: list[float],
    liquid_gamma_polar_mN_m: list[float],
) -> OwensWendtResult:
    """Decompose a SOLID's surface energy into dispersive and polar
    components (OWRK: Owens & Wendt, J. Appl. Polym. Sci. 13 (1969)
    1741-1747; Rabel, Farbe Lack 77 (1971) 997; Kaelble, J. Adhesion 2
    (1970) 66) -- a real, sourced ALTERNATIVE to work_of_adhesion's
    single scalar number, needed only when DECOMPOSING surface energy
    into dispersive/polar components matters (not a more-correct answer
    to the same question work_of_adhesion answers -- see this project's
    own literature review, Georgiev et al., Colloids and Interfaces
    8(6):62 (2024), read in full).

    Combines Young-Dupre's work of adhesion (this module's own
    work_of_adhesion, W = gamma_L*(1+cos(theta))) with the OWRK mixing
    rule:

        W_i = 2*[sqrt(gamma_S^d * gamma_L,i^d) + sqrt(gamma_S^p * gamma_L,i^p)]

    which is LINEAR in the two unknowns X=sqrt(gamma_S^d),
    Y=sqrt(gamma_S^p) for known per-liquid dispersive/polar components
    (gamma_L,i^d, gamma_L,i^p) -- solved here via ordinary least squares
    (closed form, ordinary 2-parameter normal equations) across the
    given test liquids (needs at least 2, the standard minimum for this
    method; 3+ gives a real r_squared diagnostic instead of an exact
    fit). gamma_S_total = gamma_S^d + gamma_S^p.

    contact_angles_deg, liquid_gamma_dispersive_mN_m, and
    liquid_gamma_polar_mN_m must be the same length and in matching
    order (one entry per test liquid) -- real, literature-tabulated
    values for common test liquids (e.g. water: gamma^d~21.8,
    gamma^p~51.0 mN/m; diiodomethane: gamma^d~50.8, gamma^p~0 mN/m) must
    be supplied by the caller, not guessed here.
    """
    if not (len(contact_angles_deg) == len(liquid_gamma_dispersive_mN_m) == len(liquid_gamma_polar_mN_m)):
        raise ValueError("contact_angles_deg, liquid_gamma_dispersive_mN_m, and liquid_gamma_polar_mN_m "
                          "must all be the same length")
    n = len(contact_angles_deg)
    if n < 2:
        raise ValueError("need at least 2 test liquids to solve for 2 unknowns (gamma_S^d, gamma_S^p)")
    if any(g <= 0 for g in liquid_gamma_dispersive_mN_m):
        raise ValueError("liquid_gamma_dispersive_mN_m values must all be positive")
    if any(g < 0 for g in liquid_gamma_polar_mN_m):
        raise ValueError("liquid_gamma_polar_mN_m values must all be non-negative")

    gamma_L_total = [d + p for d, p in zip(liquid_gamma_dispersive_mN_m, liquid_gamma_polar_mN_m)]
    w = [work_of_adhesion(gt, theta) for gt, theta in zip(gamma_L_total, contact_angles_deg)]
    sqrt_d = [math.sqrt(d) for d in liquid_gamma_dispersive_mN_m]
    sqrt_p = [math.sqrt(p) for p in liquid_gamma_polar_mN_m]
    # y_i = W_i / (2*sqrt_d_i) = X + Y*(sqrt_p_i/sqrt_d_i)  -- OLS on this linear form
    y = [w[i] / (2.0 * sqrt_d[i]) for i in range(n)]
    x = [sqrt_p[i] / sqrt_d[i] for i in range(n)]

    mean_x = sum(x) / n
    mean_y = sum(y) / n
    sxx = sum((xi - mean_x) ** 2 for xi in x)
    if sxx == 0:
        slope = 0.0
    else:
        sxy = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
        slope = sxy / sxx
    intercept = mean_y - slope * mean_x  # = X = sqrt(gamma_S^d)

    gamma_s_d = intercept ** 2
    gamma_s_p = slope ** 2

    fitted = [intercept + slope * xi for xi in x]
    ss_res = sum((yi - fi) ** 2 for yi, fi in zip(y, fitted))
    ss_tot = sum((yi - mean_y) ** 2 for yi in y)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

    return OwensWendtResult(
        gamma_s_dispersive_mN_m=gamma_s_d,
        gamma_s_polar_mN_m=gamma_s_p,
        gamma_s_total_mN_m=gamma_s_d + gamma_s_p,
        r_squared=r_squared,
        n_liquids=n,
    )


def _solve_3x3(matrix: list[list[float]], vector: list[float]) -> list[float]:
    """Gaussian elimination with partial pivoting for a 3x3 linear
    system -- local to this module (the identical helper already exists
    in thermodynamics.py; each module stays self-contained rather than
    sharing this small a utility, same pattern as mixed_micelle.py's own
    local golden-section helper)."""
    a = [row[:] + [vector[i]] for i, row in enumerate(matrix)]
    n = 3
    for col in range(n):
        pivot_row = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot_row][col]) < 1e-15:
            raise ValueError("the 3 test liquids are not independent enough to solve this system "
                              "(singular matrix) -- use 3 liquids with genuinely different "
                              "LW/acid/base compositions")
        a[col], a[pivot_row] = a[pivot_row], a[col]
        for r in range(col + 1, n):
            factor = a[r][col] / a[col][col]
            for c in range(col, n + 1):
                a[r][c] -= factor * a[col][c]
    x = [0.0, 0.0, 0.0]
    for i in reversed(range(n)):
        x[i] = (a[i][n] - sum(a[i][j] * x[j] for j in range(i + 1, n))) / a[i][i]
    return x


@dataclass
class VanOssChaudhuryGoodResult:
    gamma_s_lw_mN_m: float
    gamma_s_acid_mN_m: float
    gamma_s_base_mN_m: float
    gamma_s_total_mN_m: float
    n_liquids: int


def van_oss_chaudhury_good_solid_surface_energy(
    contact_angles_deg: list[float],
    liquid_gamma_lw_mN_m: list[float],
    liquid_gamma_acid_mN_m: list[float],
    liquid_gamma_base_mN_m: list[float],
) -> VanOssChaudhuryGoodResult:
    """Decompose a SOLID's surface energy into Lifshitz-van der Waals
    (LW) and Lewis acid/base components -- van Oss, Chaudhury & Good's
    three-liquid method (Chem. Rev. 88 (1988) 927-941), a real, sourced
    ALTERNATIVE to work_of_adhesion's single scalar number and to
    owens_wendt_solid_surface_energy's simpler dispersive/polar split,
    needed when the acid-base (electron-donor/electron-acceptor)
    character of the surface matters, not just polarity in general.

    Combines Young-Dupre's work of adhesion with the vOCG mixing rule:

        W_i = 2*[sqrt(gamma_S^LW*gamma_L,i^LW) + sqrt(gamma_S^a*gamma_L,i^b)
                 + sqrt(gamma_S^b*gamma_L,i^a)]

    NOTE the cross terms: the solid's acid (electron-acceptor) component
    pairs with each liquid's BASE component, and vice versa -- not
    acid-acid. This is LINEAR in the three unknowns X=sqrt(gamma_S^LW),
    Y=sqrt(gamma_S^a), Z=sqrt(gamma_S^b), solved here via a closed-form
    3x3 linear solve (exactly 3 test liquids -- the standard minimum for
    this method: one purely apolar liquid, e.g. diiodomethane, plus two
    polar liquids with known, independent acid/base character, e.g.
    water and glycerol -- do not use 3 liquids that are not genuinely
    independent in LW/acid/base composition, the solve will raise if
    they are). gamma_S_total = gamma_S^LW + 2*sqrt(gamma_S^a*gamma_S^b).

    contact_angles_deg, liquid_gamma_lw_mN_m, liquid_gamma_acid_mN_m,
    and liquid_gamma_base_mN_m must all be length exactly 3 and in
    matching order -- real, literature-tabulated per-liquid values
    (e.g. water: LW~21.8, acid~25.5, base~25.5 mN/m; diiodomethane:
    LW~50.8, acid~0, base~0 mN/m) must be supplied by the caller, not
    guessed here.
    """
    lengths = {len(contact_angles_deg), len(liquid_gamma_lw_mN_m), len(liquid_gamma_acid_mN_m), len(liquid_gamma_base_mN_m)}
    if lengths != {3}:
        raise ValueError("contact_angles_deg, liquid_gamma_lw_mN_m, liquid_gamma_acid_mN_m, and "
                          "liquid_gamma_base_mN_m must all be length exactly 3 (the standard "
                          "3-test-liquid minimum for this method)")
    if any(g <= 0 for g in liquid_gamma_lw_mN_m):
        raise ValueError("liquid_gamma_lw_mN_m values must all be positive")
    if any(g < 0 for g in liquid_gamma_acid_mN_m) or any(g < 0 for g in liquid_gamma_base_mN_m):
        raise ValueError("liquid_gamma_acid_mN_m and liquid_gamma_base_mN_m values must all be non-negative")

    gamma_L_total = [
        lw + 2.0 * math.sqrt(a * b)
        for lw, a, b in zip(liquid_gamma_lw_mN_m, liquid_gamma_acid_mN_m, liquid_gamma_base_mN_m)
    ]
    w = [work_of_adhesion(gt, theta) for gt, theta in zip(gamma_L_total, contact_angles_deg)]

    matrix = [
        [math.sqrt(liquid_gamma_lw_mN_m[i]), math.sqrt(liquid_gamma_base_mN_m[i]), math.sqrt(liquid_gamma_acid_mN_m[i])]
        for i in range(3)
    ]
    vector = [wi / 2.0 for wi in w]
    x, y, z = _solve_3x3(matrix, vector)

    gamma_s_lw = x ** 2
    gamma_s_acid = y ** 2
    gamma_s_base = z ** 2

    return VanOssChaudhuryGoodResult(
        gamma_s_lw_mN_m=gamma_s_lw,
        gamma_s_acid_mN_m=gamma_s_acid,
        gamma_s_base_mN_m=gamma_s_base,
        gamma_s_total_mN_m=gamma_s_lw + 2.0 * math.sqrt(gamma_s_acid * gamma_s_base),
        n_liquids=3,
    )


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
