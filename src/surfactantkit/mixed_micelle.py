"""Clint ideal mixing and Rubingh regular-solution theory for binary
surfactant mixtures.

Ported and generalized from a calculator originally built for a specific
cationic/anionic (amidoamine gemini / SLG) binary system, validated
against real conductometry and tensiometry data. Component 1 and
component 2 are generic here -- there is nothing cationic-specific left
in the math.
"""

from __future__ import annotations

import math

R_GAS = 8.314462618  # J / (mol K)


def clint_ideal_cmc(alpha1: float, cmc1: float, cmc2: float) -> float:
    """Clint's ideal mixed CMC for a binary surfactant mixture.

    alpha1: bulk mole fraction of component 1 in the total surfactant.
    cmc1, cmc2: pure-component CMCs (same concentration unit; the
    returned value is in that same unit).
    """
    if not (0.0 < alpha1 < 1.0):
        raise ValueError("alpha1 must be strictly between 0 and 1")
    if cmc1 <= 0 or cmc2 <= 0:
        raise ValueError("pure-component CMCs must be positive")
    return 1.0 / ((alpha1 / cmc1) + ((1.0 - alpha1) / cmc2))


def _rubingh_residual(x1: float, alpha1: float, cmc_mix: float, cmc1: float, cmc2: float) -> float:
    term1 = math.log((alpha1 * cmc_mix) / (x1 * cmc1))
    term2 = math.log(((1.0 - alpha1) * cmc_mix) / ((1.0 - x1) * cmc2))
    return (x1 * x1 * term1) - (((1.0 - x1) ** 2) * term2)


def solve_rubingh_x(
    alpha1: float,
    cmc_mix: float,
    cmc1: float,
    cmc2: float,
    grid_points: int = 4000,
    bisection_iters: int = 100,
) -> float | None:
    """Solve Rubingh's regular-solution equation for x1, the micellar mole
    fraction of component 1.

    Scans a fine grid for a sign change in the residual, then bisects to
    high precision. Returns None if no root is found in (0, 1) -- this
    happens for genuinely ideal or near-ideal mixtures, or for input data
    that isn't self-consistent (e.g. cmc_mix outside the physically
    possible range for the given alpha1).
    """
    eps = 1e-6
    xs = [eps + i * (1.0 - 2.0 * eps) / grid_points for i in range(grid_points + 1)]
    values: list[tuple[float, float]] = []
    for x in xs:
        try:
            fx = _rubingh_residual(x, alpha1, cmc_mix, cmc1, cmc2)
            if math.isfinite(fx):
                values.append((x, fx))
        except (ValueError, ZeroDivisionError):
            continue
    if not values:
        return None

    for (x_lo, f_lo), (x_hi, f_hi) in zip(values, values[1:]):
        if f_lo == 0:
            return x_lo
        if f_lo * f_hi < 0:
            lo, hi = x_lo, x_hi
            flo = f_lo
            for _ in range(bisection_iters):
                mid = (lo + hi) / 2.0
                fmid = _rubingh_residual(mid, alpha1, cmc_mix, cmc1, cmc2)
                if abs(fmid) < 1e-14:
                    return mid
                if flo * fmid <= 0:
                    hi = mid
                else:
                    lo, flo = mid, fmid
            return (lo + hi) / 2.0

    best_x, best_val = min(values, key=lambda pair: abs(pair[1]))
    return best_x if abs(best_val) < 1e-6 else None


def rubingh_beta(x1: float, alpha1: float, cmc_mix: float, cmc1: float) -> float:
    """Rubingh interaction parameter beta, given the solved micellar mole
    fraction x1 (see solve_rubingh_x)."""
    if x1 <= 0.0 or x1 >= 1.0:
        raise ValueError("x1 must be strictly between 0 and 1")
    term = math.log((alpha1 * cmc_mix) / (x1 * cmc1))
    return term / ((1.0 - x1) ** 2)


def activity_coefficients(x1: float, beta: float) -> tuple[float, float]:
    """Activity coefficients (f1, f2) from regular solution theory."""
    f1 = math.exp(beta * (1.0 - x1) ** 2)
    f2 = math.exp(beta * x1 ** 2)
    return f1, f2


def excess_free_energy(x1: float, f1: float, f2: float, temperature_k: float) -> float:
    """Excess free energy of micelle formation, kJ/mol."""
    return (R_GAS * temperature_k * (x1 * math.log(f1) + (1.0 - x1) * math.log(f2))) / 1000.0


def solve_rosen_monolayer_x(alpha1: float, c_mix_sigma: float, c1_sigma: float, c2_sigma: float, **kwargs) -> float | None:
    """Rosen's extension of Rubingh's regular-solution theory to the
    mixed ADSORBED MONOLAYER at an interface, from surface tension data
    rather than CMC data. Mathematically identical to solve_rubingh_x --
    the difference is entirely in what the inputs mean, and getting that
    right matters: c1_sigma/c2_sigma are the pure-component surfactant
    concentrations required to reach a chosen REFERENCE surface tension
    (or surface pressure) at the interface, and c_mix_sigma is the total
    mixed-surfactant concentration required to reach that SAME reference
    surface tension at bulk mole fraction alpha1 -- these are NOT CMC
    values. Using CMC values here by mistake silently produces a
    micellization-monolayer confusion, not an error."""
    return solve_rubingh_x(alpha1, c_mix_sigma, c1_sigma, c2_sigma, **kwargs)


def rosen_beta_sigma(x1_sigma: float, alpha1: float, c_mix_sigma: float, c1_sigma: float) -> float:
    """Rosen's monolayer interaction parameter beta^sigma, given the
    solved monolayer mole fraction x1_sigma (see solve_rosen_monolayer_x).
    Same formula as rubingh_beta -- see that function's docstring for
    the input-meaning distinction (surface-tension-derived concentrations,
    not CMC values)."""
    return rubingh_beta(x1_sigma, alpha1, c_mix_sigma, c1_sigma)


def corrin_harkins_predict_cmc(cmc1: float, c_counterion1: float, cmc2: float, c_counterion2: float, c_counterion_target: float) -> tuple[float, float]:
    """Predict CMC at a new counterion (salt) concentration using the
    Corrin-Harkins log-linear relation: log10(CMC) = -g*log10(C) + b.
    Fits g and b from two known (CMC, counterion concentration) points,
    then predicts CMC at c_counterion_target. Returns (predicted_cmc, g).

    g and b are SYSTEM-SPECIFIC (surfactant + salt identity), not
    universal constants -- this is why two real data points are required
    as input rather than a lookup table. All CMC and concentration
    values must be in the same, consistent unit (e.g. all mM); the
    prediction is only reliable as an interpolation between (or close
    extrapolation beyond) the two given concentrations.
    """
    if cmc1 <= 0 or cmc2 <= 0 or c_counterion1 <= 0 or c_counterion2 <= 0 or c_counterion_target <= 0:
        raise ValueError("all CMC and concentration values must be positive")
    if c_counterion1 == c_counterion2:
        raise ValueError("c_counterion1 and c_counterion2 must differ to fit a slope")
    log_cmc1, log_cmc2 = math.log10(cmc1), math.log10(cmc2)
    log_c1, log_c2 = math.log10(c_counterion1), math.log10(c_counterion2)
    g = (log_cmc2 - log_cmc1) / (log_c1 - log_c2)
    b = log_cmc1 + g * log_c1
    log_cmc_target = -g * math.log10(c_counterion_target) + b
    return 10.0 ** log_cmc_target, g
