"""AUMC (Area Under the Moment Curve) calculations."""

from typing import Optional
import numpy as np

from pynca.utils.validation import (
    validate_concentration_data,
    validate_time_data,
    validate_conc_time_match,
)
from pynca.utils.helpers import find_measurable_range


def _linear_moment_trapezoid(
    c1: float, c2: float, t1: float, t2: float
) -> float:
    """Calculate moment area using linear trapezoidal rule."""
    return (c1 * t1 + c2 * t2) * (t2 - t1) / 2


def _log_moment_trapezoid(
    c1: float, c2: float, t1: float, t2: float
) -> float:
    """
    Calculate moment area using log-linear trapezoidal rule.

    Uses linear method if concentrations are equal or one is zero.
    """
    if c1 <= 0 or c2 <= 0 or c1 == c2:
        return _linear_moment_trapezoid(c1, c2, t1, t2)

    # Log-linear moment calculation
    delta_t = t2 - t1
    log_ratio = np.log(c1) - np.log(c2)

    term1 = (c1 * t1 - c2 * t2) / log_ratio
    term2 = (c1 - c2) * delta_t / (log_ratio ** 2)

    return term1 * delta_t + term2


def calc_aumc(
    conc: np.ndarray,
    time: np.ndarray,
    method: str = "linear",
    start: Optional[float] = None,
    end: Optional[float] = None,
) -> float:
    """
    Calculate Area Under the Moment Curve (AUMC).

    AUMC = integral of (concentration * time) vs time

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    method : str
        Integration method: "linear", "log", or "linear-up/log-down"
    start : float, optional
        Start time (if None, use first time point)
    end : float, optional
        End time (if None, use last time point)

    Returns
    -------
    float
        AUMC value
    """
    conc = validate_concentration_data(conc)
    time = validate_time_data(time)
    validate_conc_time_match(conc, time)

    if len(conc) < 2:
        return 0.0

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    # Filter to interval if specified
    if start is not None or end is not None:
        start = start if start is not None else time[0]
        end = end if end is not None else time[-1]
        mask = (time >= start) & (time <= end)
        conc = conc[mask]
        time = time[mask]

    if len(conc) < 2:
        return 0.0

    # Replace NaN with 0
    conc = np.nan_to_num(conc, nan=0.0)

    aumc = 0.0

    for i in range(len(conc) - 1):
        c1, c2 = conc[i], conc[i + 1]
        t1, t2 = time[i], time[i + 1]

        if method == "linear":
            aumc += _linear_moment_trapezoid(c1, c2, t1, t2)
        elif method == "log":
            aumc += _log_moment_trapezoid(c1, c2, t1, t2)
        elif method == "linear-up/log-down":
            if c2 >= c1:
                aumc += _linear_moment_trapezoid(c1, c2, t1, t2)
            else:
                aumc += _log_moment_trapezoid(c1, c2, t1, t2)
        else:
            raise ValueError(f"Unknown AUMC method: {method}")

    return aumc


def calc_aumc_last(
    conc: np.ndarray,
    time: np.ndarray,
    method: str = "linear",
    loq: Optional[float] = None,
) -> float:
    """
    Calculate AUMC from time 0 to last measurable concentration.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    method : str
        Integration method
    loq : float, optional
        Limit of quantification

    Returns
    -------
    float
        AUMClast value
    """
    conc = validate_concentration_data(conc)
    time = validate_time_data(time)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    first_idx, last_idx = find_measurable_range(conc, time, loq)

    if first_idx is None or last_idx is None:
        return np.nan

    return calc_aumc(
        conc[: last_idx + 1],
        time[: last_idx + 1],
        method=method,
    )


def calc_aumc_inf(
    conc: np.ndarray,
    time: np.ndarray,
    lambda_z: float,
    method: str = "linear",
    loq: Optional[float] = None,
    clast_obs: Optional[float] = None,
    tlast: Optional[float] = None,
) -> float:
    """
    Calculate AUMC extrapolated to infinity.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    lambda_z : float
        Terminal elimination rate constant
    method : str
        Integration method
    loq : float, optional
        Limit of quantification
    clast_obs : float, optional
        Observed Clast
    tlast : float, optional
        Time of Clast

    Returns
    -------
    float
        AUMCinf value
    """
    if lambda_z is None or np.isnan(lambda_z) or lambda_z <= 0:
        return np.nan

    conc = validate_concentration_data(conc)
    time = validate_time_data(time)

    aumc_last = calc_aumc_last(conc, time, method=method, loq=loq)

    if np.isnan(aumc_last):
        return np.nan

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    # Get Clast and Tlast
    if clast_obs is None or tlast is None:
        _, last_idx = find_measurable_range(conc, time, loq)
        if last_idx is None:
            return np.nan
        clast_obs = conc[last_idx]
        tlast = time[last_idx]

    # Extrapolation: AUMCinf = AUMClast + (Clast*Tlast)/lambda_z + Clast/lambda_z^2
    term1 = (clast_obs * tlast) / lambda_z
    term2 = clast_obs / (lambda_z ** 2)

    return aumc_last + term1 + term2


def calc_aumc_pct_extrap(
    aumc_last: float,
    aumc_inf: float,
) -> float:
    """
    Calculate percent of AUMC extrapolated.

    Parameters
    ----------
    aumc_last : float
        AUMC to last measurable concentration
    aumc_inf : float
        AUMC extrapolated to infinity

    Returns
    -------
    float
        Percent extrapolated
    """
    if np.isnan(aumc_last) or np.isnan(aumc_inf) or aumc_inf == 0:
        return np.nan

    aumc_extrap = aumc_inf - aumc_last
    return (aumc_extrap / aumc_inf) * 100
