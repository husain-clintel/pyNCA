"""Concentration interpolation functions."""

from typing import Optional, Tuple
import numpy as np


def interpolate_conc(
    conc: np.ndarray,
    time: np.ndarray,
    target_time: float,
    method: str = "linear",
) -> float:
    """
    Interpolate concentration at a target time point.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    target_time : float
        Time point to interpolate
    method : str
        Interpolation method:
        - "linear": Linear interpolation
        - "log": Log-linear interpolation
        - "linear-up/log-down": Mixed method

    Returns
    -------
    float
        Interpolated concentration
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    if len(conc) == 0:
        return np.nan

    # Check bounds
    if target_time < time[0]:
        return np.nan
    if target_time > time[-1]:
        return np.nan

    # Check if target is an existing time point
    if target_time in time:
        idx = np.where(time == target_time)[0][0]
        return conc[idx]

    # Find surrounding points
    idx_after = np.searchsorted(time, target_time)
    idx_before = idx_after - 1

    t1, t2 = time[idx_before], time[idx_after]
    c1, c2 = conc[idx_before], conc[idx_after]

    if method == "linear":
        return interpolate_conc_linear(c1, c2, t1, t2, target_time)
    elif method == "log":
        return interpolate_conc_log(c1, c2, t1, t2, target_time)
    elif method == "linear-up/log-down":
        if c2 >= c1:
            # Ascending - use linear
            return interpolate_conc_linear(c1, c2, t1, t2, target_time)
        else:
            # Descending - use log
            return interpolate_conc_log(c1, c2, t1, t2, target_time)
    else:
        raise ValueError(f"Unknown interpolation method: {method}")


def interpolate_conc_linear(
    c1: float,
    c2: float,
    t1: float,
    t2: float,
    target_time: float,
) -> float:
    """
    Linear interpolation between two concentration-time points.

    Parameters
    ----------
    c1 : float
        Concentration at t1
    c2 : float
        Concentration at t2
    t1 : float
        First time point
    t2 : float
        Second time point
    target_time : float
        Time to interpolate

    Returns
    -------
    float
        Interpolated concentration
    """
    if t2 == t1:
        return (c1 + c2) / 2

    fraction = (target_time - t1) / (t2 - t1)
    return c1 + fraction * (c2 - c1)


def interpolate_conc_log(
    c1: float,
    c2: float,
    t1: float,
    t2: float,
    target_time: float,
) -> float:
    """
    Log-linear interpolation between two concentration-time points.

    Parameters
    ----------
    c1 : float
        Concentration at t1
    c2 : float
        Concentration at t2
    t1 : float
        First time point
    t2 : float
        Second time point
    target_time : float
        Time to interpolate

    Returns
    -------
    float
        Interpolated concentration
    """
    # Handle edge cases
    if c1 <= 0 or c2 <= 0:
        return interpolate_conc_linear(c1, c2, t1, t2, target_time)

    if t2 == t1:
        return np.sqrt(c1 * c2)  # Geometric mean

    if c1 == c2:
        return c1

    # Log-linear interpolation
    log_c1 = np.log(c1)
    log_c2 = np.log(c2)

    fraction = (target_time - t1) / (t2 - t1)
    log_result = log_c1 + fraction * (log_c2 - log_c1)

    return np.exp(log_result)


def interpolate_at_times(
    conc: np.ndarray,
    time: np.ndarray,
    target_times: np.ndarray,
    method: str = "linear",
) -> np.ndarray:
    """
    Interpolate concentrations at multiple time points.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    target_times : array-like
        Time points to interpolate
    method : str
        Interpolation method

    Returns
    -------
    np.ndarray
        Interpolated concentrations
    """
    target_times = np.asarray(target_times, dtype=float)
    result = np.zeros(len(target_times))

    for i, t in enumerate(target_times):
        result[i] = interpolate_conc(conc, time, t, method=method)

    return result


def interpolate_dose_aware(
    conc: np.ndarray,
    time: np.ndarray,
    target_time: float,
    dose_times: np.ndarray,
    method: str = "linear",
) -> float:
    """
    Dose-aware concentration interpolation.

    Does not interpolate across dose times. If target_time is in a different
    dosing interval than adjacent points, interpolation is not performed.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    target_time : float
        Time point to interpolate
    dose_times : array-like
        Times of dose administration
    method : str
        Interpolation method

    Returns
    -------
    float
        Interpolated concentration or NaN if not possible
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)
    dose_times = np.asarray(dose_times, dtype=float)

    # Sort
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]
    dose_times = np.sort(dose_times)

    if len(conc) == 0:
        return np.nan

    # Check bounds
    if target_time < time[0] or target_time > time[-1]:
        return np.nan

    # Check if target is an existing point
    if target_time in time:
        idx = np.where(time == target_time)[0][0]
        return conc[idx]

    # Find surrounding points
    idx_after = np.searchsorted(time, target_time)
    idx_before = idx_after - 1

    t1, t2 = time[idx_before], time[idx_after]

    # Check if a dose occurs between the two points
    doses_between = dose_times[(dose_times > t1) & (dose_times < t2)]
    if len(doses_between) > 0:
        # Cannot interpolate across dose
        return np.nan

    # Safe to interpolate
    return interpolate_conc(conc, time, target_time, method=method)
