"""Data imputation functions."""

from typing import Optional, Tuple
import numpy as np


def clean_na(
    conc: np.ndarray,
    time: np.ndarray,
    method: str = "drop",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Handle missing (NA/NaN) concentration values.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    method : str
        Handling method:
        - "drop": Remove rows with missing values
        - "zero": Replace missing with 0

    Returns
    -------
    tuple
        (cleaned_conc, cleaned_time)
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)

    na_mask = np.isnan(conc)

    if method == "drop":
        keep_mask = ~na_mask
        return conc[keep_mask], time[keep_mask]

    elif method == "zero":
        result_conc = conc.copy()
        result_conc[na_mask] = 0.0
        return result_conc, time

    else:
        raise ValueError(f"Unknown method: {method}")


def impute_conc(
    conc: np.ndarray,
    time: np.ndarray,
    target_time: float,
    method: str = "linear",
    lambda_z: Optional[float] = None,
) -> float:
    """
    Impute concentration at a specific time point.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    target_time : float
        Time point to impute
    method : str
        Imputation method:
        - "linear": Linear interpolation
        - "log-linear": Log-linear interpolation
        - "zero": Return 0
    lambda_z : float, optional
        Terminal rate constant for extrapolation

    Returns
    -------
    float
        Imputed concentration
    """
    from pynca.interpolation.interpolate import interpolate_conc
    from pynca.interpolation.extrapolate import extrapolate_conc

    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    if len(conc) == 0:
        return np.nan

    # Check if target is before first time
    if target_time < time[0]:
        return 0.0 if method == "zero" else np.nan

    # Check if target is after last time
    if target_time > time[-1]:
        if lambda_z is not None and lambda_z > 0:
            clast = conc[-1]
            tlast = time[-1]
            return extrapolate_conc(clast, tlast, target_time, lambda_z)
        return np.nan

    # Interpolate
    return interpolate_conc(conc, time, target_time, method=method)


def impute_predose(
    conc: np.ndarray,
    time: np.ndarray,
    dose_time: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Add predose (time 0) concentration if missing.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    dose_time : float
        Time of dose (default 0)

    Returns
    -------
    tuple
        (conc_with_predose, time_with_predose)
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    # Check if predose already exists
    if dose_time in time:
        return conc, time

    # Add predose value of 0
    new_conc = np.concatenate([[0.0], conc])
    new_time = np.concatenate([[dose_time], time])

    return new_conc, new_time


def impute_at_times(
    conc: np.ndarray,
    time: np.ndarray,
    target_times: np.ndarray,
    method: str = "linear",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Impute concentrations at specified time points.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    target_times : array-like
        Time points to add
    method : str
        Imputation method

    Returns
    -------
    tuple
        (conc_with_imputed, time_with_imputed)
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)
    target_times = np.asarray(target_times, dtype=float)

    result_conc = list(conc)
    result_time = list(time)

    for t in target_times:
        if t not in time:
            c = impute_conc(conc, time, t, method=method)
            result_conc.append(c)
            result_time.append(t)

    # Sort by time
    result_conc = np.array(result_conc)
    result_time = np.array(result_time)
    sort_idx = np.argsort(result_time)

    return result_conc[sort_idx], result_time[sort_idx]
