"""Helper utility functions."""

from typing import Tuple, Optional, Union
import numpy as np
import pandas as pd


def is_blq(
    conc: Union[float, np.ndarray],
    loq: Optional[float] = None,
) -> Union[bool, np.ndarray]:
    """
    Check if concentration value(s) are below limit of quantification.

    Parameters
    ----------
    conc : float or array-like
        Concentration value(s)
    loq : float, optional
        Limit of quantification. If None, checks for zero/NaN

    Returns
    -------
    bool or np.ndarray
        True if BLQ, False otherwise
    """
    conc = np.asarray(conc)

    if loq is not None:
        return (conc < loq) | np.isnan(conc)
    else:
        return (conc <= 0) | np.isnan(conc)


def find_measurable_range(
    conc: np.ndarray,
    time: np.ndarray,
    loq: Optional[float] = None,
) -> Tuple[int, int]:
    """
    Find the index range of measurable (non-BLQ) concentrations.

    Parameters
    ----------
    conc : np.ndarray
        Concentration values
    time : np.ndarray
        Time values
    loq : float, optional
        Limit of quantification

    Returns
    -------
    tuple
        (first_idx, last_idx) of measurable concentration range
    """
    blq_mask = is_blq(conc, loq)
    measurable_idx = np.where(~blq_mask)[0]

    if len(measurable_idx) == 0:
        return (None, None)

    return (measurable_idx[0], measurable_idx[-1])


def sort_by_time(
    conc: np.ndarray,
    time: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Sort concentration and time arrays by time.

    Parameters
    ----------
    conc : np.ndarray
        Concentration values
    time : np.ndarray
        Time values

    Returns
    -------
    tuple
        (sorted_conc, sorted_time)
    """
    sort_idx = np.argsort(time)
    return conc[sort_idx], time[sort_idx]


def get_tlast(
    conc: np.ndarray,
    time: np.ndarray,
    loq: Optional[float] = None,
) -> Optional[float]:
    """
    Get time of last measurable concentration.

    Parameters
    ----------
    conc : np.ndarray
        Concentration values
    time : np.ndarray
        Time values
    loq : float, optional
        Limit of quantification

    Returns
    -------
    float or None
        Time of last measurable concentration
    """
    _, last_idx = find_measurable_range(conc, time, loq)
    if last_idx is None:
        return None
    return time[last_idx]


def get_clast(
    conc: np.ndarray,
    time: np.ndarray,
    loq: Optional[float] = None,
) -> Optional[float]:
    """
    Get last measurable concentration.

    Parameters
    ----------
    conc : np.ndarray
        Concentration values
    time : np.ndarray
        Time values
    loq : float, optional
        Limit of quantification

    Returns
    -------
    float or None
        Last measurable concentration
    """
    _, last_idx = find_measurable_range(conc, time, loq)
    if last_idx is None:
        return None
    return conc[last_idx]


def filter_to_interval(
    conc: np.ndarray,
    time: np.ndarray,
    start: float,
    end: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Filter concentration-time data to a specific interval.

    Parameters
    ----------
    conc : np.ndarray
        Concentration values
    time : np.ndarray
        Time values
    start : float
        Interval start time
    end : float
        Interval end time

    Returns
    -------
    tuple
        (filtered_conc, filtered_time)
    """
    mask = (time >= start) & (time <= end)
    return conc[mask], time[mask]


def log_safe(
    values: np.ndarray,
    replace_zero: bool = True,
    min_value: float = 1e-10,
) -> np.ndarray:
    """
    Safely compute logarithm, handling zeros and negatives.

    Parameters
    ----------
    values : np.ndarray
        Values to take log of
    replace_zero : bool
        Whether to replace zeros with min_value
    min_value : float
        Minimum value to use for replacement

    Returns
    -------
    np.ndarray
        Log-transformed values
    """
    values = np.asarray(values, dtype=float)
    result = values.copy()

    if replace_zero:
        result[result <= 0] = min_value
    else:
        result[result <= 0] = np.nan

    return np.log(result)
