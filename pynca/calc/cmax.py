"""Peak concentration parameter calculations."""

from typing import Optional, Tuple, Union
import numpy as np

from pynca.utils.validation import (
    validate_concentration_data,
    validate_time_data,
    validate_conc_time_match,
)
from pynca.utils.helpers import find_measurable_range, is_blq


def calc_cmax(
    conc: np.ndarray,
    time: np.ndarray,
    check_tmax: bool = True,
) -> float:
    """
    Calculate maximum concentration (Cmax).

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    check_tmax : bool
        If True, verify Tmax is reasonable

    Returns
    -------
    float
        Cmax value
    """
    conc = validate_concentration_data(conc)
    time = validate_time_data(time)
    validate_conc_time_match(conc, time)

    if len(conc) == 0:
        return np.nan

    # Ignore NaN values
    valid_mask = ~np.isnan(conc)
    if not np.any(valid_mask):
        return np.nan

    return np.nanmax(conc)


def calc_tmax(
    conc: np.ndarray,
    time: np.ndarray,
    first: bool = True,
) -> float:
    """
    Calculate time of maximum concentration (Tmax).

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    first : bool
        If True and multiple maxima exist, return first occurrence

    Returns
    -------
    float
        Tmax value
    """
    conc = validate_concentration_data(conc)
    time = validate_time_data(time)
    validate_conc_time_match(conc, time)

    if len(conc) == 0:
        return np.nan

    # Find max value
    cmax = np.nanmax(conc)
    if np.isnan(cmax):
        return np.nan

    # Find indices where conc equals cmax
    max_indices = np.where(conc == cmax)[0]

    if first:
        # Return time of first maximum
        return time[max_indices[0]]
    else:
        # Return time of last maximum
        return time[max_indices[-1]]


def calc_cmax_tmax(
    conc: np.ndarray,
    time: np.ndarray,
    first: bool = True,
) -> Tuple[float, float]:
    """
    Calculate Cmax and Tmax together.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    first : bool
        If True, return first Tmax if multiple maxima exist

    Returns
    -------
    tuple
        (Cmax, Tmax)
    """
    cmax = calc_cmax(conc, time)
    tmax = calc_tmax(conc, time, first=first)
    return cmax, tmax


def calc_cmin(
    conc: np.ndarray,
    time: np.ndarray,
    exclude_zero: bool = True,
) -> float:
    """
    Calculate minimum concentration (Cmin).

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    exclude_zero : bool
        If True, exclude zero values from minimum calculation

    Returns
    -------
    float
        Cmin value
    """
    conc = validate_concentration_data(conc)
    time = validate_time_data(time)
    validate_conc_time_match(conc, time)

    if len(conc) == 0:
        return np.nan

    values = conc.copy()
    if exclude_zero:
        values = values[values > 0]

    if len(values) == 0:
        return np.nan

    return np.nanmin(values)


def calc_tlast(
    conc: np.ndarray,
    time: np.ndarray,
    loq: Optional[float] = None,
) -> float:
    """
    Calculate time of last measurable concentration (Tlast).

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    loq : float, optional
        Limit of quantification

    Returns
    -------
    float
        Tlast value
    """
    conc = validate_concentration_data(conc)
    time = validate_time_data(time)
    validate_conc_time_match(conc, time)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    _, last_idx = find_measurable_range(conc, time, loq)

    if last_idx is None:
        return np.nan

    return time[last_idx]


def calc_clast(
    conc: np.ndarray,
    time: np.ndarray,
    loq: Optional[float] = None,
) -> float:
    """
    Calculate last measurable concentration (Clast).

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    loq : float, optional
        Limit of quantification

    Returns
    -------
    float
        Clast value
    """
    conc = validate_concentration_data(conc)
    time = validate_time_data(time)
    validate_conc_time_match(conc, time)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    _, last_idx = find_measurable_range(conc, time, loq)

    if last_idx is None:
        return np.nan

    return conc[last_idx]


def calc_tfirst(
    conc: np.ndarray,
    time: np.ndarray,
    loq: Optional[float] = None,
) -> float:
    """
    Calculate time of first measurable concentration (Tfirst).

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    loq : float, optional
        Limit of quantification

    Returns
    -------
    float
        Tfirst value
    """
    conc = validate_concentration_data(conc)
    time = validate_time_data(time)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    first_idx, _ = find_measurable_range(conc, time, loq)

    if first_idx is None:
        return np.nan

    return time[first_idx]


def calc_cav(
    auc: float,
    tau: float,
) -> float:
    """
    Calculate average concentration over interval (Cav).

    Parameters
    ----------
    auc : float
        AUC over the interval
    tau : float
        Dosing interval

    Returns
    -------
    float
        Cav value
    """
    if np.isnan(auc) or tau <= 0:
        return np.nan

    return auc / tau


def calc_ctrough(
    conc: np.ndarray,
    time: np.ndarray,
    tau: Optional[float] = None,
) -> float:
    """
    Calculate trough concentration (Ctrough).

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    tau : float, optional
        Dosing interval (if None, uses last time point)

    Returns
    -------
    float
        Ctrough value
    """
    conc = validate_concentration_data(conc)
    time = validate_time_data(time)
    validate_conc_time_match(conc, time)

    if tau is not None:
        # Find concentration at tau (or closest time point)
        idx = np.argmin(np.abs(time - tau))
        return conc[idx]
    else:
        # Return last concentration
        sort_idx = np.argsort(time)
        return conc[sort_idx[-1]]


def calc_swing(
    cmax: float,
    cmin: float,
) -> float:
    """
    Calculate swing (fluctuation from Cmin).

    Swing = (Cmax - Cmin) / Cmin

    Parameters
    ----------
    cmax : float
        Maximum concentration
    cmin : float
        Minimum concentration

    Returns
    -------
    float
        Swing value
    """
    if np.isnan(cmax) or np.isnan(cmin) or cmin <= 0:
        return np.nan

    return (cmax - cmin) / cmin


def calc_ptf(
    cmax: float,
    cmin: float,
    cav: float,
) -> float:
    """
    Calculate peak-trough fluctuation (PTF).

    PTF = (Cmax - Cmin) / Cav * 100

    Parameters
    ----------
    cmax : float
        Maximum concentration
    cmin : float
        Minimum concentration
    cav : float
        Average concentration

    Returns
    -------
    float
        PTF as percentage
    """
    if np.isnan(cmax) or np.isnan(cmin) or np.isnan(cav) or cav <= 0:
        return np.nan

    return ((cmax - cmin) / cav) * 100


def calc_fluctuation(
    cmax: float,
    cmin: float,
) -> float:
    """
    Calculate fluctuation.

    Fluctuation = (Cmax - Cmin) / Cmax

    Parameters
    ----------
    cmax : float
        Maximum concentration
    cmin : float
        Minimum concentration

    Returns
    -------
    float
        Fluctuation value
    """
    if np.isnan(cmax) or np.isnan(cmin) or cmax <= 0:
        return np.nan

    return (cmax - cmin) / cmax
