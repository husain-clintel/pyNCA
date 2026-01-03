"""AUC (Area Under the Curve) calculations."""

from typing import Optional, Tuple, Union, Literal
import numpy as np

from pynca.utils.validation import (
    validate_concentration_data,
    validate_time_data,
    validate_conc_time_match,
    validate_interval,
)
from pynca.utils.helpers import find_measurable_range, log_safe


def _linear_trapezoid(c1: float, c2: float, t1: float, t2: float) -> float:
    """Calculate area of linear trapezoid."""
    return (c1 + c2) * (t2 - t1) / 2


def _log_trapezoid(c1: float, c2: float, t1: float, t2: float) -> float:
    """
    Calculate area of log-linear trapezoid.

    Uses linear trapezoidal if concentrations are equal or one is zero.
    """
    if c1 <= 0 or c2 <= 0 or c1 == c2:
        return _linear_trapezoid(c1, c2, t1, t2)

    return (c1 - c2) * (t2 - t1) / (np.log(c1) - np.log(c2))


def calc_auc(
    conc: np.ndarray,
    time: np.ndarray,
    method: str = "linear",
    start: Optional[float] = None,
    end: Optional[float] = None,
) -> float:
    """
    Calculate Area Under the Curve (AUC).

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
        AUC value

    Notes
    -----
    - "linear": Linear trapezoidal throughout
    - "log": Log-linear trapezoidal throughout
    - "linear-up/log-down": Linear trapezoidal for ascending, log for descending
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

    auc = 0.0

    for i in range(len(conc) - 1):
        c1, c2 = conc[i], conc[i + 1]
        t1, t2 = time[i], time[i + 1]

        if method == "linear":
            auc += _linear_trapezoid(c1, c2, t1, t2)
        elif method == "log":
            auc += _log_trapezoid(c1, c2, t1, t2)
        elif method == "linear-up/log-down":
            if c2 >= c1:
                # Ascending or plateau - use linear
                auc += _linear_trapezoid(c1, c2, t1, t2)
            else:
                # Descending - use log-linear
                auc += _log_trapezoid(c1, c2, t1, t2)
        else:
            raise ValueError(f"Unknown AUC method: {method}")

    return auc


def calc_auc_last(
    conc: np.ndarray,
    time: np.ndarray,
    method: str = "linear",
    loq: Optional[float] = None,
) -> float:
    """
    Calculate AUC from time 0 to the last measurable concentration (AUClast).

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
        AUClast value
    """
    conc = validate_concentration_data(conc)
    time = validate_time_data(time)

    first_idx, last_idx = find_measurable_range(conc, time, loq)

    if first_idx is None or last_idx is None:
        return np.nan

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    # Find last measurable after sorting
    first_idx, last_idx = find_measurable_range(conc, time, loq)

    # Calculate from first time to last measurable
    return calc_auc(
        conc[: last_idx + 1],
        time[: last_idx + 1],
        method=method,
    )


def calc_auc_all(
    conc: np.ndarray,
    time: np.ndarray,
    method: str = "linear",
) -> float:
    """
    Calculate AUC including all time points (AUCall).

    This includes trailing BLQ values in the calculation.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    method : str
        Integration method

    Returns
    -------
    float
        AUCall value
    """
    return calc_auc(conc, time, method=method)


def calc_auc_inf(
    conc: np.ndarray,
    time: np.ndarray,
    lambda_z: float,
    method: str = "linear",
    loq: Optional[float] = None,
    clast_obs: Optional[float] = None,
) -> float:
    """
    Calculate AUC extrapolated to infinity (AUCinf).

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    lambda_z : float
        Terminal elimination rate constant
    method : str
        Integration method for AUClast
    loq : float, optional
        Limit of quantification
    clast_obs : float, optional
        Observed Clast (if None, calculated from data)

    Returns
    -------
    float
        AUCinf value
    """
    if lambda_z is None or np.isnan(lambda_z) or lambda_z <= 0:
        return np.nan

    conc = validate_concentration_data(conc)
    time = validate_time_data(time)

    auc_last = calc_auc_last(conc, time, method=method, loq=loq)

    if np.isnan(auc_last):
        return np.nan

    # Get Clast
    if clast_obs is None:
        _, last_idx = find_measurable_range(conc, time, loq)
        if last_idx is None:
            return np.nan
        sort_idx = np.argsort(time)
        conc = conc[sort_idx]
        _, last_idx = find_measurable_range(conc, time, loq)
        clast_obs = conc[last_idx]

    # Extrapolate: AUCinf = AUClast + Clast/lambda_z
    auc_extrap = clast_obs / lambda_z
    return auc_last + auc_extrap


def calc_auc_int(
    conc: np.ndarray,
    time: np.ndarray,
    start: float,
    end: float,
    method: str = "linear",
    lambda_z: Optional[float] = None,
    clast: Optional[float] = None,
    tlast: Optional[float] = None,
) -> float:
    """
    Calculate AUC over a specific interval with interpolation/extrapolation.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    start : float
        Interval start time
    end : float
        Interval end time
    method : str
        Integration method
    lambda_z : float, optional
        Terminal rate constant for extrapolation
    clast : float, optional
        Last measurable concentration
    tlast : float, optional
        Time of last measurable concentration

    Returns
    -------
    float
        AUC over interval
    """
    from pynca.interpolation.interpolate import interpolate_conc
    from pynca.interpolation.extrapolate import extrapolate_conc

    conc = validate_concentration_data(conc)
    time = validate_time_data(time)
    validate_interval(start, end)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    # Build concentration-time data for the interval
    new_conc = []
    new_time = []

    # Handle start point
    if start < time[0]:
        # Before first observation - assume 0
        new_conc.append(0.0)
        new_time.append(start)
    elif start in time:
        idx = np.where(time == start)[0][0]
        new_conc.append(conc[idx])
        new_time.append(start)
    else:
        # Interpolate
        c_start = interpolate_conc(conc, time, start, method=method)
        new_conc.append(c_start)
        new_time.append(start)

    # Add points within interval
    mask = (time > start) & (time < end)
    new_conc.extend(conc[mask].tolist())
    new_time.extend(time[mask].tolist())

    # Handle end point
    if end <= time[-1]:
        if end in time:
            idx = np.where(time == end)[0][0]
            new_conc.append(conc[idx])
            new_time.append(end)
        else:
            # Interpolate
            c_end = interpolate_conc(conc, time, end, method=method)
            new_conc.append(c_end)
            new_time.append(end)
    else:
        # Extrapolate
        if lambda_z is not None and clast is not None and tlast is not None:
            c_end = extrapolate_conc(clast, tlast, end, lambda_z)
            new_conc.append(c_end)
            new_time.append(end)
        else:
            # Can't extrapolate - use last observation
            new_conc.append(conc[-1])
            new_time.append(end)

    return calc_auc(np.array(new_conc), np.array(new_time), method=method)


def calc_auc_pct_extrap(
    auc_last: float,
    auc_inf: float,
) -> float:
    """
    Calculate percent of AUC extrapolated.

    Parameters
    ----------
    auc_last : float
        AUC to last measurable concentration
    auc_inf : float
        AUC extrapolated to infinity

    Returns
    -------
    float
        Percent extrapolated
    """
    if np.isnan(auc_last) or np.isnan(auc_inf) or auc_inf == 0:
        return np.nan

    auc_extrap = auc_inf - auc_last
    return (auc_extrap / auc_inf) * 100


def calc_auc_partial(
    conc: np.ndarray,
    time: np.ndarray,
    start: float,
    end: float,
    method: str = "linear",
    lambda_z: Optional[float] = None,
) -> float:
    """
    Calculate partial AUC between two arbitrary time points.

    This function handles interpolation at boundaries and extrapolation
    if the end time exceeds the last observation.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    start : float
        Start time for partial AUC
    end : float
        End time for partial AUC
    method : str
        Integration method: "linear", "log", or "linear-up/log-down"
    lambda_z : float, optional
        Terminal elimination rate constant for extrapolation beyond tlast

    Returns
    -------
    float
        Partial AUC value

    Examples
    --------
    >>> conc = np.array([0, 10, 8, 6, 4, 2])
    >>> time = np.array([0, 1, 2, 4, 6, 8])
    >>> calc_auc_partial(conc, time, start=2, end=6)
    20.0
    """
    from pynca.interpolation.interpolate import interpolate_conc
    from pynca.interpolation.extrapolate import extrapolate_conc

    conc = validate_concentration_data(conc)
    time = validate_time_data(time)
    validate_interval(start, end)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    # Find tlast and clast for potential extrapolation
    first_idx, last_idx = find_measurable_range(conc, time)
    if last_idx is not None:
        tlast = time[last_idx]
        clast = conc[last_idx]
    else:
        tlast = time[-1]
        clast = conc[-1]

    # Build new time-concentration arrays for the interval
    new_time = []
    new_conc = []

    # Handle start boundary
    if start < time[0]:
        new_time.append(start)
        new_conc.append(0.0)
    elif start > time[-1]:
        # Start is beyond all data
        if lambda_z is not None and lambda_z > 0:
            c_start = extrapolate_conc(clast, tlast, start, lambda_z)
            new_time.append(start)
            new_conc.append(c_start)
        else:
            return np.nan
    elif start in time:
        idx = np.where(time == start)[0][0]
        new_time.append(start)
        new_conc.append(conc[idx])
    else:
        # Interpolate at start
        c_start = interpolate_conc(conc, time, start, method=method)
        new_time.append(start)
        new_conc.append(c_start)

    # Add all points within interval
    mask = (time > start) & (time < end)
    for t, c in zip(time[mask], conc[mask]):
        new_time.append(t)
        new_conc.append(c)

    # Handle end boundary
    if end <= time[-1]:
        if end in time:
            idx = np.where(time == end)[0][0]
            new_time.append(end)
            new_conc.append(conc[idx])
        else:
            # Interpolate at end
            c_end = interpolate_conc(conc, time, end, method=method)
            new_time.append(end)
            new_conc.append(c_end)
    else:
        # Extrapolate beyond last observation
        if lambda_z is not None and lambda_z > 0:
            c_end = extrapolate_conc(clast, tlast, end, lambda_z)
            new_time.append(end)
            new_conc.append(c_end)
        else:
            # Cannot extrapolate - use tlast as end
            if tlast > start:
                # Recalculate up to tlast only
                return calc_auc_partial(conc, time, start, tlast, method=method)
            else:
                return np.nan

    if len(new_time) < 2:
        return 0.0

    return calc_auc(np.array(new_conc), np.array(new_time), method=method)


def calc_auc_dn(
    auc: float,
    dose: float,
) -> float:
    """
    Calculate dose-normalized AUC.

    Parameters
    ----------
    auc : float
        AUC value (any type: AUClast, AUCinf, etc.)
    dose : float
        Dose amount

    Returns
    -------
    float
        Dose-normalized AUC (AUC/Dose)
    """
    if np.isnan(auc) or np.isnan(dose) or dose <= 0:
        return np.nan
    return auc / dose


def calc_cmax_dn(
    cmax: float,
    dose: float,
) -> float:
    """
    Calculate dose-normalized Cmax.

    Parameters
    ----------
    cmax : float
        Maximum concentration
    dose : float
        Dose amount

    Returns
    -------
    float
        Dose-normalized Cmax (Cmax/Dose)
    """
    if np.isnan(cmax) or np.isnan(dose) or dose <= 0:
        return np.nan
    return cmax / dose
