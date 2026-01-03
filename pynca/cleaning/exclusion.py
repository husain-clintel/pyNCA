"""Data exclusion and outlier detection functions."""

from typing import List, Optional, Tuple, Union
import numpy as np
from scipy import stats


def exclude_points(
    conc: np.ndarray,
    time: np.ndarray,
    indices: Optional[List[int]] = None,
    times: Optional[List[float]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Exclude specific data points.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    indices : list of int, optional
        Indices to exclude
    times : list of float, optional
        Time points to exclude

    Returns
    -------
    tuple
        (filtered_conc, filtered_time)
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)

    keep_mask = np.ones(len(conc), dtype=bool)

    if indices is not None:
        for idx in indices:
            if 0 <= idx < len(conc):
                keep_mask[idx] = False

    if times is not None:
        for t in times:
            time_mask = np.isclose(time, t)
            keep_mask[time_mask] = False

    return conc[keep_mask], time[keep_mask]


def flag_outliers(
    conc: np.ndarray,
    time: np.ndarray,
    method: str = "mad",
    threshold: float = 3.0,
) -> np.ndarray:
    """
    Flag potential outliers in concentration data.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    method : str
        Outlier detection method:
        - "mad": Median Absolute Deviation
        - "zscore": Z-score
        - "iqr": Interquartile Range
    threshold : float
        Threshold for flagging (default 3.0)

    Returns
    -------
    np.ndarray
        Boolean array (True = outlier)
    """
    conc = np.asarray(conc, dtype=float)

    # Remove NaN for calculation
    valid_mask = ~np.isnan(conc)
    valid_conc = conc[valid_mask]

    if len(valid_conc) < 3:
        return np.zeros(len(conc), dtype=bool)

    outlier_flags = np.zeros(len(conc), dtype=bool)

    if method == "mad":
        median = np.median(valid_conc)
        mad = np.median(np.abs(valid_conc - median))
        if mad == 0:
            return outlier_flags
        modified_z = 0.6745 * (conc - median) / mad
        outlier_flags = np.abs(modified_z) > threshold

    elif method == "zscore":
        mean = np.mean(valid_conc)
        std = np.std(valid_conc)
        if std == 0:
            return outlier_flags
        z_scores = (conc - mean) / std
        outlier_flags = np.abs(z_scores) > threshold

    elif method == "iqr":
        q1 = np.percentile(valid_conc, 25)
        q3 = np.percentile(valid_conc, 75)
        iqr = q3 - q1
        lower = q1 - threshold * iqr
        upper = q3 + threshold * iqr
        outlier_flags = (conc < lower) | (conc > upper)

    else:
        raise ValueError(f"Unknown outlier method: {method}")

    # Don't flag NaN values
    outlier_flags[~valid_mask] = False

    return outlier_flags


def flag_terminal_outliers(
    conc: np.ndarray,
    time: np.ndarray,
    residual_threshold: float = 2.0,
) -> np.ndarray:
    """
    Flag outliers in terminal phase based on regression residuals.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    residual_threshold : float
        Standardized residual threshold for flagging

    Returns
    -------
    np.ndarray
        Boolean array (True = outlier)
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)

    outlier_flags = np.zeros(len(conc), dtype=bool)

    # Find terminal phase (after Cmax)
    valid_mask = (conc > 0) & ~np.isnan(conc)
    if np.sum(valid_mask) < 3:
        return outlier_flags

    cmax_idx = np.nanargmax(conc)

    # Get terminal phase data
    terminal_mask = valid_mask.copy()
    terminal_mask[:cmax_idx] = False

    terminal_indices = np.where(terminal_mask)[0]
    if len(terminal_indices) < 3:
        return outlier_flags

    # Fit log-linear regression
    t_term = time[terminal_mask]
    log_c_term = np.log(conc[terminal_mask])

    slope, intercept, _, _, _ = stats.linregress(t_term, log_c_term)

    # Calculate standardized residuals
    predicted = intercept + slope * t_term
    residuals = log_c_term - predicted
    std_residuals = (residuals - np.mean(residuals)) / np.std(residuals)

    # Flag outliers
    term_outliers = np.abs(std_residuals) > residual_threshold

    for i, idx in enumerate(terminal_indices):
        outlier_flags[idx] = term_outliers[i]

    return outlier_flags


def exclude_outliers(
    conc: np.ndarray,
    time: np.ndarray,
    method: str = "mad",
    threshold: float = 3.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Detect and exclude outliers.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    method : str
        Outlier detection method
    threshold : float
        Detection threshold

    Returns
    -------
    tuple
        (filtered_conc, filtered_time, outlier_mask)
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)

    outlier_mask = flag_outliers(conc, time, method=method, threshold=threshold)
    keep_mask = ~outlier_mask

    return conc[keep_mask], time[keep_mask], outlier_mask
