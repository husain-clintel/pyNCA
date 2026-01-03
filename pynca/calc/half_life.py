"""Half-life and terminal phase calculations."""

from typing import Dict, Optional, Tuple, List
import numpy as np
from scipy import stats

from pynca.utils.validation import (
    validate_concentration_data,
    validate_time_data,
    validate_conc_time_match,
)
from pynca.utils.helpers import find_measurable_range, log_safe


def calc_lambda_z(
    conc: np.ndarray,
    time: np.ndarray,
    n_points: Optional[int] = None,
    min_points: int = 3,
    max_points: Optional[int] = None,
    selection_method: str = "best_fit",
    adj_r_squared_threshold: float = 0.0,
    tlast: Optional[float] = None,
) -> Dict[str, float]:
    """
    Calculate terminal elimination rate constant (lambda_z).

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    n_points : int, optional
        Fixed number of points to use (overrides automatic selection)
    min_points : int
        Minimum number of points for regression
    max_points : int, optional
        Maximum number of points to consider
    selection_method : str
        Method for automatic point selection:
        - "best_fit": maximize adjusted R-squared
        - "manual": require n_points to be specified
    adj_r_squared_threshold : float
        Minimum adjusted R-squared to accept
    tlast : float, optional
        Time of last measurable concentration

    Returns
    -------
    dict
        Dictionary containing:
        - lambda_z: terminal rate constant
        - half_life: terminal half-life
        - r_squared: R-squared of regression
        - adj_r_squared: adjusted R-squared
        - intercept: y-intercept of log-linear regression
        - n_points: number of points used
        - time_range: (first_time, last_time) used in regression
        - points_used: indices of points used
    """
    conc = validate_concentration_data(conc)
    time = validate_time_data(time)
    validate_conc_time_match(conc, time)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    # Find terminal phase data points
    _, last_idx = find_measurable_range(conc, time)
    if last_idx is None:
        return _empty_lambda_z_result()

    # Find Tmax to exclude ascending phase
    cmax_idx = np.nanargmax(conc)

    # Use only points after Tmax and up to last measurable
    terminal_mask = np.zeros(len(conc), dtype=bool)
    for i in range(cmax_idx, last_idx + 1):
        if conc[i] > 0 and not np.isnan(conc[i]):
            terminal_mask[i] = True

    terminal_indices = np.where(terminal_mask)[0]

    if len(terminal_indices) < min_points:
        return _empty_lambda_z_result()

    # Limit max points
    if max_points is not None and len(terminal_indices) > max_points:
        terminal_indices = terminal_indices[-max_points:]

    # Select points based on method
    if n_points is not None:
        # Use specified number of points from end
        if n_points > len(terminal_indices):
            n_points = len(terminal_indices)
        if n_points < min_points:
            return _empty_lambda_z_result()
        points_to_use = terminal_indices[-n_points:]
        result = _fit_lambda_z(conc, time, points_to_use)
    elif selection_method == "best_fit":
        # Try different numbers of points, find best fit
        result = _best_fit_lambda_z(
            conc, time, terminal_indices, min_points, adj_r_squared_threshold
        )
    else:
        # Default: use all terminal points
        points_to_use = terminal_indices
        if len(points_to_use) >= min_points:
            result = _fit_lambda_z(conc, time, points_to_use)
        else:
            result = _empty_lambda_z_result()

    return result


def _fit_lambda_z(
    conc: np.ndarray,
    time: np.ndarray,
    indices: np.ndarray,
) -> Dict[str, float]:
    """
    Fit lambda_z using linear regression on log-transformed concentrations.
    """
    t = time[indices]
    c = conc[indices]

    # Log transform
    log_c = np.log(c)

    # Linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(t, log_c)

    # lambda_z is negative of slope
    lambda_z = -slope

    # Calculate adjusted R-squared
    n = len(t)
    r_squared = r_value ** 2
    adj_r_squared = 1 - (1 - r_squared) * (n - 1) / (n - 2) if n > 2 else r_squared

    # Half-life
    half_life = np.log(2) / lambda_z if lambda_z > 0 else np.nan

    return {
        "lambda_z": lambda_z if lambda_z > 0 else np.nan,
        "half_life": half_life,
        "r_squared": r_squared,
        "adj_r_squared": adj_r_squared,
        "intercept": intercept,
        "n_points": n,
        "time_range": (t[0], t[-1]),
        "points_used": indices.tolist(),
        "slope": slope,
        "std_err": std_err,
    }


def _best_fit_lambda_z(
    conc: np.ndarray,
    time: np.ndarray,
    terminal_indices: np.ndarray,
    min_points: int,
    adj_r_squared_threshold: float,
) -> Dict[str, float]:
    """
    Find best lambda_z fit by maximizing adjusted R-squared.
    """
    best_result = _empty_lambda_z_result()
    best_adj_r2 = adj_r_squared_threshold

    n_terminal = len(terminal_indices)

    # Try different numbers of points, starting from last
    for n in range(min_points, n_terminal + 1):
        indices = terminal_indices[-n:]
        result = _fit_lambda_z(conc, time, indices)

        # Only accept if lambda_z is positive
        if result["lambda_z"] > 0 and result["adj_r_squared"] > best_adj_r2:
            best_result = result
            best_adj_r2 = result["adj_r_squared"]

    return best_result


def _empty_lambda_z_result() -> Dict[str, float]:
    """Return empty lambda_z result."""
    return {
        "lambda_z": np.nan,
        "half_life": np.nan,
        "r_squared": np.nan,
        "adj_r_squared": np.nan,
        "intercept": np.nan,
        "n_points": 0,
        "time_range": (np.nan, np.nan),
        "points_used": [],
        "slope": np.nan,
        "std_err": np.nan,
    }


def calc_half_life(
    conc: np.ndarray,
    time: np.ndarray,
    n_points: Optional[int] = None,
    min_points: int = 3,
    max_points: Optional[int] = None,
    selection_method: str = "best_fit",
    adj_r_squared_threshold: float = 0.0,
) -> Dict[str, float]:
    """
    Calculate terminal elimination half-life.

    This is a wrapper around calc_lambda_z that returns the same results.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    n_points : int, optional
        Fixed number of points to use
    min_points : int
        Minimum number of points for regression
    max_points : int, optional
        Maximum number of points to consider
    selection_method : str
        Point selection method: "best_fit" or "manual"
    adj_r_squared_threshold : float
        Minimum adjusted R-squared

    Returns
    -------
    dict
        Dictionary with half-life calculation results
    """
    return calc_lambda_z(
        conc=conc,
        time=time,
        n_points=n_points,
        min_points=min_points,
        max_points=max_points,
        selection_method=selection_method,
        adj_r_squared_threshold=adj_r_squared_threshold,
    )


def calc_clast_pred(
    intercept: float,
    lambda_z: float,
    tlast: float,
) -> float:
    """
    Calculate predicted Clast from lambda_z regression.

    Parameters
    ----------
    intercept : float
        Y-intercept from log-linear regression
    lambda_z : float
        Terminal rate constant
    tlast : float
        Time of last measurable concentration

    Returns
    -------
    float
        Predicted Clast
    """
    if np.isnan(intercept) or np.isnan(lambda_z) or np.isnan(tlast):
        return np.nan

    return np.exp(intercept - lambda_z * tlast)


def calc_span_ratio(
    half_life: float,
    time_range: Tuple[float, float],
) -> float:
    """
    Calculate span ratio (time range / half-life).

    Parameters
    ----------
    half_life : float
        Terminal half-life
    time_range : tuple
        (first_time, last_time) used in lambda_z calculation

    Returns
    -------
    float
        Span ratio
    """
    if np.isnan(half_life) or half_life <= 0:
        return np.nan

    time_span = time_range[1] - time_range[0]
    return time_span / half_life


def calc_effective_half_life(
    auc_tau: float,
    cmax: float,
    tau: float,
) -> float:
    """
    Calculate effective half-life at steady state.

    The effective half-life is the half-life that would be observed
    based on drug accumulation at steady state. It may differ from
    the terminal half-life due to distribution effects.

    Parameters
    ----------
    auc_tau : float
        AUC over the dosing interval (tau) at steady state
    cmax : float
        Maximum concentration at steady state
    tau : float
        Dosing interval

    Returns
    -------
    float
        Effective half-life

    Notes
    -----
    Calculated as: t1/2,eff = tau * ln(2) / ln(Cmax/(Cmax - AUC_tau/tau))

    References
    ----------
    Boxenbaum H, Battle M. Effective half-life in clinical pharmacology.
    J Clin Pharmacol. 1995;35(8):763-766.
    """
    if (np.isnan(auc_tau) or np.isnan(cmax) or np.isnan(tau) or
            tau <= 0 or cmax <= 0 or auc_tau <= 0):
        return np.nan

    # Average concentration over tau
    cav = auc_tau / tau

    # Check for valid ratio
    if cmax <= cav:
        return np.nan

    # Effective half-life formula
    ratio = cmax / (cmax - cav)
    if ratio <= 1:
        return np.nan

    return tau * np.log(2) / np.log(ratio)


def calc_accumulation_half_life(
    auc_inf_single: float,
    auc_tau_ss: float,
    tau: float,
) -> float:
    """
    Calculate accumulation half-life from single dose and steady state data.

    Parameters
    ----------
    auc_inf_single : float
        AUC from 0 to infinity after single dose
    auc_tau_ss : float
        AUC over dosing interval at steady state
    tau : float
        Dosing interval

    Returns
    -------
    float
        Accumulation half-life

    Notes
    -----
    Based on the relationship: R = AUC_tau_ss / AUC_inf = 1 / (1 - exp(-lambda*tau))
    Rearranging: lambda = -ln(1 - 1/R) / tau
    t1/2 = ln(2) / lambda
    """
    if (np.isnan(auc_inf_single) or np.isnan(auc_tau_ss) or np.isnan(tau) or
            tau <= 0 or auc_inf_single <= 0 or auc_tau_ss <= 0):
        return np.nan

    # Accumulation ratio
    R = auc_tau_ss / auc_inf_single

    if R <= 1:
        # No accumulation or invalid data
        return np.nan

    # Calculate lambda from accumulation ratio
    inner = 1 - 1/R
    if inner <= 0:
        return np.nan

    lambda_eff = -np.log(inner) / tau

    if lambda_eff <= 0:
        return np.nan

    return np.log(2) / lambda_eff


def calc_time_above_threshold(
    conc: np.ndarray,
    time: np.ndarray,
    threshold: float,
    method: str = "linear",
    lambda_z: Optional[float] = None,
) -> float:
    """
    Calculate time that concentration remains above a threshold.

    This is commonly used for:
    - Time above MIC (Minimum Inhibitory Concentration) for antibiotics
    - Time above therapeutic threshold
    - Time above EC50

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    threshold : float
        Threshold concentration (e.g., MIC, EC50)
    method : str
        Interpolation method for finding crossing times
    lambda_z : float, optional
        Terminal elimination rate constant for extrapolation

    Returns
    -------
    float
        Total time above threshold

    Examples
    --------
    >>> conc = np.array([0, 10, 8, 6, 4, 2, 1])
    >>> time = np.array([0, 1, 2, 4, 6, 8, 10])
    >>> calc_time_above_threshold(conc, time, threshold=5)
    3.0  # Time from t=1 to t=4 when conc crosses 5
    """
    conc = np.asarray(conc)
    time = np.asarray(time)

    if len(conc) == 0 or threshold <= 0:
        return np.nan

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    # Find all times when concentration crosses threshold
    above = conc >= threshold
    crossings = []

    # Track entry and exit times
    total_time_above = 0.0
    entry_time = None

    for i in range(len(conc)):
        if i == 0:
            if above[i]:
                # Starts above threshold
                entry_time = time[i]
            continue

        prev_above = above[i - 1]
        curr_above = above[i]

        if not prev_above and curr_above:
            # Crossing up - entering above threshold
            # Interpolate to find exact crossing time
            if conc[i] != conc[i-1]:
                frac = (threshold - conc[i-1]) / (conc[i] - conc[i-1])
                entry_time = time[i-1] + frac * (time[i] - time[i-1])
            else:
                entry_time = time[i]

        elif prev_above and not curr_above:
            # Crossing down - exiting above threshold
            if entry_time is not None:
                # Interpolate to find exact crossing time
                if conc[i] != conc[i-1]:
                    frac = (threshold - conc[i-1]) / (conc[i] - conc[i-1])
                    exit_time = time[i-1] + frac * (time[i] - time[i-1])
                else:
                    exit_time = time[i]
                total_time_above += exit_time - entry_time
                entry_time = None

    # Handle case where concentration is still above threshold at end
    if entry_time is not None:
        # Extrapolate to find when it drops below if lambda_z provided
        if lambda_z is not None and lambda_z > 0:
            # C(t) = Clast * exp(-lambda_z * (t - tlast)) = threshold
            # t = tlast + ln(Clast/threshold) / lambda_z
            clast = conc[-1]
            tlast = time[-1]
            if clast > threshold:
                exit_time = tlast + np.log(clast / threshold) / lambda_z
                total_time_above += exit_time - entry_time
            else:
                # Already below threshold
                total_time_above += tlast - entry_time
        else:
            # Use last observation time
            total_time_above += time[-1] - entry_time

    return total_time_above


def calc_pct_time_above_threshold(
    time_above: float,
    total_time: float,
) -> float:
    """
    Calculate percent of time above threshold.

    Parameters
    ----------
    time_above : float
        Time above threshold
    total_time : float
        Total observation/dosing interval time

    Returns
    -------
    float
        Percent of time above threshold

    Notes
    -----
    Commonly used as %T>MIC for antibiotics.
    """
    if np.isnan(time_above) or np.isnan(total_time) or total_time <= 0:
        return np.nan
    return (time_above / total_time) * 100
