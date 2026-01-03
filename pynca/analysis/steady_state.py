"""Steady-state assessment functions."""

from typing import Dict, List, Optional, Tuple
import numpy as np
from scipy import stats


def time_to_steady_state(
    conc: np.ndarray,
    time: np.ndarray,
    dose_times: np.ndarray,
    method: str = "monoexponential",
    threshold: float = 0.90,
) -> Dict[str, float]:
    """
    Estimate time to reach steady state.

    Parameters
    ----------
    conc : array-like
        Concentration values (typically trough concentrations)
    time : array-like
        Time values
    dose_times : array-like
        Times of dose administration
    method : str
        Estimation method:
        - "monoexponential": Fit approach to steady state
        - "linear_regression": Use linear regression to detect plateau
    threshold : float
        Fraction of steady state to consider as "reached" (default 0.90)

    Returns
    -------
    dict
        Dictionary with:
        - tss: Time to steady state
        - css: Estimated steady-state concentration
        - r_squared: Goodness of fit
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)

    # Remove NaN
    valid = ~np.isnan(conc)
    conc = conc[valid]
    time = time[valid]

    if len(conc) < 3:
        return {"tss": np.nan, "css": np.nan, "r_squared": np.nan}

    if method == "monoexponential":
        return _tss_monoexponential(conc, time, threshold)
    elif method == "linear_regression":
        return _tss_linear_regression(conc, time)
    else:
        raise ValueError(f"Unknown method: {method}")


def _tss_monoexponential(
    conc: np.ndarray,
    time: np.ndarray,
    threshold: float,
) -> Dict[str, float]:
    """
    Estimate TSS using monoexponential approach.

    Model: C(t) = Css * (1 - exp(-k * t))
    """
    # Estimate Css from last few points
    css_est = np.mean(conc[-3:])

    # Transform for linear regression
    # ln(Css - C) = ln(Css) - k*t
    # Only use points where C < Css
    mask = conc < css_est * 0.99
    if np.sum(mask) < 2:
        return {"tss": np.nan, "css": css_est, "r_squared": np.nan}

    y = np.log(css_est - conc[mask])
    x = time[mask]

    try:
        slope, intercept, r_value, _, _ = stats.linregress(x, y)
        k = -slope

        if k <= 0:
            return {"tss": np.nan, "css": css_est, "r_squared": r_value ** 2}

        # Time to reach threshold of steady state
        # threshold = 1 - exp(-k * tss)
        # tss = -ln(1 - threshold) / k
        tss = -np.log(1 - threshold) / k

        return {
            "tss": tss,
            "css": css_est,
            "r_squared": r_value ** 2,
            "k_approach": k,
        }
    except Exception:
        return {"tss": np.nan, "css": css_est, "r_squared": np.nan}


def _tss_linear_regression(
    conc: np.ndarray,
    time: np.ndarray,
) -> Dict[str, float]:
    """
    Estimate TSS using linear regression to detect plateau.

    Fits linear regression to sequential subsets and finds where
    slope becomes non-significant.
    """
    n = len(conc)
    window = max(3, n // 3)

    for start in range(n - window + 1):
        subset_conc = conc[start : start + window]
        subset_time = time[start : start + window]

        slope, intercept, r_value, p_value, _ = stats.linregress(
            subset_time, subset_conc
        )

        # If slope is not significantly different from 0
        if p_value > 0.05:
            return {
                "tss": time[start],
                "css": np.mean(subset_conc),
                "r_squared": r_value ** 2,
            }

    # Steady state not reached
    return {"tss": np.nan, "css": np.mean(conc[-3:]), "r_squared": np.nan}


def check_steady_state(
    trough_concs: np.ndarray,
    trough_times: np.ndarray,
    method: str = "tost",
    bioequivalence_limits: Tuple[float, float] = (0.80, 1.25),
    alpha: float = 0.05,
) -> Dict[str, any]:
    """
    Test if steady state has been achieved.

    Parameters
    ----------
    trough_concs : array-like
        Trough concentration values from multiple doses
    trough_times : array-like
        Times of trough measurements
    method : str
        Testing method:
        - "tost": Two One-Sided Tests
        - "regression": Test slope = 0
        - "cv": Coefficient of variation threshold
    bioequivalence_limits : tuple
        Lower and upper limits for TOST (default 0.80, 1.25)
    alpha : float
        Significance level

    Returns
    -------
    dict
        Dictionary with:
        - steady_state: bool indicating if SS achieved
        - p_value: p-value (method dependent)
        - details: Additional statistics
    """
    trough_concs = np.asarray(trough_concs, dtype=float)
    trough_times = np.asarray(trough_times, dtype=float)

    # Remove NaN
    valid = ~np.isnan(trough_concs)
    trough_concs = trough_concs[valid]
    trough_times = trough_times[valid]

    if len(trough_concs) < 3:
        return {
            "steady_state": False,
            "p_value": np.nan,
            "details": "Insufficient data",
        }

    if method == "tost":
        return _check_ss_tost(trough_concs, bioequivalence_limits, alpha)
    elif method == "regression":
        return _check_ss_regression(trough_concs, trough_times, alpha)
    elif method == "cv":
        return _check_ss_cv(trough_concs)
    else:
        raise ValueError(f"Unknown method: {method}")


def _check_ss_tost(
    trough_concs: np.ndarray,
    limits: Tuple[float, float],
    alpha: float,
) -> Dict[str, any]:
    """Check steady state using TOST on consecutive ratios."""
    n = len(trough_concs)

    # Calculate ratios of consecutive troughs
    ratios = trough_concs[1:] / trough_concs[:-1]

    # Log-transform for testing
    log_ratios = np.log(ratios)
    mean_log = np.mean(log_ratios)
    se_log = np.std(log_ratios, ddof=1) / np.sqrt(len(log_ratios))

    # TOST
    lower_limit = np.log(limits[0])
    upper_limit = np.log(limits[1])

    df = len(log_ratios) - 1
    t_crit = stats.t.ppf(1 - alpha, df)

    # 90% CI
    ci_lower = mean_log - t_crit * se_log
    ci_upper = mean_log + t_crit * se_log

    # Steady state if CI within limits
    steady_state = (ci_lower > lower_limit) and (ci_upper < upper_limit)

    return {
        "steady_state": steady_state,
        "p_value": np.nan,  # TOST doesn't have single p-value
        "details": {
            "geometric_mean_ratio": np.exp(mean_log),
            "ci_90_lower": np.exp(ci_lower),
            "ci_90_upper": np.exp(ci_upper),
            "n_ratios": len(ratios),
        },
    }


def _check_ss_regression(
    trough_concs: np.ndarray,
    trough_times: np.ndarray,
    alpha: float,
) -> Dict[str, any]:
    """Check steady state by testing if slope = 0."""
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        trough_times, trough_concs
    )

    steady_state = p_value > alpha

    return {
        "steady_state": steady_state,
        "p_value": p_value,
        "details": {
            "slope": slope,
            "intercept": intercept,
            "r_squared": r_value ** 2,
        },
    }


def _check_ss_cv(
    trough_concs: np.ndarray,
    cv_threshold: float = 20.0,
) -> Dict[str, any]:
    """Check steady state based on CV threshold."""
    cv = np.std(trough_concs, ddof=1) / np.mean(trough_concs) * 100
    steady_state = cv < cv_threshold

    return {
        "steady_state": steady_state,
        "p_value": np.nan,
        "details": {
            "cv_percent": cv,
            "threshold": cv_threshold,
        },
    }
