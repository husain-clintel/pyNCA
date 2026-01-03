"""Sparse NCA calculations for destructive sampling."""

from typing import Optional, Tuple, Dict
import numpy as np
from scipy import stats


def sparse_auc(
    conc: np.ndarray,
    time: np.ndarray,
    subject: np.ndarray,
    method: str = "batch",
) -> Dict[str, float]:
    """
    Calculate AUC for sparse/destructive sampling data.

    In sparse sampling, each subject contributes only one concentration
    measurement (destructive sampling). AUC is calculated from mean
    concentrations at each time point.

    Parameters
    ----------
    conc : array-like
        Concentration values (one per subject)
    time : array-like
        Time values (one per subject)
    subject : array-like
        Subject identifiers
    method : str
        Calculation method:
        - "batch": Batch means method
        - "bailer": Bailer's method for variance estimation

    Returns
    -------
    dict
        Dictionary with:
        - auc: AUC estimate
        - se: Standard error
        - cv: Coefficient of variation
        - n_subjects: Number of subjects
        - n_timepoints: Number of unique time points
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)
    subject = np.asarray(subject)

    # Get unique time points
    unique_times = np.sort(np.unique(time))
    n_times = len(unique_times)

    if n_times < 2:
        return {
            "auc": np.nan,
            "se": np.nan,
            "cv": np.nan,
            "n_subjects": len(np.unique(subject)),
            "n_timepoints": n_times,
        }

    # Calculate mean and variance at each time point
    means = []
    variances = []
    ns = []

    for t in unique_times:
        mask = time == t
        c_at_t = conc[mask]
        c_at_t = c_at_t[~np.isnan(c_at_t)]

        n = len(c_at_t)
        if n > 0:
            means.append(np.mean(c_at_t))
            variances.append(np.var(c_at_t, ddof=1) if n > 1 else 0)
            ns.append(n)
        else:
            means.append(0)
            variances.append(0)
            ns.append(0)

    means = np.array(means)
    variances = np.array(variances)
    ns = np.array(ns)

    # Calculate AUC using trapezoidal rule on mean concentrations
    auc = 0
    for i in range(n_times - 1):
        dt = unique_times[i + 1] - unique_times[i]
        auc += (means[i] + means[i + 1]) * dt / 2

    # Calculate variance of AUC
    if method == "batch":
        se = _batch_means_se(means, variances, ns, unique_times)
    elif method == "bailer":
        se = _bailer_se(means, variances, ns, unique_times)
    else:
        raise ValueError(f"Unknown method: {method}")

    cv = (se / auc * 100) if auc > 0 else np.nan

    return {
        "auc": auc,
        "se": se,
        "cv": cv,
        "n_subjects": len(np.unique(subject)),
        "n_timepoints": n_times,
    }


def _batch_means_se(
    means: np.ndarray,
    variances: np.ndarray,
    ns: np.ndarray,
    times: np.ndarray,
) -> float:
    """Calculate SE using batch means method."""
    n_times = len(times)
    var_auc = 0

    for i in range(n_times - 1):
        dt = times[i + 1] - times[i]

        # Variance contribution from each endpoint
        if ns[i] > 0:
            var_auc += (dt / 2) ** 2 * variances[i] / ns[i]
        if ns[i + 1] > 0:
            var_auc += (dt / 2) ** 2 * variances[i + 1] / ns[i + 1]

    return np.sqrt(var_auc)


def _bailer_se(
    means: np.ndarray,
    variances: np.ndarray,
    ns: np.ndarray,
    times: np.ndarray,
) -> float:
    """
    Calculate SE using Bailer's method.

    Bailer WA. Testing for the equality of area under the curves when
    using destructive measurement techniques. J Pharmacokinet Biopharm.
    1988;16(3):303-309.
    """
    n_times = len(times)
    var_auc = 0

    for i in range(n_times):
        # Weight for this time point
        if i == 0:
            w = (times[1] - times[0]) / 2
        elif i == n_times - 1:
            w = (times[-1] - times[-2]) / 2
        else:
            w = (times[i + 1] - times[i - 1]) / 2

        if ns[i] > 0:
            var_auc += w ** 2 * variances[i] / ns[i]

    return np.sqrt(var_auc)


def sparse_auc_se(
    conc: np.ndarray,
    time: np.ndarray,
    subject: np.ndarray,
) -> float:
    """
    Calculate standard error of sparse AUC.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    subject : array-like
        Subject identifiers

    Returns
    -------
    float
        Standard error of AUC
    """
    result = sparse_auc(conc, time, subject)
    return result["se"]


def sparse_cmax(
    conc: np.ndarray,
    time: np.ndarray,
    subject: np.ndarray,
) -> Dict[str, float]:
    """
    Calculate Cmax statistics for sparse data.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    subject : array-like
        Subject identifiers

    Returns
    -------
    dict
        Dictionary with cmax, tmax, se, cv
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)

    unique_times = np.sort(np.unique(time))

    means = []
    variances = []
    ns = []

    for t in unique_times:
        mask = time == t
        c_at_t = conc[mask]
        c_at_t = c_at_t[~np.isnan(c_at_t)]

        n = len(c_at_t)
        if n > 0:
            means.append(np.mean(c_at_t))
            variances.append(np.var(c_at_t, ddof=1) if n > 1 else 0)
            ns.append(n)

    means = np.array(means)
    variances = np.array(variances)
    ns = np.array(ns)

    # Find Cmax
    if len(means) == 0:
        return {"cmax": np.nan, "tmax": np.nan, "se": np.nan, "cv": np.nan}

    cmax_idx = np.argmax(means)
    cmax = means[cmax_idx]
    tmax = unique_times[cmax_idx]

    # SE of Cmax
    se = np.sqrt(variances[cmax_idx] / ns[cmax_idx]) if ns[cmax_idx] > 0 else np.nan
    cv = (se / cmax * 100) if cmax > 0 else np.nan

    return {
        "cmax": cmax,
        "tmax": tmax,
        "se": se,
        "cv": cv,
    }
