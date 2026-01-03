"""Concentration extrapolation functions."""

from typing import Optional
import numpy as np


def extrapolate_conc(
    clast: float,
    tlast: float,
    target_time: float,
    lambda_z: float,
) -> float:
    """
    Extrapolate concentration beyond last measurement.

    Uses monoexponential decay: C(t) = Clast * exp(-lambda_z * (t - tlast))

    Parameters
    ----------
    clast : float
        Last measurable concentration
    tlast : float
        Time of last measurable concentration
    target_time : float
        Time point to extrapolate to
    lambda_z : float
        Terminal elimination rate constant

    Returns
    -------
    float
        Extrapolated concentration
    """
    if np.isnan(clast) or np.isnan(tlast) or np.isnan(lambda_z):
        return np.nan

    if lambda_z <= 0:
        return np.nan

    if target_time < tlast:
        return np.nan  # Can't extrapolate backwards

    delta_t = target_time - tlast
    return clast * np.exp(-lambda_z * delta_t)


def extrapolate_conc_from_regression(
    intercept: float,
    lambda_z: float,
    target_time: float,
) -> float:
    """
    Extrapolate concentration using lambda_z regression parameters.

    Uses: C(t) = exp(intercept - lambda_z * t)

    Parameters
    ----------
    intercept : float
        Y-intercept from log-linear regression
    lambda_z : float
        Terminal elimination rate constant (negative of slope)
    target_time : float
        Time point to extrapolate to

    Returns
    -------
    float
        Extrapolated concentration
    """
    if np.isnan(intercept) or np.isnan(lambda_z) or np.isnan(target_time):
        return np.nan

    return np.exp(intercept - lambda_z * target_time)


def extrapolate_auc(
    clast: float,
    lambda_z: float,
) -> float:
    """
    Calculate the extrapolated portion of AUC to infinity.

    AUC_extrap = Clast / lambda_z

    Parameters
    ----------
    clast : float
        Last measurable concentration
    lambda_z : float
        Terminal elimination rate constant

    Returns
    -------
    float
        Extrapolated AUC portion
    """
    if np.isnan(clast) or np.isnan(lambda_z) or lambda_z <= 0:
        return np.nan

    return clast / lambda_z


def extrapolate_aumc(
    clast: float,
    tlast: float,
    lambda_z: float,
) -> float:
    """
    Calculate the extrapolated portion of AUMC to infinity.

    AUMC_extrap = (Clast * tlast) / lambda_z + Clast / lambda_z^2

    Parameters
    ----------
    clast : float
        Last measurable concentration
    tlast : float
        Time of last measurable concentration
    lambda_z : float
        Terminal elimination rate constant

    Returns
    -------
    float
        Extrapolated AUMC portion
    """
    if np.isnan(clast) or np.isnan(tlast) or np.isnan(lambda_z) or lambda_z <= 0:
        return np.nan

    term1 = (clast * tlast) / lambda_z
    term2 = clast / (lambda_z ** 2)

    return term1 + term2


def extrapolate_at_times(
    clast: float,
    tlast: float,
    target_times: np.ndarray,
    lambda_z: float,
) -> np.ndarray:
    """
    Extrapolate concentrations at multiple time points.

    Parameters
    ----------
    clast : float
        Last measurable concentration
    tlast : float
        Time of last measurable concentration
    target_times : array-like
        Time points to extrapolate to
    lambda_z : float
        Terminal elimination rate constant

    Returns
    -------
    np.ndarray
        Extrapolated concentrations
    """
    target_times = np.asarray(target_times, dtype=float)
    result = np.zeros(len(target_times))

    for i, t in enumerate(target_times):
        result[i] = extrapolate_conc(clast, tlast, t, lambda_z)

    return result


def calc_time_to_concentration(
    clast: float,
    tlast: float,
    lambda_z: float,
    target_conc: float,
) -> float:
    """
    Calculate time to reach a target concentration (extrapolation).

    Solves: target_conc = clast * exp(-lambda_z * (t - tlast))

    Parameters
    ----------
    clast : float
        Last measurable concentration
    tlast : float
        Time of last measurable concentration
    lambda_z : float
        Terminal elimination rate constant
    target_conc : float
        Target concentration to reach

    Returns
    -------
    float
        Time to reach target concentration
    """
    if np.isnan(clast) or np.isnan(tlast) or np.isnan(lambda_z) or np.isnan(target_conc):
        return np.nan

    if lambda_z <= 0 or clast <= 0 or target_conc <= 0:
        return np.nan

    if target_conc >= clast:
        return np.nan  # Can't extrapolate to higher concentration

    # t = tlast - ln(target_conc / clast) / lambda_z
    delta_t = -np.log(target_conc / clast) / lambda_z
    return tlast + delta_t
