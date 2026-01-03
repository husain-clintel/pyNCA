"""Bioavailability calculations."""

from typing import Optional
import numpy as np


def calc_f(
    auc_test: float,
    auc_ref: float,
    dose_test: float,
    dose_ref: float,
) -> float:
    """
    Calculate relative or absolute bioavailability (F).

    F = (AUC_test / AUC_ref) * (Dose_ref / Dose_test)

    Parameters
    ----------
    auc_test : float
        AUC for test formulation/route
    auc_ref : float
        AUC for reference formulation/route (IV for absolute F)
    dose_test : float
        Dose of test formulation
    dose_ref : float
        Dose of reference formulation

    Returns
    -------
    float
        Bioavailability (fractional, multiply by 100 for %)
    """
    if np.isnan(auc_test) or np.isnan(auc_ref):
        return np.nan

    if auc_ref <= 0 or dose_test <= 0:
        return np.nan

    return (auc_test / auc_ref) * (dose_ref / dose_test)


def calc_f_from_cl(
    cl_test: float,
    cl_ref: float,
) -> float:
    """
    Calculate bioavailability from clearance values.

    F = CL_ref / CL_test

    This assumes the same true clearance for both routes.

    Parameters
    ----------
    cl_test : float
        Apparent clearance for test route (CL/F)
    cl_ref : float
        Clearance for reference route (typically IV)

    Returns
    -------
    float
        Bioavailability
    """
    if np.isnan(cl_test) or np.isnan(cl_ref):
        return np.nan

    if cl_test <= 0:
        return np.nan

    return cl_ref / cl_test


def calc_f_percent(f: float) -> float:
    """
    Convert fractional bioavailability to percentage.

    Parameters
    ----------
    f : float
        Fractional bioavailability

    Returns
    -------
    float
        Bioavailability as percentage
    """
    if np.isnan(f):
        return np.nan

    return f * 100


def calc_relative_f(
    auc_test: float,
    auc_ref: float,
    dose_test: Optional[float] = None,
    dose_ref: Optional[float] = None,
) -> float:
    """
    Calculate relative bioavailability (Frel).

    If doses are equal, Frel = AUC_test / AUC_ref
    Otherwise, Frel = (AUC_test / AUC_ref) * (Dose_ref / Dose_test)

    Parameters
    ----------
    auc_test : float
        AUC for test formulation
    auc_ref : float
        AUC for reference formulation
    dose_test : float, optional
        Dose of test formulation
    dose_ref : float, optional
        Dose of reference formulation

    Returns
    -------
    float
        Relative bioavailability
    """
    if np.isnan(auc_test) or np.isnan(auc_ref) or auc_ref <= 0:
        return np.nan

    if dose_test is None or dose_ref is None:
        # Assume equal doses
        return auc_test / auc_ref

    if dose_test <= 0:
        return np.nan

    return (auc_test / auc_ref) * (dose_ref / dose_test)


def calc_accumulation_index(
    half_life: float,
    tau: float,
) -> float:
    """
    Calculate accumulation index (Rac).

    Rac = 1 / (1 - exp(-lambda_z * tau))
    or equivalently:
    Rac = 1 / (1 - exp(-ln(2) * tau / t1/2))

    Parameters
    ----------
    half_life : float
        Terminal half-life
    tau : float
        Dosing interval

    Returns
    -------
    float
        Accumulation index
    """
    if np.isnan(half_life) or np.isnan(tau):
        return np.nan

    if half_life <= 0 or tau <= 0:
        return np.nan

    lambda_z = np.log(2) / half_life
    exp_term = np.exp(-lambda_z * tau)

    if exp_term >= 1:
        return np.nan

    return 1 / (1 - exp_term)


def calc_accumulation_index_from_lambda_z(
    lambda_z: float,
    tau: float,
) -> float:
    """
    Calculate accumulation index from lambda_z.

    Rac = 1 / (1 - exp(-lambda_z * tau))

    Parameters
    ----------
    lambda_z : float
        Terminal elimination rate constant
    tau : float
        Dosing interval

    Returns
    -------
    float
        Accumulation index
    """
    if np.isnan(lambda_z) or np.isnan(tau):
        return np.nan

    if lambda_z <= 0 or tau <= 0:
        return np.nan

    exp_term = np.exp(-lambda_z * tau)

    if exp_term >= 1:
        return np.nan

    return 1 / (1 - exp_term)


def calc_linearity_factor(
    auc_dose_n: float,
    auc_dose_1: float,
    dose_n: float,
    dose_1: float,
) -> float:
    """
    Calculate dose linearity factor.

    Linearity = (AUC_n / AUC_1) / (Dose_n / Dose_1)

    A value of 1.0 indicates dose-linear pharmacokinetics.

    Parameters
    ----------
    auc_dose_n : float
        AUC at higher dose
    auc_dose_1 : float
        AUC at lower/reference dose
    dose_n : float
        Higher dose
    dose_1 : float
        Lower/reference dose

    Returns
    -------
    float
        Linearity factor (1.0 = linear)
    """
    if np.isnan(auc_dose_n) or np.isnan(auc_dose_1):
        return np.nan

    if auc_dose_1 <= 0 or dose_1 <= 0 or dose_n <= 0:
        return np.nan

    auc_ratio = auc_dose_n / auc_dose_1
    dose_ratio = dose_n / dose_1

    return auc_ratio / dose_ratio
