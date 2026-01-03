"""Clearance and volume of distribution calculations."""

from typing import Optional
import numpy as np

from pynca.utils.validation import validate_positive


def calc_cl(
    dose: float,
    auc_inf: float,
    f: float = 1.0,
) -> float:
    """
    Calculate total clearance (CL or CL/F).

    CL = Dose / AUCinf (for IV)
    CL/F = Dose / AUCinf (for extravascular, where F is bioavailability)

    Parameters
    ----------
    dose : float
        Dose amount
    auc_inf : float
        AUC extrapolated to infinity
    f : float
        Bioavailability factor (default 1.0 for IV)

    Returns
    -------
    float
        Clearance value (CL or CL/F depending on route)
    """
    if np.isnan(dose) or np.isnan(auc_inf) or auc_inf <= 0:
        return np.nan

    if dose <= 0:
        return np.nan

    return (dose * f) / auc_inf


def calc_cl_last(
    dose: float,
    auc_last: float,
    f: float = 1.0,
) -> float:
    """
    Calculate clearance based on AUClast.

    Parameters
    ----------
    dose : float
        Dose amount
    auc_last : float
        AUC to last measurable concentration
    f : float
        Bioavailability factor

    Returns
    -------
    float
        Clearance value
    """
    if np.isnan(dose) or np.isnan(auc_last) or auc_last <= 0:
        return np.nan

    if dose <= 0:
        return np.nan

    return (dose * f) / auc_last


def calc_vz(
    cl: float,
    lambda_z: float,
) -> float:
    """
    Calculate volume of distribution based on terminal phase (Vz or Vz/F).

    Vz = CL / lambda_z

    Parameters
    ----------
    cl : float
        Clearance (CL or CL/F)
    lambda_z : float
        Terminal elimination rate constant

    Returns
    -------
    float
        Volume of distribution (Vz or Vz/F)
    """
    if np.isnan(cl) or np.isnan(lambda_z) or lambda_z <= 0:
        return np.nan

    return cl / lambda_z


def calc_vz_alt(
    dose: float,
    auc_inf: float,
    lambda_z: float,
    f: float = 1.0,
) -> float:
    """
    Calculate Vz directly from dose, AUC, and lambda_z.

    Vz = Dose / (AUCinf * lambda_z)

    Parameters
    ----------
    dose : float
        Dose amount
    auc_inf : float
        AUC extrapolated to infinity
    lambda_z : float
        Terminal elimination rate constant
    f : float
        Bioavailability factor

    Returns
    -------
    float
        Volume of distribution
    """
    if np.isnan(dose) or np.isnan(auc_inf) or np.isnan(lambda_z):
        return np.nan

    if auc_inf <= 0 or lambda_z <= 0 or dose <= 0:
        return np.nan

    return (dose * f) / (auc_inf * lambda_z)


def calc_vss(
    dose: float,
    aumc_inf: float,
    auc_inf: float,
    f: float = 1.0,
) -> float:
    """
    Calculate volume of distribution at steady state (Vss or Vss/F).

    Vss = Dose * AUMC_inf / (AUC_inf)^2
    or equivalently:
    Vss = MRT * CL

    Parameters
    ----------
    dose : float
        Dose amount
    aumc_inf : float
        AUMC extrapolated to infinity
    auc_inf : float
        AUC extrapolated to infinity
    f : float
        Bioavailability factor

    Returns
    -------
    float
        Vss (or Vss/F for extravascular)
    """
    if np.isnan(dose) or np.isnan(aumc_inf) or np.isnan(auc_inf):
        return np.nan

    if auc_inf <= 0 or dose <= 0:
        return np.nan

    return (dose * f * aumc_inf) / (auc_inf ** 2)


def calc_vss_from_mrt(
    mrt: float,
    cl: float,
) -> float:
    """
    Calculate Vss from MRT and clearance.

    Vss = MRT * CL

    Parameters
    ----------
    mrt : float
        Mean residence time
    cl : float
        Clearance

    Returns
    -------
    float
        Volume of distribution at steady state
    """
    if np.isnan(mrt) or np.isnan(cl):
        return np.nan

    if mrt <= 0 or cl <= 0:
        return np.nan

    return mrt * cl


def calc_mrt(
    aumc_inf: float,
    auc_inf: float,
) -> float:
    """
    Calculate mean residence time (MRT).

    MRT = AUMC_inf / AUC_inf

    Parameters
    ----------
    aumc_inf : float
        AUMC extrapolated to infinity
    auc_inf : float
        AUC extrapolated to infinity

    Returns
    -------
    float
        Mean residence time
    """
    if np.isnan(aumc_inf) or np.isnan(auc_inf) or auc_inf <= 0:
        return np.nan

    return aumc_inf / auc_inf


def calc_mrt_last(
    aumc_last: float,
    auc_last: float,
) -> float:
    """
    Calculate MRT based on AUMClast and AUClast.

    Parameters
    ----------
    aumc_last : float
        AUMC to last measurable concentration
    auc_last : float
        AUC to last measurable concentration

    Returns
    -------
    float
        Mean residence time to Tlast
    """
    if np.isnan(aumc_last) or np.isnan(auc_last) or auc_last <= 0:
        return np.nan

    return aumc_last / auc_last


def calc_mrt_iv(
    aumc_inf: float,
    auc_inf: float,
    duration: float = 0.0,
) -> float:
    """
    Calculate MRT corrected for IV infusion duration.

    MRT_iv = (AUMC_inf / AUC_inf) - (duration / 2)

    Parameters
    ----------
    aumc_inf : float
        AUMC extrapolated to infinity
    auc_inf : float
        AUC extrapolated to infinity
    duration : float
        Infusion duration (0 for bolus)

    Returns
    -------
    float
        Corrected mean residence time
    """
    mrt = calc_mrt(aumc_inf, auc_inf)

    if np.isnan(mrt):
        return np.nan

    return mrt - (duration / 2)


def calc_mat(
    mrt_po: float,
    mrt_iv: float,
) -> float:
    """
    Calculate mean absorption time (MAT).

    MAT = MRT_po - MRT_iv

    Parameters
    ----------
    mrt_po : float
        Mean residence time after oral administration
    mrt_iv : float
        Mean residence time after IV administration

    Returns
    -------
    float
        Mean absorption time
    """
    if np.isnan(mrt_po) or np.isnan(mrt_iv):
        return np.nan

    return mrt_po - mrt_iv


def calc_kel(
    cl: float,
    vz: float,
) -> float:
    """
    Calculate elimination rate constant from CL and Vz.

    kel = CL / Vz

    Note: This should equal lambda_z if the model is monoexponential.

    Parameters
    ----------
    cl : float
        Clearance
    vz : float
        Volume of distribution

    Returns
    -------
    float
        Elimination rate constant
    """
    if np.isnan(cl) or np.isnan(vz) or vz <= 0:
        return np.nan

    return cl / vz
