"""Superposition calculations for predicting steady-state concentrations."""

from typing import Optional, Tuple
import numpy as np


def superposition(
    conc: np.ndarray,
    time: np.ndarray,
    tau: float,
    n_doses: Optional[int] = None,
    method: str = "concentration",
    lambda_z: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Predict concentrations at steady state using superposition principle.

    The superposition principle assumes linear pharmacokinetics and
    predicts steady-state concentrations by summing contributions from
    multiple doses.

    Parameters
    ----------
    conc : array-like
        Single-dose concentration values
    time : array-like
        Time values (relative to dose time)
    tau : float
        Dosing interval
    n_doses : int, optional
        Number of doses to superimpose. If None, calculated to reach
        approximate steady state (5 half-lives)
    method : str
        Superposition method:
        - "concentration": Sum concentrations at each time point
        - "auc": Calculate steady-state AUC
    lambda_z : float, optional
        Terminal elimination rate constant (required if n_doses is None)

    Returns
    -------
    tuple
        (steady_state_conc, time_within_tau)

    Notes
    -----
    Assumes:
    - Linear pharmacokinetics (dose-proportional)
    - Time-invariant kinetics
    - Complete elimination between doses not required
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    # Determine number of doses needed
    if n_doses is None:
        if lambda_z is not None and lambda_z > 0:
            half_life = np.log(2) / lambda_z
            # 5-7 half-lives to reach ~97-99% of steady state
            n_doses = max(1, int(np.ceil(7 * half_life / tau)))
        else:
            n_doses = 10  # Default

    # Create time points within tau
    mask = time <= tau
    time_tau = time[mask]
    conc_tau = conc[mask]

    if len(time_tau) == 0:
        return np.array([]), np.array([])

    # Initialize steady-state concentrations
    ss_conc = np.zeros_like(conc_tau)

    # Superimpose contributions from each dose
    for dose_num in range(n_doses):
        dose_time = dose_num * tau

        for i, t in enumerate(time_tau):
            # Time since this dose
            t_since_dose = t + dose_time

            # Interpolate/extrapolate concentration at this time
            if t_since_dose <= time[-1]:
                # Interpolate from single-dose data
                c = np.interp(t_since_dose, time, conc)
            elif lambda_z is not None and lambda_z > 0:
                # Extrapolate using terminal phase
                clast = conc[-1]
                tlast = time[-1]
                c = clast * np.exp(-lambda_z * (t_since_dose - tlast))
            else:
                c = 0

            ss_conc[i] += c

    return ss_conc, time_tau


def superposition_auc(
    auc_single: float,
    tau: float,
    lambda_z: float,
) -> float:
    """
    Calculate steady-state AUC using superposition.

    At steady state: AUC_tau_ss = AUC_inf_single * (1 / (1 - exp(-lambda_z * tau)))
    For tau -> infinity: AUC_tau_ss = AUC_inf_single

    Parameters
    ----------
    auc_single : float
        AUC from single dose (AUCinf)
    tau : float
        Dosing interval
    lambda_z : float
        Terminal elimination rate constant

    Returns
    -------
    float
        Steady-state AUC over one dosing interval
    """
    if np.isnan(auc_single) or np.isnan(tau) or np.isnan(lambda_z):
        return np.nan

    if lambda_z <= 0 or tau <= 0:
        return np.nan

    # Accumulation factor
    accumulation = 1 / (1 - np.exp(-lambda_z * tau))

    return auc_single * accumulation


def superposition_cmax(
    conc: np.ndarray,
    time: np.ndarray,
    tau: float,
    lambda_z: float,
    n_doses: Optional[int] = None,
) -> float:
    """
    Predict Cmax at steady state.

    Parameters
    ----------
    conc : array-like
        Single-dose concentration values
    time : array-like
        Time values
    tau : float
        Dosing interval
    lambda_z : float
        Terminal elimination rate constant
    n_doses : int, optional
        Number of doses to simulate

    Returns
    -------
    float
        Predicted steady-state Cmax
    """
    ss_conc, _ = superposition(
        conc, time, tau,
        n_doses=n_doses,
        lambda_z=lambda_z,
    )

    if len(ss_conc) == 0:
        return np.nan

    return np.max(ss_conc)


def superposition_cmin(
    conc: np.ndarray,
    time: np.ndarray,
    tau: float,
    lambda_z: float,
    n_doses: Optional[int] = None,
) -> float:
    """
    Predict Cmin (trough) at steady state.

    Parameters
    ----------
    conc : array-like
        Single-dose concentration values
    time : array-like
        Time values
    tau : float
        Dosing interval
    lambda_z : float
        Terminal elimination rate constant
    n_doses : int, optional
        Number of doses to simulate

    Returns
    -------
    float
        Predicted steady-state Cmin
    """
    ss_conc, ss_time = superposition(
        conc, time, tau,
        n_doses=n_doses,
        lambda_z=lambda_z,
    )

    if len(ss_conc) == 0:
        return np.nan

    # Cmin is typically at end of interval
    return ss_conc[-1]


def accumulation_ratio(
    tau: float,
    lambda_z: float,
) -> float:
    """
    Calculate theoretical accumulation ratio.

    R = 1 / (1 - exp(-lambda_z * tau))

    Parameters
    ----------
    tau : float
        Dosing interval
    lambda_z : float
        Terminal elimination rate constant

    Returns
    -------
    float
        Accumulation ratio
    """
    if np.isnan(tau) or np.isnan(lambda_z):
        return np.nan

    if lambda_z <= 0 or tau <= 0:
        return np.nan

    return 1 / (1 - np.exp(-lambda_z * tau))


def time_to_steady_state_doses(
    lambda_z: float,
    tau: float,
    fraction: float = 0.90,
) -> int:
    """
    Calculate number of doses to reach a fraction of steady state.

    Parameters
    ----------
    lambda_z : float
        Terminal elimination rate constant
    tau : float
        Dosing interval
    fraction : float
        Fraction of steady state (default 0.90 = 90%)

    Returns
    -------
    int
        Number of doses required
    """
    if np.isnan(lambda_z) or np.isnan(tau):
        return np.nan

    if lambda_z <= 0 or tau <= 0 or fraction <= 0 or fraction >= 1:
        return np.nan

    # Fraction achieved after n doses: f = 1 - exp(-n * lambda_z * tau)
    # Solving for n: n = -ln(1 - fraction) / (lambda_z * tau)
    n = -np.log(1 - fraction) / (lambda_z * tau)

    return int(np.ceil(n))
