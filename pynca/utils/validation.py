"""Input validation utilities."""

from typing import Union, Optional
import numpy as np
import pandas as pd


def validate_concentration_data(
    conc: Union[np.ndarray, pd.Series, list],
    allow_negative: bool = False,
) -> np.ndarray:
    """
    Validate concentration data.

    Parameters
    ----------
    conc : array-like
        Concentration values
    allow_negative : bool
        Whether to allow negative concentrations

    Returns
    -------
    np.ndarray
        Validated concentration array

    Raises
    ------
    ValueError
        If validation fails
    """
    conc = np.asarray(conc, dtype=float)

    if conc.ndim != 1:
        raise ValueError("Concentration data must be 1-dimensional")

    if len(conc) == 0:
        raise ValueError("Concentration data cannot be empty")

    if not allow_negative and np.any(conc[~np.isnan(conc)] < 0):
        raise ValueError("Concentration values cannot be negative")

    return conc


def validate_time_data(
    time: Union[np.ndarray, pd.Series, list],
    allow_negative: bool = False,
) -> np.ndarray:
    """
    Validate time data.

    Parameters
    ----------
    time : array-like
        Time values
    allow_negative : bool
        Whether to allow negative times (pre-dose)

    Returns
    -------
    np.ndarray
        Validated time array

    Raises
    ------
    ValueError
        If validation fails
    """
    time = np.asarray(time, dtype=float)

    if time.ndim != 1:
        raise ValueError("Time data must be 1-dimensional")

    if len(time) == 0:
        raise ValueError("Time data cannot be empty")

    if np.any(np.isnan(time)):
        raise ValueError("Time values cannot contain NaN")

    if not allow_negative and np.any(time < 0):
        raise ValueError("Time values cannot be negative")

    return time


def validate_conc_time_match(
    conc: np.ndarray,
    time: np.ndarray,
) -> None:
    """
    Validate that concentration and time arrays have matching lengths.

    Parameters
    ----------
    conc : np.ndarray
        Concentration values
    time : np.ndarray
        Time values

    Raises
    ------
    ValueError
        If arrays have different lengths
    """
    if len(conc) != len(time):
        raise ValueError(
            f"Concentration and time arrays must have same length. "
            f"Got {len(conc)} and {len(time)}"
        )


def validate_positive(
    value: Union[float, np.ndarray],
    name: str = "value",
) -> None:
    """
    Validate that a value is positive.

    Parameters
    ----------
    value : float or array-like
        Value to validate
    name : str
        Name of the value for error messages

    Raises
    ------
    ValueError
        If value is not positive
    """
    value = np.asarray(value)
    if np.any(value <= 0):
        raise ValueError(f"{name} must be positive")


def validate_non_negative(
    value: Union[float, np.ndarray],
    name: str = "value",
) -> None:
    """
    Validate that a value is non-negative.

    Parameters
    ----------
    value : float or array-like
        Value to validate
    name : str
        Name of the value for error messages

    Raises
    ------
    ValueError
        If value is negative
    """
    value = np.asarray(value)
    if np.any(value < 0):
        raise ValueError(f"{name} cannot be negative")


def validate_interval(
    start: float,
    end: float,
) -> None:
    """
    Validate a time interval.

    Parameters
    ----------
    start : float
        Interval start time
    end : float
        Interval end time

    Raises
    ------
    ValueError
        If interval is invalid
    """
    if start >= end:
        raise ValueError(f"Interval start ({start}) must be less than end ({end})")


def validate_lambda_z(lambda_z: float) -> None:
    """
    Validate terminal rate constant.

    Parameters
    ----------
    lambda_z : float
        Terminal elimination rate constant

    Raises
    ------
    ValueError
        If lambda_z is invalid
    """
    if lambda_z is None or np.isnan(lambda_z):
        raise ValueError("lambda_z is not available")
    if lambda_z <= 0:
        raise ValueError("lambda_z must be positive")
