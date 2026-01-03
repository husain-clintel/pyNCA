"""Interval selection and management."""

from typing import Any, List, Optional, Union
from dataclasses import dataclass, field
import numpy as np


@dataclass
class Interval:
    """
    Specification for a calculation interval.

    Parameters
    ----------
    start : float
        Interval start time
    end : float
        Interval end time
    parameters : list of str, optional
        Parameters to calculate in this interval
    subject : any, optional
        Subject this interval applies to (None = all subjects)
    name : str, optional
        Name for this interval
    """

    start: float
    end: float
    parameters: Optional[List[str]] = None
    subject: Optional[Any] = None
    name: Optional[str] = None

    def __post_init__(self):
        if self.start >= self.end:
            raise ValueError(
                f"Interval start ({self.start}) must be less than end ({self.end})"
            )

    def contains(self, time: float) -> bool:
        """Check if a time point is within the interval."""
        return self.start <= time <= self.end

    def overlaps(self, other: "Interval") -> bool:
        """Check if this interval overlaps with another."""
        return self.start < other.end and other.start < self.end

    @property
    def duration(self) -> float:
        """Get interval duration."""
        return self.end - self.start

    def __repr__(self) -> str:
        name_str = f", name='{self.name}'" if self.name else ""
        subj_str = f", subject={self.subject}" if self.subject is not None else ""
        return f"Interval(start={self.start}, end={self.end}{name_str}{subj_str})"


def create_intervals(
    start: Union[float, List[float]],
    end: Union[float, List[float]],
    parameters: Optional[List[str]] = None,
    subjects: Optional[List] = None,
    names: Optional[List[str]] = None,
) -> List[Interval]:
    """
    Create calculation intervals.

    Parameters
    ----------
    start : float or list of float
        Interval start time(s)
    end : float or list of float
        Interval end time(s)
    parameters : list of str, optional
        Parameters to calculate
    subjects : list, optional
        Subjects to apply intervals to
    names : list of str, optional
        Names for intervals

    Returns
    -------
    list of Interval
        Created intervals
    """
    # Normalize to lists
    if not isinstance(start, list):
        start = [start]
    if not isinstance(end, list):
        end = [end]

    if len(start) != len(end):
        raise ValueError("start and end must have same length")

    if names is not None and len(names) != len(start):
        raise ValueError("names must have same length as start/end")

    intervals = []

    for i, (s, e) in enumerate(zip(start, end)):
        name = names[i] if names else None

        if subjects is None:
            intervals.append(
                Interval(start=s, end=e, parameters=parameters, name=name)
            )
        else:
            for subj in subjects:
                intervals.append(
                    Interval(
                        start=s, end=e, parameters=parameters, subject=subj, name=name
                    )
                )

    return intervals


def auto_select_intervals(
    time: np.ndarray,
    dose_times: Optional[np.ndarray] = None,
    tau: Optional[float] = None,
) -> List[Interval]:
    """
    Automatically select calculation intervals based on data.

    Parameters
    ----------
    time : np.ndarray
        Time values
    dose_times : np.ndarray, optional
        Dose administration times
    tau : float, optional
        Dosing interval for multiple doses

    Returns
    -------
    list of Interval
        Auto-selected intervals
    """
    if len(time) == 0:
        return []

    intervals = []

    if dose_times is None or len(dose_times) == 0:
        # Single interval from first to last time
        intervals.append(
            Interval(start=time.min(), end=time.max(), name="0-tlast")
        )
    elif len(dose_times) == 1:
        # Single dose
        intervals.append(
            Interval(start=dose_times[0], end=time.max(), name="single_dose")
        )
    else:
        # Multiple doses
        dose_times = np.sort(dose_times)

        # Calculate tau if not provided
        if tau is None:
            tau = np.median(np.diff(dose_times))

        # Create interval for each dose
        for i, dt in enumerate(dose_times):
            end_time = dt + tau
            if i < len(dose_times) - 1:
                end_time = min(end_time, dose_times[i + 1])

            intervals.append(
                Interval(start=dt, end=end_time, name=f"dose_{i + 1}")
            )

    return intervals


def intervals_for_subject(
    intervals: List[Interval],
    subject: Any,
) -> List[Interval]:
    """
    Get intervals applicable to a specific subject.

    Parameters
    ----------
    intervals : list of Interval
        All intervals
    subject : any
        Subject identifier

    Returns
    -------
    list of Interval
        Applicable intervals
    """
    return [i for i in intervals if i.subject is None or i.subject == subject]
