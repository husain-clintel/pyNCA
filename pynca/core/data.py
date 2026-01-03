"""Combined NCA data class."""

from typing import Dict, List, Optional, Union, Iterator, Tuple, TYPE_CHECKING
import numpy as np
import pandas as pd

from pynca.core.concentration import NCAConcentration
from pynca.core.dose import NCADose

if TYPE_CHECKING:
    from pynca.analysis.intervals import Interval


class NCAData:
    """
    Combined concentration and dose data for NCA analysis.

    This class combines concentration-time data with dosing information
    and determines calculation intervals, analogous to PKNCAdata in R.

    Parameters
    ----------
    conc : NCAConcentration
        Concentration-time data
    dose : NCADose, optional
        Dosing information
    intervals : list of Interval, optional
        Calculation intervals. If None, automatically determined.

    Attributes
    ----------
    conc : NCAConcentration
        Concentration data
    dose : NCADose
        Dose data
    intervals : list
        Calculation intervals

    Examples
    --------
    >>> import pynca as nca
    >>> data = nca.NCAData(conc=conc_obj, dose=dose_obj)
    """

    def __init__(
        self,
        conc: NCAConcentration,
        dose: Optional[NCADose] = None,
        intervals: Optional[List["Interval"]] = None,
    ):
        self._validate_inputs(conc, dose)

        self.conc = conc
        self.dose = dose
        self._intervals = intervals

        # Auto-select intervals if not provided
        if self._intervals is None:
            self._intervals = self._auto_intervals()

    def _validate_inputs(
        self,
        conc: NCAConcentration,
        dose: Optional[NCADose],
    ) -> None:
        """Validate inputs."""
        if not isinstance(conc, NCAConcentration):
            raise TypeError("conc must be an NCAConcentration object")

        if dose is not None and not isinstance(dose, NCADose):
            raise TypeError("dose must be an NCADose object")

        # Check subject alignment if both have subjects
        if dose is not None:
            if conc.subject_col and dose.subject_col:
                conc_subjects = set(conc.subjects)
                dose_subjects = set(dose.subjects)
                if not conc_subjects.issubset(dose_subjects):
                    missing = conc_subjects - dose_subjects
                    raise ValueError(
                        f"Some concentration subjects not found in dose data: {missing}"
                    )

    def _auto_intervals(self) -> List["Interval"]:
        """Automatically determine calculation intervals."""
        from pynca.analysis.intervals import Interval

        intervals = []

        for subject in self.conc.subjects:
            conc_array, time_array = self.conc.get_conc_time(subject)

            if len(time_array) == 0:
                continue

            start = time_array.min()
            end = time_array.max()

            # Check for dosing information
            if self.dose is not None:
                dose_time, _ = self.dose.get_first_dose(subject)
                if dose_time is not None:
                    start = min(start, dose_time)

            intervals.append(
                Interval(
                    start=start,
                    end=end,
                    subject=subject,
                )
            )

        return intervals

    @property
    def intervals(self) -> List["Interval"]:
        """Get calculation intervals."""
        return self._intervals

    def set_intervals(
        self,
        intervals: List["Interval"],
    ) -> "NCAData":
        """
        Set calculation intervals.

        Parameters
        ----------
        intervals : list of Interval
            New intervals

        Returns
        -------
        NCAData
            Self for chaining
        """
        self._intervals = intervals
        return self

    def add_interval(
        self,
        start: float,
        end: float,
        parameters: Optional[List[str]] = None,
        subject: Optional[Union[str, int, List]] = None,
        name: Optional[str] = None,
    ) -> "NCAData":
        """
        Add a calculation interval.

        Parameters
        ----------
        start : float
            Interval start time
        end : float
            Interval end time
        parameters : list of str, optional
            Parameters to calculate
        subject : str, int, or list, optional
            Subject(s) to apply interval to. None = all subjects
        name : str, optional
            Interval name

        Returns
        -------
        NCAData
            Self for chaining
        """
        from pynca.analysis.intervals import Interval

        subjects = subject
        if subjects is None:
            subjects = self.conc.subjects
        elif not isinstance(subjects, list):
            subjects = [subjects]

        for subj in subjects:
            self._intervals.append(
                Interval(
                    start=start,
                    end=end,
                    parameters=parameters,
                    subject=subj,
                    name=name,
                )
            )

        return self

    def get_subject_data(
        self,
        subject: Optional[Union[str, int]] = None,
    ) -> Tuple[np.ndarray, np.ndarray, Optional[float], Optional[float]]:
        """
        Get concentration, time, and dose data for a subject.

        Parameters
        ----------
        subject : str or int, optional
            Subject identifier

        Returns
        -------
        tuple
            (conc_array, time_array, dose_amount, dose_time)
        """
        conc_array, time_array = self.conc.get_conc_time(subject)

        dose_amount = None
        dose_time = None
        if self.dose is not None:
            dose_time, dose_amount = self.dose.get_first_dose(subject)

        return conc_array, time_array, dose_amount, dose_time

    def get_intervals_for_subject(
        self,
        subject: Optional[Union[str, int]] = None,
    ) -> List["Interval"]:
        """
        Get intervals applicable to a subject.

        Parameters
        ----------
        subject : str or int, optional
            Subject identifier

        Returns
        -------
        list of Interval
            Applicable intervals
        """
        return [
            i for i in self._intervals
            if i.subject is None or i.subject == subject
        ]

    def iter_analysis_units(self) -> Iterator[Tuple[any, "Interval", Dict]]:
        """
        Iterate over analysis units (subject-interval combinations).

        Yields
        ------
        tuple
            (subject, interval, data_dict)
        """
        for subject in self.conc.subjects:
            conc_array, time_array, dose_amount, dose_time = self.get_subject_data(
                subject
            )

            for interval in self.get_intervals_for_subject(subject):
                yield subject, interval, {
                    "conc": conc_array,
                    "time": time_array,
                    "dose": dose_amount,
                    "dose_time": dose_time,
                    "duration": (
                        self.dose.get_duration(subject)
                        if self.dose is not None
                        else 0.0
                    ),
                    "route": (
                        self.dose.get_route(subject)
                        if self.dose is not None
                        else None
                    ),
                }

    @property
    def subjects(self) -> List:
        """Get list of subjects."""
        return self.conc.subjects

    @property
    def n_subjects(self) -> int:
        """Get number of subjects."""
        return self.conc.n_subjects

    def filter(
        self,
        subjects: Optional[List] = None,
    ) -> "NCAData":
        """
        Filter to specific subjects.

        Parameters
        ----------
        subjects : list, optional
            List of subjects to include

        Returns
        -------
        NCAData
            Filtered data
        """
        new_conc = self.conc.filter(subjects=subjects)
        new_dose = self.dose.filter(subjects=subjects) if self.dose else None

        # Filter intervals
        new_intervals = [
            i for i in self._intervals
            if i.subject is None or i.subject in subjects
        ]

        return NCAData(
            conc=new_conc,
            dose=new_dose,
            intervals=new_intervals,
        )

    def __repr__(self) -> str:
        return (
            f"NCAData(n_subjects={self.n_subjects}, "
            f"n_intervals={len(self._intervals)}, "
            f"has_dose={self.dose is not None})"
        )
