"""NCA Dose data class."""

from typing import Dict, List, Optional, Union, Iterator, Tuple
import numpy as np
import pandas as pd


class NCADose:
    """
    Container for dosing information.

    This class organizes dose data with subject/group identifiers,
    analogous to PKNCAdose in the R package.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing dose data
    dose_col : str
        Name of the dose amount column
    time_col : str
        Name of the dose time column
    subject_col : str, optional
        Name of the subject identifier column
    route_col : str, optional
        Name of the route of administration column
    duration_col : str, optional
        Name of the infusion duration column
    group_cols : list of str, optional
        Names of additional grouping columns

    Attributes
    ----------
    data : pd.DataFrame
        The underlying data
    dose_col : str
        Dose column name
    time_col : str
        Time column name
    subject_col : str
        Subject column name

    Examples
    --------
    >>> import pandas as pd
    >>> import pynca as nca
    >>> df = pd.DataFrame({
    ...     'subject': [1, 2],
    ...     'time': [0, 0],
    ...     'dose': [100, 100]
    ... })
    >>> dose = nca.NCADose(df, dose_col='dose', time_col='time',
    ...                     subject_col='subject')
    """

    def __init__(
        self,
        data: pd.DataFrame,
        dose_col: str,
        time_col: str,
        subject_col: Optional[str] = None,
        route_col: Optional[str] = None,
        duration_col: Optional[str] = None,
        group_cols: Optional[List[str]] = None,
    ):
        self._validate_inputs(data, dose_col, time_col, subject_col, group_cols)

        self.data = data.copy()
        self.dose_col = dose_col
        self.time_col = time_col
        self.subject_col = subject_col
        self.route_col = route_col
        self.duration_col = duration_col
        self.group_cols = group_cols or []

        # Sort data
        sort_cols = []
        if self.subject_col:
            sort_cols.append(self.subject_col)
        sort_cols.extend(self.group_cols)
        sort_cols.append(self.time_col)

        self.data = self.data.sort_values(sort_cols).reset_index(drop=True)

    def _validate_inputs(
        self,
        data: pd.DataFrame,
        dose_col: str,
        time_col: str,
        subject_col: Optional[str],
        group_cols: Optional[List[str]],
    ) -> None:
        """Validate input parameters."""
        if not isinstance(data, pd.DataFrame):
            raise TypeError("data must be a pandas DataFrame")

        if dose_col not in data.columns:
            raise ValueError(f"Dose column '{dose_col}' not found in data")

        if time_col not in data.columns:
            raise ValueError(f"Time column '{time_col}' not found in data")

        if subject_col and subject_col not in data.columns:
            raise ValueError(f"Subject column '{subject_col}' not found in data")

        if group_cols:
            for col in group_cols:
                if col not in data.columns:
                    raise ValueError(f"Group column '{col}' not found in data")

    @property
    def subjects(self) -> List:
        """Get list of unique subjects."""
        if self.subject_col is None:
            return [None]
        return self.data[self.subject_col].unique().tolist()

    @property
    def n_subjects(self) -> int:
        """Get number of subjects."""
        return len(self.subjects)

    @property
    def grouping_columns(self) -> List[str]:
        """Get all grouping columns including subject."""
        cols = []
        if self.subject_col:
            cols.append(self.subject_col)
        cols.extend(self.group_cols)
        return cols

    def get_subject_dose(
        self,
        subject: Optional[Union[str, int]] = None,
        **group_filters,
    ) -> pd.DataFrame:
        """
        Get dose data for a specific subject.

        Parameters
        ----------
        subject : str or int, optional
            Subject identifier
        **group_filters
            Additional filters

        Returns
        -------
        pd.DataFrame
            Filtered dose data
        """
        mask = pd.Series(True, index=self.data.index)

        if subject is not None and self.subject_col:
            mask &= self.data[self.subject_col] == subject

        for col, value in group_filters.items():
            if col in self.data.columns:
                mask &= self.data[col] == value

        return self.data[mask].copy()

    def get_first_dose(
        self,
        subject: Optional[Union[str, int]] = None,
        **group_filters,
    ) -> Tuple[float, float]:
        """
        Get the first dose time and amount for a subject.

        Parameters
        ----------
        subject : str or int, optional
            Subject identifier
        **group_filters
            Additional filters

        Returns
        -------
        tuple
            (dose_time, dose_amount)
        """
        df = self.get_subject_dose(subject, **group_filters)
        if len(df) == 0:
            return None, None
        first = df.iloc[0]
        return first[self.time_col], first[self.dose_col]

    def get_total_dose(
        self,
        subject: Optional[Union[str, int]] = None,
        **group_filters,
    ) -> float:
        """
        Get total dose for a subject.

        Parameters
        ----------
        subject : str or int, optional
            Subject identifier
        **group_filters
            Additional filters

        Returns
        -------
        float
            Total dose amount
        """
        df = self.get_subject_dose(subject, **group_filters)
        return df[self.dose_col].sum()

    def get_route(
        self,
        subject: Optional[Union[str, int]] = None,
        **group_filters,
    ) -> Optional[str]:
        """
        Get route of administration for a subject.

        Parameters
        ----------
        subject : str or int, optional
            Subject identifier
        **group_filters
            Additional filters

        Returns
        -------
        str or None
            Route of administration
        """
        if self.route_col is None:
            return None
        df = self.get_subject_dose(subject, **group_filters)
        if len(df) == 0:
            return None
        return df[self.route_col].iloc[0]

    def get_duration(
        self,
        subject: Optional[Union[str, int]] = None,
        **group_filters,
    ) -> Optional[float]:
        """
        Get infusion duration for a subject.

        Parameters
        ----------
        subject : str or int, optional
            Subject identifier
        **group_filters
            Additional filters

        Returns
        -------
        float or None
            Infusion duration
        """
        if self.duration_col is None:
            return 0.0
        df = self.get_subject_dose(subject, **group_filters)
        if len(df) == 0:
            return None
        return df[self.duration_col].iloc[0]

    def is_single_dose(
        self,
        subject: Optional[Union[str, int]] = None,
        **group_filters,
    ) -> bool:
        """
        Check if subject received single dose.

        Parameters
        ----------
        subject : str or int, optional
            Subject identifier
        **group_filters
            Additional filters

        Returns
        -------
        bool
            True if single dose
        """
        df = self.get_subject_dose(subject, **group_filters)
        return len(df) == 1

    def iter_subjects(self) -> Iterator[Tuple[any, pd.DataFrame]]:
        """
        Iterate over subjects.

        Yields
        ------
        tuple
            (subject_id, subject_dose_data)
        """
        if self.subject_col is None:
            yield None, self.data
        else:
            for subject in self.subjects:
                yield subject, self.get_subject_dose(subject)

    def filter(
        self,
        subjects: Optional[List] = None,
        **conditions,
    ) -> "NCADose":
        """
        Filter dose data.

        Parameters
        ----------
        subjects : list, optional
            List of subjects to include
        **conditions
            Column conditions for filtering

        Returns
        -------
        NCADose
            Filtered dose object
        """
        mask = pd.Series(True, index=self.data.index)

        if subjects is not None and self.subject_col:
            mask &= self.data[self.subject_col].isin(subjects)

        for col, value in conditions.items():
            if col in self.data.columns:
                if isinstance(value, (list, tuple)):
                    mask &= self.data[col].isin(value)
                else:
                    mask &= self.data[col] == value

        return NCADose(
            data=self.data[mask],
            dose_col=self.dose_col,
            time_col=self.time_col,
            subject_col=self.subject_col,
            route_col=self.route_col,
            duration_col=self.duration_col,
            group_cols=self.group_cols,
        )

    def to_dataframe(self) -> pd.DataFrame:
        """Return the underlying DataFrame."""
        return self.data.copy()

    def __len__(self) -> int:
        """Return number of dose records."""
        return len(self.data)

    def __repr__(self) -> str:
        return (
            f"NCADose(n_doses={len(self)}, n_subjects={self.n_subjects}, "
            f"dose_col='{self.dose_col}', time_col='{self.time_col}')"
        )
