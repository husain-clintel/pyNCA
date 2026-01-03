"""NCA Concentration data class."""

from typing import Dict, List, Optional, Union, Iterator, Tuple
import numpy as np
import pandas as pd

from pynca.utils.validation import validate_concentration_data, validate_time_data


class NCAConcentration:
    """
    Container for concentration-time data.

    This class organizes concentration-time data with subject/group identifiers,
    analogous to PKNCAconc in the R package.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing concentration-time data
    conc_col : str
        Name of the concentration column
    time_col : str
        Name of the time column
    subject_col : str, optional
        Name of the subject identifier column
    group_cols : list of str, optional
        Names of additional grouping columns (e.g., treatment, period)
    analyte_col : str, optional
        Name of the analyte column (for multi-analyte data)

    Attributes
    ----------
    data : pd.DataFrame
        The underlying data
    conc_col : str
        Concentration column name
    time_col : str
        Time column name
    subject_col : str
        Subject column name
    group_cols : list
        Grouping column names

    Examples
    --------
    >>> import pandas as pd
    >>> import pynca as nca
    >>> df = pd.DataFrame({
    ...     'subject': [1, 1, 1, 2, 2, 2],
    ...     'time': [0, 1, 2, 0, 1, 2],
    ...     'conc': [0, 10, 5, 0, 12, 6]
    ... })
    >>> conc = nca.NCAConcentration(df, conc_col='conc', time_col='time',
    ...                              subject_col='subject')
    """

    def __init__(
        self,
        data: pd.DataFrame,
        conc_col: str,
        time_col: str,
        subject_col: Optional[str] = None,
        group_cols: Optional[List[str]] = None,
        analyte_col: Optional[str] = None,
    ):
        self._validate_inputs(data, conc_col, time_col, subject_col, group_cols)

        self.data = data.copy()
        self.conc_col = conc_col
        self.time_col = time_col
        self.subject_col = subject_col
        self.group_cols = group_cols or []
        self.analyte_col = analyte_col

        # Sort data by subject (if present) and time
        sort_cols = []
        if self.subject_col:
            sort_cols.append(self.subject_col)
        sort_cols.extend(self.group_cols)
        if self.analyte_col:
            sort_cols.append(self.analyte_col)
        sort_cols.append(self.time_col)

        self.data = self.data.sort_values(sort_cols).reset_index(drop=True)

    def _validate_inputs(
        self,
        data: pd.DataFrame,
        conc_col: str,
        time_col: str,
        subject_col: Optional[str],
        group_cols: Optional[List[str]],
    ) -> None:
        """Validate input parameters."""
        if not isinstance(data, pd.DataFrame):
            raise TypeError("data must be a pandas DataFrame")

        if conc_col not in data.columns:
            raise ValueError(f"Concentration column '{conc_col}' not found in data")

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
        if self.analyte_col:
            cols.append(self.analyte_col)
        return cols

    def get_subject_data(
        self,
        subject: Optional[Union[str, int]] = None,
        **group_filters,
    ) -> pd.DataFrame:
        """
        Get data for a specific subject.

        Parameters
        ----------
        subject : str or int, optional
            Subject identifier
        **group_filters
            Additional filters for group columns

        Returns
        -------
        pd.DataFrame
            Filtered data
        """
        mask = pd.Series(True, index=self.data.index)

        if subject is not None and self.subject_col:
            mask &= self.data[self.subject_col] == subject

        for col, value in group_filters.items():
            if col in self.data.columns:
                mask &= self.data[col] == value

        return self.data[mask].copy()

    def get_conc_time(
        self,
        subject: Optional[Union[str, int]] = None,
        **group_filters,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get concentration and time arrays for a subject.

        Parameters
        ----------
        subject : str or int, optional
            Subject identifier
        **group_filters
            Additional filters

        Returns
        -------
        tuple
            (concentration_array, time_array)
        """
        df = self.get_subject_data(subject, **group_filters)
        return df[self.conc_col].values, df[self.time_col].values

    def iter_subjects(self) -> Iterator[Tuple[any, pd.DataFrame]]:
        """
        Iterate over subjects.

        Yields
        ------
        tuple
            (subject_id, subject_data)
        """
        if self.subject_col is None:
            yield None, self.data
        else:
            for subject in self.subjects:
                yield subject, self.get_subject_data(subject)

    def iter_groups(self) -> Iterator[Tuple[tuple, pd.DataFrame]]:
        """
        Iterate over all grouping combinations.

        Yields
        ------
        tuple
            (group_key, group_data)
        """
        if not self.grouping_columns:
            yield (), self.data
        else:
            for group_key, group_data in self.data.groupby(self.grouping_columns):
                if not isinstance(group_key, tuple):
                    group_key = (group_key,)
                yield group_key, group_data

    def filter(
        self,
        subjects: Optional[List] = None,
        **conditions,
    ) -> "NCAConcentration":
        """
        Filter concentration data.

        Parameters
        ----------
        subjects : list, optional
            List of subjects to include
        **conditions
            Column conditions for filtering

        Returns
        -------
        NCAConcentration
            Filtered concentration object
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

        return NCAConcentration(
            data=self.data[mask],
            conc_col=self.conc_col,
            time_col=self.time_col,
            subject_col=self.subject_col,
            group_cols=self.group_cols,
            analyte_col=self.analyte_col,
        )

    def to_dataframe(self) -> pd.DataFrame:
        """Return the underlying DataFrame."""
        return self.data.copy()

    def __len__(self) -> int:
        """Return number of observations."""
        return len(self.data)

    def __repr__(self) -> str:
        return (
            f"NCAConcentration(n_obs={len(self)}, n_subjects={self.n_subjects}, "
            f"conc_col='{self.conc_col}', time_col='{self.time_col}')"
        )
