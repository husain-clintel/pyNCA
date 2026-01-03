"""NCA Results class."""

from typing import Any, Dict, List, Optional, Union
import numpy as np
import pandas as pd


class NCAResults:
    """
    Container for NCA calculation results.

    This class stores and manages NCA results, providing methods for
    filtering, summarizing, and exporting results.

    Parameters
    ----------
    results : list of dict
        List of result dictionaries, each containing:
        - subject: subject identifier
        - parameter: parameter name
        - value: calculated value
        - interval_start: interval start time
        - interval_end: interval end time
        - Additional metadata

    Attributes
    ----------
    data : pd.DataFrame
        Results as a DataFrame

    Examples
    --------
    >>> results = nca.run_nca(data)
    >>> results.to_dataframe()
    >>> results.summary()
    """

    def __init__(
        self,
        results: Optional[List[Dict[str, Any]]] = None,
    ):
        self._results = results or []
        self._data = None

    def add_result(
        self,
        subject: Any,
        parameter: str,
        value: float,
        interval_start: float,
        interval_end: float,
        **metadata,
    ) -> "NCAResults":
        """
        Add a single result.

        Parameters
        ----------
        subject : any
            Subject identifier
        parameter : str
            Parameter name
        value : float
            Calculated value
        interval_start : float
            Interval start time
        interval_end : float
            Interval end time
        **metadata
            Additional metadata (units, etc.)

        Returns
        -------
        NCAResults
            Self for chaining
        """
        result = {
            "subject": subject,
            "parameter": parameter,
            "value": value,
            "interval_start": interval_start,
            "interval_end": interval_end,
            **metadata,
        }
        self._results.append(result)
        self._data = None  # Invalidate cached DataFrame
        return self

    def add_results(
        self,
        results: List[Dict[str, Any]],
    ) -> "NCAResults":
        """
        Add multiple results.

        Parameters
        ----------
        results : list of dict
            Results to add

        Returns
        -------
        NCAResults
            Self for chaining
        """
        self._results.extend(results)
        self._data = None
        return self

    @property
    def data(self) -> pd.DataFrame:
        """Get results as DataFrame."""
        if self._data is None:
            self._data = pd.DataFrame(self._results)
        return self._data

    def to_dataframe(
        self,
        wide: bool = False,
    ) -> pd.DataFrame:
        """
        Convert results to DataFrame.

        Parameters
        ----------
        wide : bool
            If True, pivot to wide format (parameters as columns)

        Returns
        -------
        pd.DataFrame
            Results DataFrame
        """
        df = self.data.copy()

        if wide and len(df) > 0:
            # Pivot to wide format
            index_cols = ["subject", "interval_start", "interval_end"]
            index_cols = [c for c in index_cols if c in df.columns]

            df = df.pivot_table(
                index=index_cols,
                columns="parameter",
                values="value",
                aggfunc="first",
            ).reset_index()
            df.columns.name = None

        return df

    def filter(
        self,
        subjects: Optional[List] = None,
        parameters: Optional[List[str]] = None,
        interval_start: Optional[float] = None,
        interval_end: Optional[float] = None,
    ) -> "NCAResults":
        """
        Filter results.

        Parameters
        ----------
        subjects : list, optional
            Subjects to include
        parameters : list of str, optional
            Parameters to include
        interval_start : float, optional
            Filter by interval start
        interval_end : float, optional
            Filter by interval end

        Returns
        -------
        NCAResults
            Filtered results
        """
        filtered = self._results.copy()

        if subjects is not None:
            filtered = [r for r in filtered if r.get("subject") in subjects]

        if parameters is not None:
            filtered = [r for r in filtered if r.get("parameter") in parameters]

        if interval_start is not None:
            filtered = [
                r for r in filtered if r.get("interval_start") == interval_start
            ]

        if interval_end is not None:
            filtered = [r for r in filtered if r.get("interval_end") == interval_end]

        return NCAResults(filtered)

    def get_parameter(
        self,
        parameter: str,
        subject: Optional[Any] = None,
    ) -> Union[float, pd.Series]:
        """
        Get values for a specific parameter.

        Parameters
        ----------
        parameter : str
            Parameter name
        subject : any, optional
            Subject identifier

        Returns
        -------
        float or pd.Series
            Parameter value(s)
        """
        df = self.data
        mask = df["parameter"] == parameter

        if subject is not None:
            mask &= df["subject"] == subject
            values = df.loc[mask, "value"]
            return values.iloc[0] if len(values) > 0 else np.nan

        return df.loc[mask, ["subject", "value"]].set_index("subject")["value"]

    def exclude(
        self,
        subjects: Optional[List] = None,
        parameters: Optional[List[str]] = None,
    ) -> "NCAResults":
        """
        Exclude specific subjects or parameters.

        Parameters
        ----------
        subjects : list, optional
            Subjects to exclude
        parameters : list of str, optional
            Parameters to exclude

        Returns
        -------
        NCAResults
            Results with exclusions
        """
        filtered = self._results.copy()

        if subjects is not None:
            filtered = [r for r in filtered if r.get("subject") not in subjects]

        if parameters is not None:
            filtered = [r for r in filtered if r.get("parameter") not in parameters]

        return NCAResults(filtered)

    def summary(
        self,
        parameters: Optional[List[str]] = None,
        stats: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Summarize results across subjects.

        Parameters
        ----------
        parameters : list of str, optional
            Parameters to summarize (None = all)
        stats : list of str, optional
            Statistics to calculate

        Returns
        -------
        pd.DataFrame
            Summary statistics
        """
        from pynca.summary.summarize import summarize_results

        return summarize_results(self, parameters=parameters, stats=stats)

    def to_csv(
        self,
        path: str,
        wide: bool = True,
        **kwargs,
    ) -> None:
        """
        Export results to CSV.

        Parameters
        ----------
        path : str
            Output file path
        wide : bool
            Use wide format
        **kwargs
            Additional arguments for DataFrame.to_csv
        """
        df = self.to_dataframe(wide=wide)
        df.to_csv(path, index=False, **kwargs)

    def to_excel(
        self,
        path: str,
        wide: bool = True,
        **kwargs,
    ) -> None:
        """
        Export results to Excel.

        Parameters
        ----------
        path : str
            Output file path
        wide : bool
            Use wide format
        **kwargs
            Additional arguments for DataFrame.to_excel
        """
        df = self.to_dataframe(wide=wide)
        df.to_excel(path, index=False, **kwargs)

    @property
    def parameters(self) -> List[str]:
        """Get list of unique parameters."""
        return self.data["parameter"].unique().tolist() if len(self.data) > 0 else []

    @property
    def subjects(self) -> List:
        """Get list of unique subjects."""
        return self.data["subject"].unique().tolist() if len(self.data) > 0 else []

    @property
    def n_subjects(self) -> int:
        """Get number of subjects."""
        return len(self.subjects)

    def __len__(self) -> int:
        """Return number of results."""
        return len(self._results)

    def __repr__(self) -> str:
        n_params = len(self.parameters)
        return (
            f"NCAResults(n_results={len(self)}, n_subjects={self.n_subjects}, "
            f"n_parameters={n_params})"
        )
