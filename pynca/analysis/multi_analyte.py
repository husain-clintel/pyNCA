"""Multi-analyte NCA analysis and metabolite ratio calculations."""

from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd

from pynca.core.concentration import NCAConcentration
from pynca.core.dose import NCADose
from pynca.core.data import NCAData
from pynca.core.results import NCAResults


def run_multi_analyte_nca(
    data: NCAData,
    analytes: Optional[List[str]] = None,
    **kwargs,
) -> Dict[str, NCAResults]:
    """
    Run NCA analysis for multiple analytes.

    Parameters
    ----------
    data : NCAData
        NCA data object with analyte_col specified
    analytes : list of str, optional
        Specific analytes to analyze. If None, analyze all.
    **kwargs
        Additional arguments passed to run_nca

    Returns
    -------
    dict
        Dictionary mapping analyte name to NCAResults

    Examples
    --------
    >>> # Data with parent drug and metabolite
    >>> conc = NCAConcentration(df, conc_col='conc', time_col='time',
    ...                         subject_col='subject', analyte_col='analyte')
    >>> data = NCAData(conc=conc, dose=dose)
    >>> results = run_multi_analyte_nca(data)
    >>> parent_results = results['parent']
    >>> metabolite_results = results['metabolite']
    """
    from pynca.analysis.nca import run_nca

    if data.conc.analyte_col is None:
        # Single analyte - just run regular NCA
        return {"default": run_nca(data, **kwargs)}

    # Get list of analytes
    available_analytes = data.conc.data[data.conc.analyte_col].unique().tolist()

    if analytes is None:
        analytes = available_analytes
    else:
        # Validate requested analytes exist
        missing = set(analytes) - set(available_analytes)
        if missing:
            raise ValueError(f"Analytes not found in data: {missing}")

    results = {}

    for analyte in analytes:
        # Filter concentration data for this analyte
        analyte_conc = data.conc.filter(**{data.conc.analyte_col: analyte})

        # Create new NCAData for this analyte
        analyte_data = NCAData(
            conc=analyte_conc,
            dose=data.dose,
            intervals=data.intervals,
        )

        # Run NCA
        results[analyte] = run_nca(analyte_data, **kwargs)

    return results


def calc_metabolite_ratio(
    parent_value: float,
    metabolite_value: float,
    ratio_type: str = "metabolite_to_parent",
) -> float:
    """
    Calculate metabolite-to-parent ratio.

    Parameters
    ----------
    parent_value : float
        Parameter value for parent drug
    metabolite_value : float
        Parameter value for metabolite
    ratio_type : str
        Type of ratio:
        - "metabolite_to_parent": M/P ratio
        - "parent_to_metabolite": P/M ratio
        - "molar": Molar ratio (requires MW adjustment externally)

    Returns
    -------
    float
        Ratio value

    Examples
    --------
    >>> calc_metabolite_ratio(parent_auc=100, metabolite_auc=50)
    0.5
    """
    if np.isnan(parent_value) or np.isnan(metabolite_value):
        return np.nan

    if ratio_type == "metabolite_to_parent":
        if parent_value == 0:
            return np.nan
        return metabolite_value / parent_value
    elif ratio_type == "parent_to_metabolite":
        if metabolite_value == 0:
            return np.nan
        return parent_value / metabolite_value
    elif ratio_type == "molar":
        # Same as M/P for molar (MW adjustment should be done externally)
        if parent_value == 0:
            return np.nan
        return metabolite_value / parent_value
    else:
        raise ValueError(f"Unknown ratio_type: {ratio_type}")


def calc_metabolite_ratios(
    parent_results: NCAResults,
    metabolite_results: NCAResults,
    parameters: Optional[List[str]] = None,
    ratio_type: str = "metabolite_to_parent",
) -> pd.DataFrame:
    """
    Calculate metabolite-to-parent ratios for multiple parameters.

    Parameters
    ----------
    parent_results : NCAResults
        NCA results for parent drug
    metabolite_results : NCAResults
        NCA results for metabolite
    parameters : list of str, optional
        Parameters to calculate ratios for. Default: AUC and Cmax
    ratio_type : str
        Type of ratio calculation

    Returns
    -------
    pd.DataFrame
        DataFrame with ratio values for each subject and parameter

    Examples
    --------
    >>> ratios = calc_metabolite_ratios(parent_results, metabolite_results)
    >>> print(ratios[['subject', 'parameter', 'ratio']])
    """
    if parameters is None:
        parameters = ["cmax", "auc.last", "auc.inf.obs"]

    # Get results as DataFrames
    parent_df = parent_results.to_dataframe()
    metabolite_df = metabolite_results.to_dataframe()

    ratio_rows = []

    for param in parameters:
        # Get parameter values for each subject
        parent_param = parent_df[parent_df["parameter"] == param]
        metabolite_param = metabolite_df[metabolite_df["parameter"] == param]

        if len(parent_param) == 0 or len(metabolite_param) == 0:
            continue

        # Merge on subject
        merged = parent_param.merge(
            metabolite_param,
            on="subject",
            suffixes=("_parent", "_metabolite"),
        )

        for _, row in merged.iterrows():
            ratio = calc_metabolite_ratio(
                row["value_parent"],
                row["value_metabolite"],
                ratio_type=ratio_type,
            )

            ratio_rows.append({
                "subject": row["subject"],
                "parameter": param,
                "parent_value": row["value_parent"],
                "metabolite_value": row["value_metabolite"],
                "ratio": ratio,
                "ratio_type": ratio_type,
            })

    return pd.DataFrame(ratio_rows)


def summarize_metabolite_ratios(
    ratios_df: pd.DataFrame,
    stats: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Summarize metabolite ratios across subjects.

    Parameters
    ----------
    ratios_df : pd.DataFrame
        DataFrame from calc_metabolite_ratios
    stats : list of str, optional
        Statistics to calculate

    Returns
    -------
    pd.DataFrame
        Summary statistics for each parameter ratio
    """
    from pynca.summary.summarize import calculate_stats

    if stats is None:
        stats = ["n", "mean", "sd", "cv", "median", "geomean", "geocv"]

    summary_rows = []

    for param in ratios_df["parameter"].unique():
        param_ratios = ratios_df[ratios_df["parameter"] == param]["ratio"].values
        param_ratios = param_ratios[~np.isnan(param_ratios)]

        if len(param_ratios) == 0:
            continue

        row = {"parameter": f"{param}_ratio"}
        row.update(calculate_stats(param_ratios, stats))
        summary_rows.append(row)

    return pd.DataFrame(summary_rows)


def molar_ratio(
    mass_ratio: float,
    parent_mw: float,
    metabolite_mw: float,
) -> float:
    """
    Convert mass-based ratio to molar ratio.

    Parameters
    ----------
    mass_ratio : float
        Mass-based metabolite/parent ratio
    parent_mw : float
        Molecular weight of parent drug
    metabolite_mw : float
        Molecular weight of metabolite

    Returns
    -------
    float
        Molar ratio

    Notes
    -----
    Molar ratio = Mass ratio * (Parent MW / Metabolite MW)
    """
    if np.isnan(mass_ratio) or parent_mw <= 0 or metabolite_mw <= 0:
        return np.nan

    return mass_ratio * (parent_mw / metabolite_mw)


class MultiAnalyteResults:
    """
    Container for multi-analyte NCA results.

    Parameters
    ----------
    results : dict
        Dictionary mapping analyte names to NCAResults objects

    Attributes
    ----------
    analytes : list
        List of analyte names
    results : dict
        Results by analyte
    """

    def __init__(self, results: Dict[str, NCAResults]):
        self.results = results
        self.analytes = list(results.keys())

    def get_analyte(self, analyte: str) -> NCAResults:
        """Get results for a specific analyte."""
        if analyte not in self.results:
            raise KeyError(f"Analyte '{analyte}' not found. Available: {self.analytes}")
        return self.results[analyte]

    def to_dataframe(self, wide: bool = False) -> pd.DataFrame:
        """
        Combine all analyte results into a single DataFrame.

        Parameters
        ----------
        wide : bool
            If True, return wide format

        Returns
        -------
        pd.DataFrame
            Combined results with analyte column
        """
        dfs = []
        for analyte, result in self.results.items():
            df = result.to_dataframe(wide=wide)
            df["analyte"] = analyte
            dfs.append(df)

        if not dfs:
            return pd.DataFrame()

        return pd.concat(dfs, ignore_index=True)

    def calc_ratios(
        self,
        parent_analyte: str,
        metabolite_analyte: str,
        parameters: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Calculate metabolite/parent ratios.

        Parameters
        ----------
        parent_analyte : str
            Name of parent analyte
        metabolite_analyte : str
            Name of metabolite analyte
        parameters : list of str, optional
            Parameters to calculate ratios for

        Returns
        -------
        pd.DataFrame
            Ratio results
        """
        return calc_metabolite_ratios(
            self.get_analyte(parent_analyte),
            self.get_analyte(metabolite_analyte),
            parameters=parameters,
        )

    def __repr__(self) -> str:
        return f"MultiAnalyteResults(analytes={self.analytes})"
