"""Result summarization functions."""

from typing import Any, Callable, Dict, List, Optional, Union, Tuple
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats


def summarize_results(
    results: "NCAResults",
    parameters: Optional[List[str]] = None,
    stats: Optional[List[str]] = None,
    by: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Summarize NCA results across subjects.

    Parameters
    ----------
    results : NCAResults
        NCA calculation results
    parameters : list of str, optional
        Parameters to summarize (None = all)
    stats : list of str, optional
        Statistics to calculate. Options:
        - "n": Count
        - "mean": Arithmetic mean
        - "sd": Standard deviation
        - "cv": Coefficient of variation (%)
        - "median": Median
        - "min": Minimum
        - "max": Maximum
        - "geomean": Geometric mean
        - "geocv": Geometric CV (%)
        - "se": Standard error
        - "ci90": 90% confidence interval
        - "ci95": 95% confidence interval
    by : list of str, optional
        Grouping columns

    Returns
    -------
    pd.DataFrame
        Summary statistics
    """
    if stats is None:
        stats = ["n", "mean", "sd", "cv", "median", "min", "max"]

    df = results.data

    if len(df) == 0:
        return pd.DataFrame()

    # Filter parameters if specified
    if parameters is not None:
        df = df[df["parameter"].isin(parameters)]

    # Calculate statistics for each parameter
    summary_rows = []

    for param in df["parameter"].unique():
        param_df = df[df["parameter"] == param]
        values = param_df["value"].dropna().values

        row = {"parameter": param}
        row.update(calculate_stats(values, stats))
        summary_rows.append(row)

    return pd.DataFrame(summary_rows)


def calculate_stats(
    values: np.ndarray,
    stats: List[str],
) -> Dict[str, float]:
    """
    Calculate summary statistics for a set of values.

    Parameters
    ----------
    values : array-like
        Values to summarize
    stats : list of str
        Statistics to calculate

    Returns
    -------
    dict
        Dictionary of statistic names and values
    """
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]

    result = {}

    for stat in stats:
        if stat == "n":
            result["n"] = len(values)
        elif stat == "mean":
            result["mean"] = np.mean(values) if len(values) > 0 else np.nan
        elif stat == "sd":
            result["sd"] = np.std(values, ddof=1) if len(values) > 1 else np.nan
        elif stat == "cv":
            if len(values) > 1 and np.mean(values) != 0:
                result["cv"] = (np.std(values, ddof=1) / np.mean(values)) * 100
            else:
                result["cv"] = np.nan
        elif stat == "median":
            result["median"] = np.median(values) if len(values) > 0 else np.nan
        elif stat == "min":
            result["min"] = np.min(values) if len(values) > 0 else np.nan
        elif stat == "max":
            result["max"] = np.max(values) if len(values) > 0 else np.nan
        elif stat == "geomean":
            result["geomean"] = geometric_mean(values)
        elif stat == "geocv":
            result["geocv"] = geometric_cv(values)
        elif stat == "se":
            if len(values) > 1:
                result["se"] = np.std(values, ddof=1) / np.sqrt(len(values))
            else:
                result["se"] = np.nan
        elif stat == "ci90":
            ci = confidence_interval(values, 0.90)
            result["ci90_lower"] = ci[0]
            result["ci90_upper"] = ci[1]
        elif stat == "ci95":
            ci = confidence_interval(values, 0.95)
            result["ci95_lower"] = ci[0]
            result["ci95_upper"] = ci[1]
        else:
            result[stat] = np.nan

    return result


def geometric_mean(values: np.ndarray) -> float:
    """
    Calculate geometric mean.

    Parameters
    ----------
    values : array-like
        Values (must be positive)

    Returns
    -------
    float
        Geometric mean
    """
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]
    values = values[values > 0]  # Geometric mean only for positive values

    if len(values) == 0:
        return np.nan

    return np.exp(np.mean(np.log(values)))


def geometric_cv(values: np.ndarray) -> float:
    """
    Calculate geometric coefficient of variation (%).

    GeoCV = sqrt(exp(var(log(x))) - 1) * 100

    Parameters
    ----------
    values : array-like
        Values (must be positive)

    Returns
    -------
    float
        Geometric CV as percentage
    """
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]
    values = values[values > 0]

    if len(values) < 2:
        return np.nan

    log_values = np.log(values)
    var_log = np.var(log_values, ddof=1)

    return np.sqrt(np.exp(var_log) - 1) * 100


def confidence_interval(
    values: np.ndarray,
    confidence: float = 0.95,
) -> tuple:
    """
    Calculate confidence interval for the mean.

    Parameters
    ----------
    values : array-like
        Values
    confidence : float
        Confidence level (0.90, 0.95, etc.)

    Returns
    -------
    tuple
        (lower_bound, upper_bound)
    """
    from scipy import stats

    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]

    if len(values) < 2:
        return (np.nan, np.nan)

    mean = np.mean(values)
    se = np.std(values, ddof=1) / np.sqrt(len(values))
    df = len(values) - 1

    t_crit = stats.t.ppf((1 + confidence) / 2, df)

    return (mean - t_crit * se, mean + t_crit * se)


def geometric_confidence_interval(
    values: np.ndarray,
    confidence: float = 0.90,
) -> tuple:
    """
    Calculate confidence interval for geometric mean.

    Parameters
    ----------
    values : array-like
        Values (must be positive)
    confidence : float
        Confidence level

    Returns
    -------
    tuple
        (lower_bound, upper_bound)
    """
    from scipy import stats

    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]
    values = values[values > 0]

    if len(values) < 2:
        return (np.nan, np.nan)

    log_values = np.log(values)
    mean_log = np.mean(log_values)
    se_log = np.std(log_values, ddof=1) / np.sqrt(len(log_values))
    df = len(values) - 1

    t_crit = stats.t.ppf((1 + confidence) / 2, df)

    lower = np.exp(mean_log - t_crit * se_log)
    upper = np.exp(mean_log + t_crit * se_log)

    return (lower, upper)


def format_summary(
    summary_df: pd.DataFrame,
    decimal_places: int = 3,
) -> pd.DataFrame:
    """
    Format summary table for display/export.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary statistics
    decimal_places : int
        Number of decimal places

    Returns
    -------
    pd.DataFrame
        Formatted summary
    """
    formatted = summary_df.copy()

    # Round numeric columns
    numeric_cols = formatted.select_dtypes(include=[np.number]).columns
    formatted[numeric_cols] = formatted[numeric_cols].round(decimal_places)

    return formatted


def pivot_summary(
    summary_df: pd.DataFrame,
    stat: str = "mean",
) -> pd.DataFrame:
    """
    Pivot summary to wide format.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary statistics
    stat : str
        Statistic to show in pivoted table

    Returns
    -------
    pd.DataFrame
        Pivoted summary
    """
    if stat not in summary_df.columns:
        raise ValueError(f"Statistic '{stat}' not found in summary")

    return summary_df.pivot_table(
        values=stat,
        index="parameter",
        aggfunc="first",
    )


def bootstrap_ci(
    values: np.ndarray,
    statistic: str = "mean",
    confidence: float = 0.95,
    n_bootstrap: int = 1000,
    method: str = "percentile",
    random_state: Optional[int] = None,
) -> Dict[str, float]:
    """
    Calculate bootstrap confidence interval.

    Parameters
    ----------
    values : array-like
        Values to bootstrap
    statistic : str
        Statistic to calculate: "mean", "median", "geomean"
    confidence : float
        Confidence level (default 0.95)
    n_bootstrap : int
        Number of bootstrap samples (default 1000)
    method : str
        Bootstrap CI method: "percentile", "bca" (bias-corrected accelerated)
    random_state : int, optional
        Random seed for reproducibility

    Returns
    -------
    dict
        Dictionary with 'estimate', 'lower', 'upper', 'se'

    Examples
    --------
    >>> values = np.array([10, 12, 15, 11, 14, 13])
    >>> result = bootstrap_ci(values, statistic="mean", n_bootstrap=2000)
    >>> print(f"Mean: {result['estimate']:.2f} ({result['lower']:.2f}, {result['upper']:.2f})")
    """
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]

    if len(values) < 2:
        return {
            "estimate": np.nan,
            "lower": np.nan,
            "upper": np.nan,
            "se": np.nan,
        }

    if random_state is not None:
        np.random.seed(random_state)

    n = len(values)

    # Define statistic function
    if statistic == "mean":
        stat_func = np.mean
    elif statistic == "median":
        stat_func = np.median
    elif statistic == "geomean":
        def stat_func(x):
            x = x[x > 0]
            if len(x) == 0:
                return np.nan
            return np.exp(np.mean(np.log(x)))
    else:
        raise ValueError(f"Unknown statistic: {statistic}")

    # Original estimate
    original_estimate = stat_func(values)

    # Generate bootstrap samples
    bootstrap_stats = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        sample = np.random.choice(values, size=n, replace=True)
        bootstrap_stats[i] = stat_func(sample)

    # Remove NaN values
    bootstrap_stats = bootstrap_stats[~np.isnan(bootstrap_stats)]

    if len(bootstrap_stats) < 10:
        return {
            "estimate": original_estimate,
            "lower": np.nan,
            "upper": np.nan,
            "se": np.nan,
        }

    # Calculate CI based on method
    alpha = 1 - confidence

    if method == "percentile":
        lower = np.percentile(bootstrap_stats, alpha / 2 * 100)
        upper = np.percentile(bootstrap_stats, (1 - alpha / 2) * 100)
    elif method == "bca":
        # Bias-corrected and accelerated (BCa) method
        lower, upper = _bca_ci(values, bootstrap_stats, original_estimate,
                               stat_func, confidence)
    else:
        raise ValueError(f"Unknown method: {method}")

    return {
        "estimate": original_estimate,
        "lower": lower,
        "upper": upper,
        "se": np.std(bootstrap_stats),
    }


def _bca_ci(
    values: np.ndarray,
    bootstrap_stats: np.ndarray,
    original_estimate: float,
    stat_func: Callable,
    confidence: float,
) -> tuple:
    """
    Calculate BCa (bias-corrected and accelerated) confidence interval.

    Parameters
    ----------
    values : array-like
        Original values
    bootstrap_stats : array-like
        Bootstrap statistics
    original_estimate : float
        Statistic from original data
    stat_func : callable
        Function to calculate statistic
    confidence : float
        Confidence level

    Returns
    -------
    tuple
        (lower, upper) bounds
    """
    from scipy import stats as scipy_stats

    n = len(values)
    alpha = 1 - confidence

    # Bias correction factor (z0)
    prop_less = np.mean(bootstrap_stats < original_estimate)
    if prop_less == 0:
        prop_less = 1 / (2 * len(bootstrap_stats))
    elif prop_less == 1:
        prop_less = 1 - 1 / (2 * len(bootstrap_stats))
    z0 = scipy_stats.norm.ppf(prop_less)

    # Acceleration factor (a) using jackknife
    jackknife_stats = np.zeros(n)
    for i in range(n):
        jack_sample = np.delete(values, i)
        jackknife_stats[i] = stat_func(jack_sample)

    jack_mean = np.mean(jackknife_stats)
    numerator = np.sum((jack_mean - jackknife_stats) ** 3)
    denominator = 6 * (np.sum((jack_mean - jackknife_stats) ** 2) ** 1.5)

    if denominator == 0:
        a = 0
    else:
        a = numerator / denominator

    # Adjusted percentiles
    z_alpha_lower = scipy_stats.norm.ppf(alpha / 2)
    z_alpha_upper = scipy_stats.norm.ppf(1 - alpha / 2)

    # BCa adjustment
    def bca_percentile(z_alpha):
        num = z0 + z_alpha
        denom = 1 - a * num
        if denom == 0:
            return 0.5
        adjusted_z = z0 + num / denom
        return scipy_stats.norm.cdf(adjusted_z)

    p_lower = bca_percentile(z_alpha_lower)
    p_upper = bca_percentile(z_alpha_upper)

    # Clip to valid range
    p_lower = np.clip(p_lower, 0.001, 0.999)
    p_upper = np.clip(p_upper, 0.001, 0.999)

    lower = np.percentile(bootstrap_stats, p_lower * 100)
    upper = np.percentile(bootstrap_stats, p_upper * 100)

    return lower, upper


def bootstrap_summary(
    results: "NCAResults",
    parameters: Optional[List[str]] = None,
    statistic: str = "mean",
    confidence: float = 0.95,
    n_bootstrap: int = 1000,
    method: str = "percentile",
    random_state: Optional[int] = None,
) -> pd.DataFrame:
    """
    Calculate bootstrap confidence intervals for NCA parameters.

    Parameters
    ----------
    results : NCAResults
        NCA calculation results
    parameters : list of str, optional
        Parameters to summarize (None = all)
    statistic : str
        Statistic to bootstrap: "mean", "median", "geomean"
    confidence : float
        Confidence level (default 0.95)
    n_bootstrap : int
        Number of bootstrap samples
    method : str
        Bootstrap CI method: "percentile", "bca"
    random_state : int, optional
        Random seed for reproducibility

    Returns
    -------
    pd.DataFrame
        Bootstrap summary with estimate, lower, upper, and SE

    Examples
    --------
    >>> summary = bootstrap_summary(results, parameters=["cmax", "auc.last"])
    >>> print(summary)
    """
    df = results.data

    if len(df) == 0:
        return pd.DataFrame()

    # Filter parameters if specified
    if parameters is not None:
        df = df[df["parameter"].isin(parameters)]

    summary_rows = []

    for param in df["parameter"].unique():
        param_df = df[df["parameter"] == param]
        values = param_df["value"].dropna().values

        boot_result = bootstrap_ci(
            values,
            statistic=statistic,
            confidence=confidence,
            n_bootstrap=n_bootstrap,
            method=method,
            random_state=random_state,
        )

        summary_rows.append({
            "parameter": param,
            "n": len(values),
            "statistic": statistic,
            "estimate": boot_result["estimate"],
            "lower": boot_result["lower"],
            "upper": boot_result["upper"],
            "se": boot_result["se"],
            "confidence": confidence,
        })

    return pd.DataFrame(summary_rows)


# =============================================================================
# Population NCA Summary Functions
# =============================================================================


def summarize_by_group(
    results: "NCAResults",
    group_col: str,
    parameters: Optional[List[str]] = None,
    stats: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Summarize NCA results by group (e.g., treatment, dose level).

    Parameters
    ----------
    results : NCAResults
        NCA calculation results
    group_col : str
        Column name for grouping variable
    parameters : list of str, optional
        Parameters to summarize (None = all)
    stats : list of str, optional
        Statistics to calculate

    Returns
    -------
    pd.DataFrame
        Summary statistics by group

    Examples
    --------
    >>> summary = summarize_by_group(results, group_col='treatment',
    ...                              parameters=['cmax', 'auc.last'])
    """
    if stats is None:
        stats = ["n", "mean", "sd", "cv", "geomean", "geocv"]

    df = results.data

    if len(df) == 0:
        return pd.DataFrame()

    if group_col not in df.columns:
        raise ValueError(f"Group column '{group_col}' not found in results")

    # Filter parameters if specified
    if parameters is not None:
        df = df[df["parameter"].isin(parameters)]

    summary_rows = []

    for group in df[group_col].unique():
        group_df = df[df[group_col] == group]

        for param in group_df["parameter"].unique():
            param_df = group_df[group_df["parameter"] == param]
            values = param_df["value"].dropna().values

            row = {
                group_col: group,
                "parameter": param,
            }
            row.update(calculate_stats(values, stats))
            summary_rows.append(row)

    return pd.DataFrame(summary_rows)


def compare_groups(
    results: "NCAResults",
    group_col: str,
    test_group: str,
    ref_group: str,
    parameters: Optional[List[str]] = None,
    confidence: float = 0.90,
    log_transform: bool = True,
) -> pd.DataFrame:
    """
    Compare two groups with ratio and confidence interval.

    Parameters
    ----------
    results : NCAResults
        NCA calculation results
    group_col : str
        Column name for grouping variable
    test_group : str
        Name of test group
    ref_group : str
        Name of reference group
    parameters : list of str, optional
        Parameters to compare (None = all)
    confidence : float
        Confidence level (default 0.90 for bioequivalence)
    log_transform : bool
        If True, calculate geometric mean ratio (default True)

    Returns
    -------
    pd.DataFrame
        Comparison results with ratio and CI

    Examples
    --------
    >>> comparison = compare_groups(results, group_col='treatment',
    ...                             test_group='test', ref_group='reference',
    ...                             parameters=['cmax', 'auc.last'])
    """
    df = results.data

    if len(df) == 0:
        return pd.DataFrame()

    if group_col not in df.columns:
        raise ValueError(f"Group column '{group_col}' not found in results")

    # Filter parameters if specified
    if parameters is not None:
        df = df[df["parameter"].isin(parameters)]

    comparison_rows = []

    for param in df["parameter"].unique():
        param_df = df[df["parameter"] == param]

        test_values = param_df[param_df[group_col] == test_group]["value"].dropna().values
        ref_values = param_df[param_df[group_col] == ref_group]["value"].dropna().values

        if len(test_values) == 0 or len(ref_values) == 0:
            continue

        if log_transform:
            # Geometric mean ratio
            test_values_pos = test_values[test_values > 0]
            ref_values_pos = ref_values[ref_values > 0]

            if len(test_values_pos) < 2 or len(ref_values_pos) < 2:
                continue

            test_log = np.log(test_values_pos)
            ref_log = np.log(ref_values_pos)

            test_geomean = np.exp(np.mean(test_log))
            ref_geomean = np.exp(np.mean(ref_log))
            ratio = test_geomean / ref_geomean

            # CI on log scale, then back-transform
            diff_mean = np.mean(test_log) - np.mean(ref_log)
            se_diff = np.sqrt(
                np.var(test_log, ddof=1) / len(test_log) +
                np.var(ref_log, ddof=1) / len(ref_log)
            )

            # Degrees of freedom (Welch-Satterthwaite approximation)
            v1 = np.var(test_log, ddof=1) / len(test_log)
            v2 = np.var(ref_log, ddof=1) / len(ref_log)
            df_welch = (v1 + v2) ** 2 / (
                v1 ** 2 / (len(test_log) - 1) + v2 ** 2 / (len(ref_log) - 1)
            )

            t_crit = scipy_stats.t.ppf((1 + confidence) / 2, df_welch)
            lower = np.exp(diff_mean - t_crit * se_diff)
            upper = np.exp(diff_mean + t_crit * se_diff)

            comparison_rows.append({
                "parameter": param,
                "test_group": test_group,
                "ref_group": ref_group,
                "test_n": len(test_values_pos),
                "ref_n": len(ref_values_pos),
                "test_geomean": test_geomean,
                "ref_geomean": ref_geomean,
                "ratio": ratio,
                "ratio_pct": ratio * 100,
                "lower": lower,
                "upper": upper,
                "lower_pct": lower * 100,
                "upper_pct": upper * 100,
                "confidence": confidence,
            })
        else:
            # Arithmetic mean ratio
            test_mean = np.mean(test_values)
            ref_mean = np.mean(ref_values)
            ratio = test_mean / ref_mean if ref_mean != 0 else np.nan

            comparison_rows.append({
                "parameter": param,
                "test_group": test_group,
                "ref_group": ref_group,
                "test_n": len(test_values),
                "ref_n": len(ref_values),
                "test_mean": test_mean,
                "ref_mean": ref_mean,
                "ratio": ratio,
                "ratio_pct": ratio * 100 if not np.isnan(ratio) else np.nan,
            })

    return pd.DataFrame(comparison_rows)


def bioequivalence_analysis(
    results: "NCAResults",
    group_col: str,
    test_group: str,
    ref_group: str,
    parameters: Optional[List[str]] = None,
    be_limits: Tuple[float, float] = (0.80, 1.25),
    confidence: float = 0.90,
) -> pd.DataFrame:
    """
    Perform bioequivalence analysis with 90% CI.

    Parameters
    ----------
    results : NCAResults
        NCA calculation results
    group_col : str
        Column name for grouping variable
    test_group : str
        Name of test formulation group
    ref_group : str
        Name of reference formulation group
    parameters : list of str, optional
        Parameters to analyze (default: cmax, auc.last, auc.inf.obs)
    be_limits : tuple
        Bioequivalence limits (default: 0.80-1.25)
    confidence : float
        Confidence level (default 0.90)

    Returns
    -------
    pd.DataFrame
        Bioequivalence results with pass/fail assessment

    Examples
    --------
    >>> be_results = bioequivalence_analysis(
    ...     results,
    ...     group_col='formulation',
    ...     test_group='test',
    ...     ref_group='reference'
    ... )
    """
    if parameters is None:
        parameters = ["cmax", "auc.last", "auc.inf.obs"]

    comparison = compare_groups(
        results,
        group_col=group_col,
        test_group=test_group,
        ref_group=ref_group,
        parameters=parameters,
        confidence=confidence,
        log_transform=True,
    )

    if len(comparison) == 0:
        return pd.DataFrame()

    lower_limit, upper_limit = be_limits

    # Add BE assessment
    comparison["be_lower_limit"] = lower_limit
    comparison["be_upper_limit"] = upper_limit
    comparison["within_limits"] = (
        (comparison["lower"] >= lower_limit) &
        (comparison["upper"] <= upper_limit)
    )
    comparison["be_conclusion"] = comparison["within_limits"].apply(
        lambda x: "Bioequivalent" if x else "Not bioequivalent"
    )

    return comparison


def population_pk_summary(
    results: "NCAResults",
    parameters: Optional[List[str]] = None,
    include_individual: bool = False,
) -> pd.DataFrame:
    """
    Generate comprehensive population PK summary.

    Parameters
    ----------
    results : NCAResults
        NCA calculation results
    parameters : list of str, optional
        Parameters to include (None = all)
    include_individual : bool
        Include individual subject values

    Returns
    -------
    pd.DataFrame
        Population summary table

    Examples
    --------
    >>> pop_summary = population_pk_summary(results)
    """
    stats = [
        "n", "mean", "sd", "cv",
        "median", "min", "max",
        "geomean", "geocv",
    ]

    summary = summarize_results(results, parameters=parameters, stats=stats)

    if len(summary) == 0:
        return pd.DataFrame()

    # Add 95% CI for mean
    df = results.data

    if parameters is not None:
        df = df[df["parameter"].isin(parameters)]

    ci_data = []
    for param in summary["parameter"]:
        values = df[df["parameter"] == param]["value"].dropna().values
        ci = confidence_interval(values, 0.95)
        geo_ci = geometric_confidence_interval(values, 0.90)
        ci_data.append({
            "parameter": param,
            "ci95_lower": ci[0],
            "ci95_upper": ci[1],
            "geo_ci90_lower": geo_ci[0],
            "geo_ci90_upper": geo_ci[1],
        })

    ci_df = pd.DataFrame(ci_data)
    summary = summary.merge(ci_df, on="parameter")

    return summary


def inter_subject_variability(
    results: "NCAResults",
    parameters: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Calculate inter-subject variability metrics.

    Parameters
    ----------
    results : NCAResults
        NCA calculation results
    parameters : list of str, optional
        Parameters to analyze

    Returns
    -------
    pd.DataFrame
        Variability metrics (CV%, GeoCV%, range ratio)

    Examples
    --------
    >>> isv = inter_subject_variability(results)
    """
    df = results.data

    if len(df) == 0:
        return pd.DataFrame()

    if parameters is not None:
        df = df[df["parameter"].isin(parameters)]

    variability_rows = []

    for param in df["parameter"].unique():
        values = df[df["parameter"] == param]["value"].dropna().values

        if len(values) < 2:
            continue

        cv = (np.std(values, ddof=1) / np.mean(values)) * 100 if np.mean(values) != 0 else np.nan
        geocv = geometric_cv(values)
        range_ratio = np.max(values) / np.min(values) if np.min(values) != 0 else np.nan

        variability_rows.append({
            "parameter": param,
            "n": len(values),
            "cv_pct": cv,
            "geocv_pct": geocv,
            "min": np.min(values),
            "max": np.max(values),
            "range_ratio": range_ratio,
            "variability_category": _categorize_variability(geocv),
        })

    return pd.DataFrame(variability_rows)


def _categorize_variability(geocv: float) -> str:
    """Categorize variability based on GeoCV%."""
    if np.isnan(geocv):
        return "Unknown"
    elif geocv <= 30:
        return "Low"
    elif geocv <= 60:
        return "Moderate"
    else:
        return "High"


def dose_proportionality(
    results_by_dose: Dict[float, "NCAResults"],
    parameter: str = "auc.last",
    method: str = "power_model",
) -> Dict[str, Any]:
    """
    Assess dose proportionality across multiple dose levels.

    Parameters
    ----------
    results_by_dose : dict
        Dictionary mapping dose levels to NCAResults
    parameter : str
        Parameter to assess (default: auc.last)
    method : str
        Analysis method: "power_model" or "anova"

    Returns
    -------
    dict
        Dose proportionality assessment results

    Notes
    -----
    Power model: PK = A * Dose^beta
    - beta = 1 indicates dose proportionality
    - 90% CI for beta containing 1 supports proportionality

    Examples
    --------
    >>> results_10mg = run_nca(data_10mg)
    >>> results_20mg = run_nca(data_20mg)
    >>> results_40mg = run_nca(data_40mg)
    >>> dp = dose_proportionality(
    ...     {10: results_10mg, 20: results_20mg, 40: results_40mg},
    ...     parameter='auc.last'
    ... )
    """
    doses = []
    values = []

    for dose, results in results_by_dose.items():
        df = results.data
        param_values = df[df["parameter"] == parameter]["value"].dropna().values

        for val in param_values:
            doses.append(dose)
            values.append(val)

    doses = np.array(doses)
    values = np.array(values)

    if len(doses) < 3:
        return {
            "method": method,
            "parameter": parameter,
            "error": "Insufficient data (need at least 3 dose-value pairs)",
        }

    if method == "power_model":
        # Log-linear regression: log(PK) = log(A) + beta * log(Dose)
        log_doses = np.log(doses)
        log_values = np.log(values)

        # Filter out non-positive values
        valid = (values > 0) & (doses > 0)
        log_doses = log_doses[valid]
        log_values = log_values[valid]

        if len(log_doses) < 3:
            return {
                "method": method,
                "parameter": parameter,
                "error": "Insufficient positive values",
            }

        slope, intercept, r_value, p_value, std_err = scipy_stats.linregress(
            log_doses, log_values
        )

        # 90% CI for slope (beta)
        n = len(log_doses)
        t_crit = scipy_stats.t.ppf(0.95, n - 2)
        beta_lower = slope - t_crit * std_err
        beta_upper = slope + t_crit * std_err

        # Dose proportionality if 90% CI contains 1
        is_proportional = beta_lower <= 1.0 <= beta_upper

        return {
            "method": method,
            "parameter": parameter,
            "n": n,
            "n_doses": len(results_by_dose),
            "doses": sorted(results_by_dose.keys()),
            "beta": slope,
            "beta_se": std_err,
            "beta_lower_90": beta_lower,
            "beta_upper_90": beta_upper,
            "intercept": np.exp(intercept),
            "r_squared": r_value ** 2,
            "p_value": p_value,
            "is_proportional": is_proportional,
            "conclusion": (
                "Dose proportional" if is_proportional
                else "Not dose proportional"
            ),
        }

    elif method == "anova":
        # ANOVA-based approach: compare dose-normalized values
        dn_values = values / doses
        dose_groups = [str(d) for d in doses]
        unique_doses = sorted(set(doses))

        groups = [dn_values[doses == d] for d in unique_doses]

        if any(len(g) < 2 for g in groups):
            return {
                "method": method,
                "parameter": parameter,
                "error": "Need at least 2 subjects per dose group",
            }

        f_stat, p_value = scipy_stats.f_oneway(*groups)

        # If p > 0.05, dose-normalized values are similar = dose proportional
        is_proportional = p_value > 0.05

        return {
            "method": method,
            "parameter": parameter,
            "n": len(values),
            "n_doses": len(unique_doses),
            "doses": unique_doses,
            "f_statistic": f_stat,
            "p_value": p_value,
            "is_proportional": is_proportional,
            "conclusion": (
                "Dose proportional" if is_proportional
                else "Not dose proportional"
            ),
        }

    else:
        raise ValueError(f"Unknown method: {method}")
