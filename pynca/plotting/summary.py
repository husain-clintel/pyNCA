"""Summary plots for NCA results."""

from typing import Optional, List, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    import matplotlib.axes
    from pynca.core.results import NCAResults


def _get_matplotlib():
    """Import matplotlib lazily."""
    try:
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        raise ImportError(
            "matplotlib is required for plotting. "
            "Install with: pip install matplotlib"
        )


def plot_parameter_summary(
    results: "NCAResults",
    parameters: Optional[List[str]] = None,
    plot_type: str = "box",
    figsize: Optional[tuple] = None,
    title: Optional[str] = None,
    **kwargs,
) -> "plt.Figure":
    """
    Plot summary of NCA parameters across subjects.

    Parameters
    ----------
    results : NCAResults
        NCA results object
    parameters : list of str, optional
        Parameters to plot. None = select common parameters.
    plot_type : str
        Type of plot: "box", "violin", "strip", or "bar"
    figsize : tuple, optional
        Figure size. If None, calculated automatically.
    title : str, optional
        Plot title
    **kwargs
        Additional arguments passed to plot function

    Returns
    -------
    matplotlib.figure.Figure
        The figure object
    """
    plt = _get_matplotlib()

    # Get results dataframe
    df = results.to_dataframe()

    # Select parameters
    if parameters is None:
        # Default: select common parameters
        common_params = ["cmax", "tmax", "auc.last", "auc.inf.obs", "half.life", "cl.obs"]
        available = df["parameter"].unique()
        parameters = [p for p in common_params if p in available]

    n_params = len(parameters)

    # Calculate figure size
    if figsize is None:
        figsize = (max(8, n_params * 1.5), 6)

    fig, ax = plt.subplots(figsize=figsize)

    # Prepare data for plotting
    data_to_plot = []
    labels = []
    for param in parameters:
        param_data = df[df["parameter"] == param]["value"].dropna()
        if len(param_data) > 0:
            data_to_plot.append(param_data.values)
            labels.append(param)

    if len(data_to_plot) == 0:
        ax.text(0.5, 0.5, "No data to plot", ha="center", va="center",
                transform=ax.transAxes)
        return fig

    # Create plot based on type
    positions = range(1, len(data_to_plot) + 1)

    if plot_type == "box":
        bp = ax.boxplot(data_to_plot, positions=positions, **kwargs)
        ax.set_xticklabels(labels)

    elif plot_type == "violin":
        vp = ax.violinplot(data_to_plot, positions=positions, **kwargs)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels)

    elif plot_type == "strip":
        for i, (data, label) in enumerate(zip(data_to_plot, labels), 1):
            # Add jitter
            x = np.random.normal(i, 0.04, size=len(data))
            ax.scatter(x, data, alpha=0.7, s=50, **kwargs)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels)

    elif plot_type == "bar":
        means = [np.mean(d) for d in data_to_plot]
        stds = [np.std(d) for d in data_to_plot]
        ax.bar(positions, means, yerr=stds, capsize=5, **kwargs)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels)

    else:
        raise ValueError(f"Unknown plot_type: {plot_type}")

    ax.set_ylabel("Value")
    ax.set_xlabel("Parameter")

    if title:
        ax.set_title(title)
    else:
        ax.set_title("NCA Parameter Summary")

    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()

    return fig


def plot_forest(
    results: "NCAResults",
    parameter: str,
    reference_value: Optional[float] = None,
    ci_level: float = 0.95,
    figsize: tuple = (8, 6),
    title: Optional[str] = None,
    **kwargs,
) -> "matplotlib.axes.Axes":
    """
    Create a forest plot for a single parameter.

    Shows individual subject values with overall mean and confidence interval.

    Parameters
    ----------
    results : NCAResults
        NCA results object
    parameter : str
        Parameter to plot (e.g., "cmax", "auc.last")
    reference_value : float, optional
        Reference value to show as vertical line
    ci_level : float
        Confidence interval level (default 0.95 for 95% CI)
    figsize : tuple
        Figure size
    title : str, optional
        Plot title
    **kwargs
        Additional arguments

    Returns
    -------
    matplotlib.axes.Axes
        The axes object
    """
    plt = _get_matplotlib()
    from scipy import stats

    # Get parameter values
    values = results.get_parameter(parameter)
    if len(values) == 0:
        raise ValueError(f"No data for parameter: {parameter}")

    # Get subject identifiers
    df = results.to_dataframe()
    param_df = df[df["parameter"] == parameter]
    subjects = param_df["subject"].values

    n = len(values)

    # Calculate statistics
    mean_val = np.mean(values)
    std_val = np.std(values, ddof=1)
    sem = std_val / np.sqrt(n)

    # Calculate confidence interval
    t_crit = stats.t.ppf((1 + ci_level) / 2, n - 1)
    ci_low = mean_val - t_crit * sem
    ci_high = mean_val + t_crit * sem

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot individual values
    y_positions = range(n)
    ax.scatter(values, y_positions, color="C0", s=80, zorder=3, label="Individual")

    # Plot mean and CI
    ax.axvline(mean_val, color="C1", linewidth=2, label=f"Mean: {mean_val:.3f}")
    ax.axvspan(ci_low, ci_high, color="C1", alpha=0.2, label=f"{ci_level*100:.0f}% CI")

    # Plot reference value if provided
    if reference_value is not None:
        ax.axvline(reference_value, color="gray", linestyle="--",
                   label=f"Reference: {reference_value}")

    # Configure axes
    ax.set_yticks(y_positions)
    ax.set_yticklabels([f"Subject {s}" for s in subjects])
    ax.set_xlabel(parameter)

    if title:
        ax.set_title(title)
    else:
        ax.set_title(f"Forest Plot: {parameter}")

    ax.legend(loc="best")
    ax.grid(True, axis="x", alpha=0.3)

    fig.tight_layout()

    return ax


def plot_pk_profile(
    results: "NCAResults",
    figsize: tuple = (12, 8),
    title: Optional[str] = None,
) -> "plt.Figure":
    """
    Create a multi-panel PK profile summary.

    Parameters
    ----------
    results : NCAResults
        NCA results object
    figsize : tuple
        Figure size
    title : str, optional
        Overall title

    Returns
    -------
    matplotlib.figure.Figure
        The figure object
    """
    plt = _get_matplotlib()

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    df = results.to_dataframe()

    # Panel 1: Cmax distribution
    ax1 = axes[0, 0]
    if "cmax" in df["parameter"].values:
        cmax_vals = df[df["parameter"] == "cmax"]["value"].dropna()
        ax1.hist(cmax_vals, bins="auto", edgecolor="black", alpha=0.7)
        ax1.axvline(cmax_vals.mean(), color="red", linestyle="--", label="Mean")
        ax1.set_xlabel("Cmax")
        ax1.set_ylabel("Frequency")
        ax1.set_title("Cmax Distribution")
        ax1.legend()

    # Panel 2: AUC distribution
    ax2 = axes[0, 1]
    auc_col = "auc.inf.obs" if "auc.inf.obs" in df["parameter"].values else "auc.last"
    if auc_col in df["parameter"].values:
        auc_vals = df[df["parameter"] == auc_col]["value"].dropna()
        ax2.hist(auc_vals, bins="auto", edgecolor="black", alpha=0.7)
        ax2.axvline(auc_vals.mean(), color="red", linestyle="--", label="Mean")
        ax2.set_xlabel(auc_col)
        ax2.set_ylabel("Frequency")
        ax2.set_title(f"{auc_col} Distribution")
        ax2.legend()

    # Panel 3: Half-life distribution
    ax3 = axes[1, 0]
    if "half.life" in df["parameter"].values:
        hl_vals = df[df["parameter"] == "half.life"]["value"].dropna()
        ax3.hist(hl_vals, bins="auto", edgecolor="black", alpha=0.7)
        ax3.axvline(hl_vals.mean(), color="red", linestyle="--", label="Mean")
        ax3.set_xlabel("Half-life")
        ax3.set_ylabel("Frequency")
        ax3.set_title("Half-life Distribution")
        ax3.legend()

    # Panel 4: Clearance distribution
    ax4 = axes[1, 1]
    if "cl.obs" in df["parameter"].values:
        cl_vals = df[df["parameter"] == "cl.obs"]["value"].dropna()
        ax4.hist(cl_vals, bins="auto", edgecolor="black", alpha=0.7)
        ax4.axvline(cl_vals.mean(), color="red", linestyle="--", label="Mean")
        ax4.set_xlabel("CL/F")
        ax4.set_ylabel("Frequency")
        ax4.set_title("Clearance Distribution")
        ax4.legend()

    if title:
        fig.suptitle(title, fontsize=14)

    fig.tight_layout()
    return fig
