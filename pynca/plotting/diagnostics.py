"""Diagnostic plots for NCA analysis."""

from typing import Optional, Dict, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    import matplotlib.axes


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


def plot_lambda_z(
    conc: np.ndarray,
    time: np.ndarray,
    lambda_z_result: Optional[Dict] = None,
    log_y: bool = True,
    title: Optional[str] = None,
    figsize: tuple = (8, 6),
    ax: Optional["matplotlib.axes.Axes"] = None,
    show_all_points: bool = True,
    **kwargs,
) -> "matplotlib.axes.Axes":
    """
    Plot lambda_z regression diagnostic.

    Shows concentration-time data with terminal phase regression line.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    lambda_z_result : dict, optional
        Result from calc_lambda_z(). If None, calculates automatically.
    log_y : bool
        Use logarithmic y-axis (default True for terminal phase)
    title : str, optional
        Plot title
    figsize : tuple
        Figure size (width, height)
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    show_all_points : bool
        Show all data points, not just terminal phase
    **kwargs
        Additional arguments passed to plt.plot()

    Returns
    -------
    matplotlib.axes.Axes
        The axes object
    """
    plt = _get_matplotlib()
    from pynca.calc.half_life import calc_lambda_z

    conc = np.asarray(conc)
    time = np.asarray(time)

    # Calculate lambda_z if not provided
    if lambda_z_result is None:
        lambda_z_result = calc_lambda_z(conc, time)

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot all data points
    if show_all_points:
        ax.scatter(time, conc, color="C0", alpha=0.5, label="All data", zorder=2)

    # Plot terminal phase points
    if lambda_z_result and "points_used" in lambda_z_result:
        points_used = lambda_z_result["points_used"]
        if len(points_used) > 0:
            terminal_time = time[points_used]
            terminal_conc = conc[points_used]
            ax.scatter(
                terminal_time, terminal_conc,
                color="C1", s=100, marker="o",
                label="Terminal phase", zorder=3
            )

            # Plot regression line
            if not np.isnan(lambda_z_result.get("lambda_z", np.nan)):
                lambda_z = lambda_z_result["lambda_z"]
                intercept = lambda_z_result["intercept"]

                t_line = np.linspace(terminal_time.min(), terminal_time.max(), 100)
                c_line = np.exp(intercept - lambda_z * t_line)
                ax.plot(
                    t_line, c_line,
                    color="C1", linestyle="--", linewidth=2,
                    label=f"Regression (λz={lambda_z:.4f})",
                    zorder=1
                )

    # Configure axes
    if log_y:
        ax.set_yscale("log")

    ax.set_xlabel("Time")
    ax.set_ylabel("Concentration")

    if title:
        ax.set_title(title)
    else:
        # Add regression info to title
        if lambda_z_result and not np.isnan(lambda_z_result.get("r_squared", np.nan)):
            r2 = lambda_z_result["r_squared"]
            n = lambda_z_result["n_points"]
            hl = lambda_z_result.get("half_life", np.nan)
            ax.set_title(f"Lambda_z Regression\nR²={r2:.4f}, n={n}, t½={hl:.2f}")

    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax


def plot_residuals(
    conc: np.ndarray,
    time: np.ndarray,
    lambda_z_result: Optional[Dict] = None,
    title: Optional[str] = None,
    figsize: tuple = (10, 4),
    **kwargs,
) -> "plt.Figure":
    """
    Plot residuals from lambda_z regression.

    Creates a two-panel figure with:
    1. Residuals vs time
    2. Q-Q plot of residuals

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    lambda_z_result : dict, optional
        Result from calc_lambda_z(). If None, calculates automatically.
    title : str, optional
        Overall plot title
    figsize : tuple
        Figure size (width, height)
    **kwargs
        Additional arguments

    Returns
    -------
    matplotlib.figure.Figure
        The figure object
    """
    plt = _get_matplotlib()
    from scipy import stats
    from pynca.calc.half_life import calc_lambda_z

    conc = np.asarray(conc)
    time = np.asarray(time)

    # Calculate lambda_z if not provided
    if lambda_z_result is None:
        lambda_z_result = calc_lambda_z(conc, time)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    if lambda_z_result and "points_used" in lambda_z_result and len(lambda_z_result["points_used"]) > 0:
        points_used = lambda_z_result["points_used"]
        terminal_time = time[points_used]
        terminal_conc = conc[points_used]

        # Calculate predicted values
        lambda_z = lambda_z_result["lambda_z"]
        intercept = lambda_z_result["intercept"]
        predicted = np.exp(intercept - lambda_z * terminal_time)

        # Calculate residuals (on log scale)
        residuals = np.log(terminal_conc) - np.log(predicted)

        # Plot 1: Residuals vs time
        ax1.scatter(terminal_time, residuals, color="C0", s=50)
        ax1.axhline(0, color="gray", linestyle="--")
        ax1.set_xlabel("Time")
        ax1.set_ylabel("Residuals (log scale)")
        ax1.set_title("Residuals vs Time")
        ax1.grid(True, alpha=0.3)

        # Plot 2: Q-Q plot
        stats.probplot(residuals, dist="norm", plot=ax2)
        ax2.set_title("Q-Q Plot")
        ax2.grid(True, alpha=0.3)

    else:
        ax1.text(0.5, 0.5, "No regression data", ha="center", va="center",
                 transform=ax1.transAxes)
        ax2.text(0.5, 0.5, "No regression data", ha="center", va="center",
                 transform=ax2.transAxes)

    if title:
        fig.suptitle(title)

    fig.tight_layout()
    return fig
