"""Concentration-time plots."""

from typing import Optional, List, Union, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    import matplotlib.axes
    from pynca.core.concentration import NCAConcentration
    from pynca.core.data import NCAData


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


def plot_conc_time(
    data: Union["NCAConcentration", "NCAData"],
    subjects: Optional[List] = None,
    log_y: bool = False,
    show_dose: bool = True,
    title: Optional[str] = None,
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: tuple = (10, 6),
    ax: Optional["matplotlib.axes.Axes"] = None,
    **kwargs,
) -> "matplotlib.axes.Axes":
    """
    Plot concentration-time data.

    Parameters
    ----------
    data : NCAConcentration or NCAData
        Concentration data to plot
    subjects : list, optional
        List of subjects to include. None = all subjects
    log_y : bool
        Use logarithmic y-axis
    show_dose : bool
        Show dose markers (if NCAData with dose info)
    title : str, optional
        Plot title
    xlabel : str
        X-axis label
    ylabel : str
        Y-axis label
    figsize : tuple
        Figure size (width, height)
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    **kwargs
        Additional arguments passed to plt.plot()

    Returns
    -------
    matplotlib.axes.Axes
        The axes object
    """
    plt = _get_matplotlib()

    # Import here to avoid circular imports
    from pynca.core.concentration import NCAConcentration
    from pynca.core.data import NCAData

    # Get concentration object
    if isinstance(data, NCAData):
        conc_obj = data.conc
        dose_obj = data.dose
    elif isinstance(data, NCAConcentration):
        conc_obj = data
        dose_obj = None
    else:
        raise TypeError("data must be NCAConcentration or NCAData")

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Get subjects
    if subjects is None:
        subjects = conc_obj.subjects

    # Default line style
    default_kwargs = {"marker": "o", "linestyle": "-", "markersize": 5}
    default_kwargs.update(kwargs)

    # Plot each subject
    for subject in subjects:
        conc_array, time_array = conc_obj.get_conc_time(subject)
        label = f"Subject {subject}" if len(subjects) > 1 else None
        ax.plot(time_array, conc_array, label=label, **default_kwargs)

        # Show dose markers if available
        if show_dose and dose_obj is not None:
            dose_time, dose_amount = dose_obj.get_first_dose(subject)
            if dose_time is not None:
                ax.axvline(dose_time, color="gray", linestyle="--", alpha=0.5)

    # Configure axes
    if log_y:
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title)
    elif len(subjects) == 1:
        ax.set_title(f"Subject {subjects[0]}")

    if len(subjects) > 1:
        ax.legend()

    ax.grid(True, alpha=0.3)

    return ax


def plot_conc_time_by_subject(
    data: Union["NCAConcentration", "NCAData"],
    subjects: Optional[List] = None,
    ncols: int = 3,
    log_y: bool = False,
    figsize: Optional[tuple] = None,
    sharex: bool = True,
    sharey: bool = True,
    **kwargs,
) -> "plt.Figure":
    """
    Plot concentration-time data with one subplot per subject.

    Parameters
    ----------
    data : NCAConcentration or NCAData
        Concentration data to plot
    subjects : list, optional
        List of subjects to include. None = all subjects
    ncols : int
        Number of columns in subplot grid
    log_y : bool
        Use logarithmic y-axis
    figsize : tuple, optional
        Figure size. If None, calculated automatically.
    sharex : bool
        Share x-axis across subplots
    sharey : bool
        Share y-axis across subplots
    **kwargs
        Additional arguments passed to plt.plot()

    Returns
    -------
    matplotlib.figure.Figure
        The figure object
    """
    plt = _get_matplotlib()

    # Import here to avoid circular imports
    from pynca.core.concentration import NCAConcentration
    from pynca.core.data import NCAData

    # Get concentration object
    if isinstance(data, NCAData):
        conc_obj = data.conc
    elif isinstance(data, NCAConcentration):
        conc_obj = data
    else:
        raise TypeError("data must be NCAConcentration or NCAData")

    # Get subjects
    if subjects is None:
        subjects = conc_obj.subjects

    n_subjects = len(subjects)
    nrows = int(np.ceil(n_subjects / ncols))

    # Calculate figure size
    if figsize is None:
        figsize = (4 * ncols, 3 * nrows)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=figsize,
        sharex=sharex,
        sharey=sharey,
        squeeze=False
    )

    # Default line style
    default_kwargs = {"marker": "o", "linestyle": "-", "markersize": 4, "color": "C0"}
    default_kwargs.update(kwargs)

    # Plot each subject
    for idx, subject in enumerate(subjects):
        row = idx // ncols
        col = idx % ncols
        ax = axes[row, col]

        conc_array, time_array = conc_obj.get_conc_time(subject)
        ax.plot(time_array, conc_array, **default_kwargs)

        if log_y:
            ax.set_yscale("log")

        ax.set_title(f"Subject {subject}")
        ax.grid(True, alpha=0.3)

    # Hide empty subplots
    for idx in range(n_subjects, nrows * ncols):
        row = idx // ncols
        col = idx % ncols
        axes[row, col].set_visible(False)

    # Add common labels
    fig.supxlabel("Time")
    fig.supylabel("Concentration")
    fig.tight_layout()

    return fig
