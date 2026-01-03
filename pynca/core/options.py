"""Global options and configuration for pyNCA."""

from typing import Any, Dict, List, Optional, Union
from dataclasses import dataclass, field
import copy


@dataclass
class NCAOptions:
    """
    Configuration options for NCA calculations.

    Attributes
    ----------
    auc_method : str
        Method for AUC calculation: "linear", "log", "linear-up/log-down"
    lambda_z_method : str
        Method for lambda_z estimation: "ols", "wls"
    lambda_z_n_min : int
        Minimum number of points for lambda_z calculation
    lambda_z_n_max : int or None
        Maximum number of points for lambda_z calculation (None = no limit)
    lambda_z_selection : str
        Method for automatic point selection: "best_fit", "manual"
    blq_first_handling : str
        Handling of BLQ values before first measurable: "zero", "drop"
    blq_last_handling : str
        Handling of BLQ values after last measurable: "zero", "drop"
    blq_middle_handling : str
        Handling of BLQ values between measurables: "zero", "drop", "interpolate"
    loq : float or None
        Limit of quantification
    conc_units : str or None
        Concentration units
    time_units : str or None
        Time units
    dose_units : str or None
        Dose units
    adj_r_squared_threshold : float
        Minimum adjusted R-squared for lambda_z acceptance
    extrapolation_max_pct : float
        Maximum allowed percent extrapolation for AUCinf
    summary_stats : list
        Default statistics for summary tables
    """

    # AUC calculation options
    auc_method: str = "linear"

    # Lambda Z / Half-life options
    lambda_z_method: str = "ols"
    lambda_z_n_min: int = 3
    lambda_z_n_max: Optional[int] = None
    lambda_z_selection: str = "best_fit"
    adj_r_squared_threshold: float = 0.0

    # BLQ handling options
    blq_first_handling: str = "zero"
    blq_last_handling: str = "drop"
    blq_middle_handling: str = "zero"
    loq: Optional[float] = None

    # Unit options
    conc_units: Optional[str] = None
    time_units: Optional[str] = None
    dose_units: Optional[str] = None

    # Extrapolation options
    extrapolation_max_pct: float = 20.0

    # Summary options
    summary_stats: List[str] = field(
        default_factory=lambda: ["n", "mean", "sd", "cv", "median", "min", "max"]
    )

    # Interval options
    default_parameters: List[str] = field(
        default_factory=lambda: [
            "cmax",
            "tmax",
            "tlast",
            "clast",
            "auc.last",
            "auc.inf.obs",
            "auc.pct.extrap.obs",
            "aumc.last",
            "aumc.inf.obs",
            "lambda.z",
            "half.life",
            "cl.obs",
            "vz.obs",
            "vss.obs",
            "mrt.last",
        ]
    )

    def set(self, **kwargs) -> "NCAOptions":
        """
        Set multiple options at once.

        Parameters
        ----------
        **kwargs
            Option names and values

        Returns
        -------
        NCAOptions
            Self for chaining
        """
        for key, value in kwargs.items():
            if not hasattr(self, key):
                raise ValueError(f"Unknown option: {key}")
            setattr(self, key, value)
        return self

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get an option value.

        Parameters
        ----------
        key : str
            Option name
        default : any
            Default value if option not set

        Returns
        -------
        any
            Option value
        """
        return getattr(self, key, default)

    def copy(self) -> "NCAOptions":
        """Create a copy of options."""
        return copy.deepcopy(self)

    def to_dict(self) -> Dict[str, Any]:
        """Convert options to dictionary."""
        return {
            "auc_method": self.auc_method,
            "lambda_z_method": self.lambda_z_method,
            "lambda_z_n_min": self.lambda_z_n_min,
            "lambda_z_n_max": self.lambda_z_n_max,
            "lambda_z_selection": self.lambda_z_selection,
            "adj_r_squared_threshold": self.adj_r_squared_threshold,
            "blq_first_handling": self.blq_first_handling,
            "blq_last_handling": self.blq_last_handling,
            "blq_middle_handling": self.blq_middle_handling,
            "loq": self.loq,
            "conc_units": self.conc_units,
            "time_units": self.time_units,
            "dose_units": self.dose_units,
            "extrapolation_max_pct": self.extrapolation_max_pct,
            "summary_stats": self.summary_stats,
            "default_parameters": self.default_parameters,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "NCAOptions":
        """Create options from dictionary."""
        return cls(**d)


# Global options instance
_global_options = NCAOptions()


class OptionsProxy:
    """
    Proxy class for accessing and modifying global options.

    This allows using `options.auc_method` syntax for both getting and setting.
    """

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            return object.__getattribute__(self, name)
        return getattr(_global_options, name)

    def __setattr__(self, name: str, value: Any) -> None:
        if name.startswith("_"):
            object.__setattr__(self, name, value)
        else:
            setattr(_global_options, name, value)

    def set(self, **kwargs) -> "OptionsProxy":
        """Set multiple options at once."""
        _global_options.set(**kwargs)
        return self

    def get(self, key: str, default: Any = None) -> Any:
        """Get an option value."""
        return _global_options.get(key, default)

    def reset(self) -> "OptionsProxy":
        """Reset all options to defaults."""
        global _global_options
        _global_options = NCAOptions()
        return self

    def copy(self) -> NCAOptions:
        """Get a copy of current options."""
        return _global_options.copy()

    def to_dict(self) -> Dict[str, Any]:
        """Convert options to dictionary."""
        return _global_options.to_dict()


# Global options proxy
options = OptionsProxy()


def get_default_options() -> NCAOptions:
    """Get a new instance with default options."""
    return NCAOptions()


def get_current_options() -> NCAOptions:
    """Get a copy of current global options."""
    return _global_options.copy()
