"""Unit handling and conversion."""

from typing import Dict, Optional, Tuple
import numpy as np

# Try to import pint, make optional
try:
    import pint
    ureg = pint.UnitRegistry()
    PINT_AVAILABLE = True
except ImportError:
    PINT_AVAILABLE = False
    ureg = None


# Unit type mappings for NCA parameters
PARAMETER_UNIT_TYPES = {
    "cmax": "concentration",
    "cmin": "concentration",
    "clast": "concentration",
    "clast.pred": "concentration",
    "cav": "concentration",
    "ctrough": "concentration",
    "tmax": "time",
    "tlast": "time",
    "tfirst": "time",
    "auc.last": "concentration*time",
    "auc.all": "concentration*time",
    "auc.inf.obs": "concentration*time",
    "auc.inf.pred": "concentration*time",
    "aumc.last": "concentration*time^2",
    "aumc.inf.obs": "concentration*time^2",
    "aumc.inf.pred": "concentration*time^2",
    "lambda.z": "1/time",
    "half.life": "time",
    "cl.obs": "volume/time",
    "cl.pred": "volume/time",
    "cl.last": "volume/time",
    "vz.obs": "volume",
    "vz.pred": "volume",
    "vss.obs": "volume",
    "vss.pred": "volume",
    "mrt.last": "time",
    "mrt.inf.obs": "time",
    "mrt.inf.pred": "time",
    "auc.pct.extrap.obs": "percent",
    "auc.pct.extrap.pred": "percent",
    "r.squared": "dimensionless",
    "adj.r.squared": "dimensionless",
    "lambda.z.n.points": "count",
    "f": "dimensionless",
    "swing": "dimensionless",
    "ptf": "percent",
}


def derive_parameter_units(
    param_name: str,
    conc_unit: str,
    time_unit: str,
    dose_unit: Optional[str] = None,
) -> str:
    """
    Derive appropriate units for an NCA parameter.

    Parameters
    ----------
    param_name : str
        Parameter name
    conc_unit : str
        Concentration unit (e.g., "ng/mL")
    time_unit : str
        Time unit (e.g., "h")
    dose_unit : str, optional
        Dose unit (e.g., "mg")

    Returns
    -------
    str
        Derived unit string
    """
    unit_type = PARAMETER_UNIT_TYPES.get(param_name, "unknown")

    if unit_type == "concentration":
        return conc_unit
    elif unit_type == "time":
        return time_unit
    elif unit_type == "concentration*time":
        return f"{conc_unit}*{time_unit}"
    elif unit_type == "concentration*time^2":
        return f"{conc_unit}*{time_unit}^2"
    elif unit_type == "1/time":
        return f"1/{time_unit}"
    elif unit_type == "volume":
        # Derive from dose/concentration
        if dose_unit and conc_unit:
            return _derive_volume_unit(dose_unit, conc_unit)
        return "L"  # Default
    elif unit_type == "volume/time":
        vol_unit = _derive_volume_unit(dose_unit, conc_unit) if dose_unit else "L"
        return f"{vol_unit}/{time_unit}"
    elif unit_type == "percent":
        return "%"
    elif unit_type == "dimensionless":
        return ""
    elif unit_type == "count":
        return ""
    else:
        return ""


def _derive_volume_unit(dose_unit: str, conc_unit: str) -> str:
    """Derive volume unit from dose and concentration units."""
    # Simple heuristics for common cases
    if "ng" in conc_unit.lower() and "mg" in dose_unit.lower():
        return "L"
    elif "ug" in conc_unit.lower() or "µg" in conc_unit.lower():
        if "mg" in dose_unit.lower():
            return "L"
    return "L"  # Default


def assign_units(
    results: "NCAResults",
    conc_unit: str,
    time_unit: str,
    dose_unit: Optional[str] = None,
) -> "NCAResults":
    """
    Assign units to NCA results.

    Parameters
    ----------
    results : NCAResults
        NCA results to annotate
    conc_unit : str
        Concentration unit
    time_unit : str
        Time unit
    dose_unit : str, optional
        Dose unit

    Returns
    -------
    NCAResults
        Results with units assigned
    """
    # Add unit column to results
    for result in results._results:
        param = result.get("parameter", "")
        unit = derive_parameter_units(param, conc_unit, time_unit, dose_unit)
        result["unit"] = unit

    results._data = None  # Invalidate cache
    return results


def convert_units(
    value: float,
    from_unit: str,
    to_unit: str,
) -> float:
    """
    Convert a value between units.

    Parameters
    ----------
    value : float
        Value to convert
    from_unit : str
        Source unit
    to_unit : str
        Target unit

    Returns
    -------
    float
        Converted value

    Raises
    ------
    ImportError
        If pint is not installed
    ValueError
        If conversion is not possible
    """
    if not PINT_AVAILABLE:
        raise ImportError(
            "Unit conversion requires the 'pint' package. "
            "Install with: pip install pint"
        )

    if np.isnan(value):
        return np.nan

    try:
        quantity = value * ureg(from_unit)
        return quantity.to(to_unit).magnitude
    except Exception as e:
        raise ValueError(f"Cannot convert from '{from_unit}' to '{to_unit}': {e}")


def get_unit_type(param_name: str) -> str:
    """
    Get the unit type for a parameter.

    Parameters
    ----------
    param_name : str
        Parameter name

    Returns
    -------
    str
        Unit type
    """
    return PARAMETER_UNIT_TYPES.get(param_name, "unknown")


# Common unit conversions
CONCENTRATION_UNITS = {
    "ng/mL": 1.0,
    "ug/mL": 1000.0,
    "µg/mL": 1000.0,
    "mg/mL": 1e6,
    "pg/mL": 0.001,
    "ng/L": 0.001,
    "ug/L": 1.0,
    "mg/L": 1000.0,
}

TIME_UNITS = {
    "h": 1.0,
    "hr": 1.0,
    "hour": 1.0,
    "hours": 1.0,
    "min": 1 / 60,
    "minute": 1 / 60,
    "minutes": 1 / 60,
    "d": 24.0,
    "day": 24.0,
    "days": 24.0,
    "wk": 168.0,
    "week": 168.0,
    "weeks": 168.0,
}


def simple_convert(
    value: float,
    from_unit: str,
    to_unit: str,
    unit_type: str = "concentration",
) -> float:
    """
    Simple unit conversion without pint dependency.

    Parameters
    ----------
    value : float
        Value to convert
    from_unit : str
        Source unit
    to_unit : str
        Target unit
    unit_type : str
        Type of unit ("concentration" or "time")

    Returns
    -------
    float
        Converted value
    """
    if np.isnan(value):
        return np.nan

    if from_unit == to_unit:
        return value

    if unit_type == "concentration":
        units = CONCENTRATION_UNITS
    elif unit_type == "time":
        units = TIME_UNITS
    else:
        raise ValueError(f"Unknown unit type: {unit_type}")

    if from_unit not in units:
        raise ValueError(f"Unknown {unit_type} unit: {from_unit}")
    if to_unit not in units:
        raise ValueError(f"Unknown {unit_type} unit: {to_unit}")

    # Convert to base unit, then to target
    base_value = value * units[from_unit]
    return base_value / units[to_unit]
