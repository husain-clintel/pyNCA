"""NCA parameter registry and metadata."""

from typing import Any, Callable, Dict, List, Optional


# Parameter registry - defines all available NCA parameters
PARAMETERS: Dict[str, Dict[str, Any]] = {
    # Peak parameters
    "cmax": {
        "description": "Maximum observed concentration",
        "unit_type": "concentration",
        "requires": ["conc", "time"],
        "module": "cmax",
        "function": "calc_cmax",
    },
    "tmax": {
        "description": "Time of maximum observed concentration",
        "unit_type": "time",
        "requires": ["conc", "time"],
        "module": "cmax",
        "function": "calc_tmax",
    },
    "cmin": {
        "description": "Minimum observed concentration",
        "unit_type": "concentration",
        "requires": ["conc", "time"],
        "module": "cmax",
        "function": "calc_cmin",
    },
    "tlast": {
        "description": "Time of last measurable concentration",
        "unit_type": "time",
        "requires": ["conc", "time"],
        "module": "cmax",
        "function": "calc_tlast",
    },
    "clast": {
        "description": "Last measurable concentration",
        "unit_type": "concentration",
        "requires": ["conc", "time"],
        "module": "cmax",
        "function": "calc_clast",
    },
    "cav": {
        "description": "Average concentration over interval",
        "unit_type": "concentration",
        "requires": ["auc", "tau"],
        "module": "cmax",
        "function": "calc_cav",
    },
    "ctrough": {
        "description": "Trough concentration",
        "unit_type": "concentration",
        "requires": ["conc", "time"],
        "module": "cmax",
        "function": "calc_ctrough",
    },
    "swing": {
        "description": "Swing (Cmax - Cmin) / Cmin",
        "unit_type": "dimensionless",
        "requires": ["cmax", "cmin"],
        "module": "cmax",
        "function": "calc_swing",
    },
    "ptf": {
        "description": "Peak-trough fluctuation (%)",
        "unit_type": "percent",
        "requires": ["cmax", "cmin", "cav"],
        "module": "cmax",
        "function": "calc_ptf",
    },
    # AUC parameters
    "auc.last": {
        "description": "AUC from time 0 to last measurable concentration",
        "unit_type": "concentration*time",
        "requires": ["conc", "time"],
        "module": "auc",
        "function": "calc_auc_last",
    },
    "auc.all": {
        "description": "AUC from time 0 to last time point",
        "unit_type": "concentration*time",
        "requires": ["conc", "time"],
        "module": "auc",
        "function": "calc_auc_all",
    },
    "auc.inf.obs": {
        "description": "AUC extrapolated to infinity (observed Clast)",
        "unit_type": "concentration*time",
        "requires": ["conc", "time", "lambda.z"],
        "module": "auc",
        "function": "calc_auc_inf",
    },
    "auc.inf.pred": {
        "description": "AUC extrapolated to infinity (predicted Clast)",
        "unit_type": "concentration*time",
        "requires": ["conc", "time", "lambda.z", "clast.pred"],
        "module": "auc",
        "function": "calc_auc_inf",
    },
    "auc.pct.extrap.obs": {
        "description": "Percent of AUC extrapolated (observed)",
        "unit_type": "percent",
        "requires": ["auc.last", "auc.inf.obs"],
        "module": "auc",
        "function": "calc_auc_pct_extrap",
    },
    "auc.pct.extrap.pred": {
        "description": "Percent of AUC extrapolated (predicted)",
        "unit_type": "percent",
        "requires": ["auc.last", "auc.inf.pred"],
        "module": "auc",
        "function": "calc_auc_pct_extrap",
    },
    # AUMC parameters
    "aumc.last": {
        "description": "AUMC from time 0 to last measurable concentration",
        "unit_type": "concentration*time^2",
        "requires": ["conc", "time"],
        "module": "aumc",
        "function": "calc_aumc_last",
    },
    "aumc.inf.obs": {
        "description": "AUMC extrapolated to infinity (observed)",
        "unit_type": "concentration*time^2",
        "requires": ["conc", "time", "lambda.z"],
        "module": "aumc",
        "function": "calc_aumc_inf",
    },
    "aumc.inf.pred": {
        "description": "AUMC extrapolated to infinity (predicted)",
        "unit_type": "concentration*time^2",
        "requires": ["conc", "time", "lambda.z", "clast.pred"],
        "module": "aumc",
        "function": "calc_aumc_inf",
    },
    # Half-life parameters
    "lambda.z": {
        "description": "Terminal elimination rate constant",
        "unit_type": "1/time",
        "requires": ["conc", "time"],
        "module": "half_life",
        "function": "calc_lambda_z",
    },
    "half.life": {
        "description": "Terminal elimination half-life",
        "unit_type": "time",
        "requires": ["lambda.z"],
        "module": "half_life",
        "function": "calc_half_life",
    },
    "r.squared": {
        "description": "R-squared of lambda_z regression",
        "unit_type": "dimensionless",
        "requires": ["lambda.z"],
        "module": "half_life",
        "function": "calc_lambda_z",
    },
    "adj.r.squared": {
        "description": "Adjusted R-squared of lambda_z regression",
        "unit_type": "dimensionless",
        "requires": ["lambda.z"],
        "module": "half_life",
        "function": "calc_lambda_z",
    },
    "lambda.z.n.points": {
        "description": "Number of points used for lambda_z",
        "unit_type": "count",
        "requires": ["lambda.z"],
        "module": "half_life",
        "function": "calc_lambda_z",
    },
    "clast.pred": {
        "description": "Predicted Clast from lambda_z regression",
        "unit_type": "concentration",
        "requires": ["lambda.z", "tlast"],
        "module": "half_life",
        "function": "calc_clast_pred",
    },
    # Clearance parameters
    "cl.obs": {
        "description": "Clearance (observed, CL or CL/F)",
        "unit_type": "volume/time",
        "requires": ["dose", "auc.inf.obs"],
        "module": "clearance",
        "function": "calc_cl",
    },
    "cl.pred": {
        "description": "Clearance (predicted, CL or CL/F)",
        "unit_type": "volume/time",
        "requires": ["dose", "auc.inf.pred"],
        "module": "clearance",
        "function": "calc_cl",
    },
    "cl.last": {
        "description": "Clearance based on AUClast",
        "unit_type": "volume/time",
        "requires": ["dose", "auc.last"],
        "module": "clearance",
        "function": "calc_cl_last",
    },
    # Volume parameters
    "vz.obs": {
        "description": "Volume of distribution (observed, Vz or Vz/F)",
        "unit_type": "volume",
        "requires": ["cl.obs", "lambda.z"],
        "module": "clearance",
        "function": "calc_vz",
    },
    "vz.pred": {
        "description": "Volume of distribution (predicted)",
        "unit_type": "volume",
        "requires": ["cl.pred", "lambda.z"],
        "module": "clearance",
        "function": "calc_vz",
    },
    "vss.obs": {
        "description": "Volume of distribution at steady state (observed)",
        "unit_type": "volume",
        "requires": ["dose", "aumc.inf.obs", "auc.inf.obs"],
        "module": "clearance",
        "function": "calc_vss",
    },
    "vss.pred": {
        "description": "Volume of distribution at steady state (predicted)",
        "unit_type": "volume",
        "requires": ["dose", "aumc.inf.pred", "auc.inf.pred"],
        "module": "clearance",
        "function": "calc_vss",
    },
    # MRT parameters
    "mrt.last": {
        "description": "Mean residence time to Tlast",
        "unit_type": "time",
        "requires": ["aumc.last", "auc.last"],
        "module": "clearance",
        "function": "calc_mrt_last",
    },
    "mrt.inf.obs": {
        "description": "Mean residence time (observed)",
        "unit_type": "time",
        "requires": ["aumc.inf.obs", "auc.inf.obs"],
        "module": "clearance",
        "function": "calc_mrt",
    },
    "mrt.inf.pred": {
        "description": "Mean residence time (predicted)",
        "unit_type": "time",
        "requires": ["aumc.inf.pred", "auc.inf.pred"],
        "module": "clearance",
        "function": "calc_mrt",
    },
    # Bioavailability parameters
    "f": {
        "description": "Bioavailability",
        "unit_type": "dimensionless",
        "requires": ["auc.test", "auc.ref", "dose.test", "dose.ref"],
        "module": "bioavailability",
        "function": "calc_f",
    },
    "accumulation.index": {
        "description": "Accumulation index (Rac)",
        "unit_type": "dimensionless",
        "requires": ["half.life", "tau"],
        "module": "bioavailability",
        "function": "calc_accumulation_index",
    },
}


def get_parameter_info(param_name: str) -> Optional[Dict[str, Any]]:
    """
    Get information about a specific parameter.

    Parameters
    ----------
    param_name : str
        Parameter name

    Returns
    -------
    dict or None
        Parameter metadata if found, None otherwise
    """
    return PARAMETERS.get(param_name)


def list_parameters(
    category: Optional[str] = None,
) -> List[str]:
    """
    List available parameters.

    Parameters
    ----------
    category : str, optional
        Filter by category (e.g., "auc", "clearance")

    Returns
    -------
    list of str
        Parameter names
    """
    if category is None:
        return list(PARAMETERS.keys())

    return [
        name
        for name, info in PARAMETERS.items()
        if info.get("module") == category
    ]


def get_parameter_description(param_name: str) -> str:
    """
    Get description for a parameter.

    Parameters
    ----------
    param_name : str
        Parameter name

    Returns
    -------
    str
        Description
    """
    info = PARAMETERS.get(param_name, {})
    return info.get("description", "Unknown parameter")


def get_parameter_unit_type(param_name: str) -> str:
    """
    Get unit type for a parameter.

    Parameters
    ----------
    param_name : str
        Parameter name

    Returns
    -------
    str
        Unit type
    """
    info = PARAMETERS.get(param_name, {})
    return info.get("unit_type", "unknown")


def get_parameter_requirements(param_name: str) -> List[str]:
    """
    Get requirements for calculating a parameter.

    Parameters
    ----------
    param_name : str
        Parameter name

    Returns
    -------
    list of str
        Required inputs
    """
    info = PARAMETERS.get(param_name, {})
    return info.get("requires", [])
