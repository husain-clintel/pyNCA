"""NCA calculation functions."""

from pynca.calc.auc import (
    calc_auc,
    calc_auc_last,
    calc_auc_inf,
    calc_auc_all,
    calc_auc_int,
    calc_auc_pct_extrap,
)
from pynca.calc.aumc import calc_aumc, calc_aumc_last, calc_aumc_inf
from pynca.calc.cmax import (
    calc_cmax,
    calc_cmin,
    calc_tmax,
    calc_tlast,
    calc_clast,
    calc_cav,
    calc_ctrough,
    calc_swing,
    calc_ptf,
)
from pynca.calc.half_life import calc_half_life, calc_lambda_z
from pynca.calc.clearance import calc_cl, calc_vz, calc_vss, calc_mrt
from pynca.calc.bioavailability import calc_f
from pynca.calc.parameters import PARAMETERS, get_parameter_info, list_parameters

__all__ = [
    # AUC
    "calc_auc",
    "calc_auc_last",
    "calc_auc_inf",
    "calc_auc_all",
    "calc_auc_int",
    "calc_auc_pct_extrap",
    # AUMC
    "calc_aumc",
    "calc_aumc_last",
    "calc_aumc_inf",
    # Peak
    "calc_cmax",
    "calc_cmin",
    "calc_tmax",
    "calc_tlast",
    "calc_clast",
    "calc_cav",
    "calc_ctrough",
    "calc_swing",
    "calc_ptf",
    # Half-life
    "calc_half_life",
    "calc_lambda_z",
    # Clearance
    "calc_cl",
    "calc_vz",
    "calc_vss",
    "calc_mrt",
    # Bioavailability
    "calc_f",
    # Registry
    "PARAMETERS",
    "get_parameter_info",
    "list_parameters",
]
