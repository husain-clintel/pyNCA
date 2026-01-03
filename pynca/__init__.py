"""
pyNCA: Python Non-Compartmental Analysis Package

A Python package for Pharmacokinetic Non-Compartmental Analysis (NCA),
inspired by the PKNCA R package.
"""

__version__ = "0.1.0"

from pynca.core.concentration import NCAConcentration
from pynca.core.dose import NCADose
from pynca.core.data import NCAData
from pynca.core.results import NCAResults
from pynca.core.options import options, NCAOptions, get_default_options

from pynca.analysis.nca import run_nca, NCA
from pynca.analysis.multi_analyte import (
    run_multi_analyte_nca,
    calc_metabolite_ratio,
    calc_metabolite_ratios,
    summarize_metabolite_ratios,
    molar_ratio,
    MultiAnalyteResults,
)

from pynca.calc.auc import (
    calc_auc,
    calc_auc_last,
    calc_auc_inf,
    calc_auc_all,
    calc_auc_int,
    calc_auc_pct_extrap,
    calc_auc_partial,
    calc_auc_dn,
    calc_cmax_dn,
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
from pynca.calc.half_life import (
    calc_half_life,
    calc_lambda_z,
    calc_effective_half_life,
    calc_time_above_threshold,
    calc_pct_time_above_threshold,
)
from pynca.calc.clearance import calc_cl, calc_vz, calc_vss, calc_mrt
from pynca.calc.bioavailability import calc_f

from pynca.cleaning.blq import clean_blq
from pynca.cleaning.imputation import impute_conc

from pynca.interpolation.interpolate import interpolate_conc
from pynca.interpolation.extrapolate import extrapolate_conc

# Plotting (optional - requires matplotlib)
try:
    from pynca.plotting import (
        plot_conc_time,
        plot_conc_time_by_subject,
        plot_lambda_z,
        plot_residuals,
        plot_parameter_summary,
        plot_forest,
    )
    _HAS_PLOTTING = True
except ImportError:
    _HAS_PLOTTING = False

__all__ = [
    # Version
    "__version__",
    # Core classes
    "NCAConcentration",
    "NCADose",
    "NCAData",
    "NCAResults",
    "NCAOptions",
    "options",
    "get_default_options",
    # Analysis
    "run_nca",
    "NCA",
    # Multi-analyte
    "run_multi_analyte_nca",
    "calc_metabolite_ratio",
    "calc_metabolite_ratios",
    "summarize_metabolite_ratios",
    "molar_ratio",
    "MultiAnalyteResults",
    # AUC calculations
    "calc_auc",
    "calc_auc_last",
    "calc_auc_inf",
    "calc_auc_all",
    "calc_auc_int",
    "calc_auc_pct_extrap",
    "calc_auc_partial",
    "calc_auc_dn",
    "calc_cmax_dn",
    # AUMC calculations
    "calc_aumc",
    "calc_aumc_last",
    "calc_aumc_inf",
    # Peak parameters
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
    "calc_effective_half_life",
    "calc_time_above_threshold",
    "calc_pct_time_above_threshold",
    # Clearance and volume
    "calc_cl",
    "calc_vz",
    "calc_vss",
    "calc_mrt",
    # Bioavailability
    "calc_f",
    # Cleaning
    "clean_blq",
    "impute_conc",
    # Interpolation
    "interpolate_conc",
    "extrapolate_conc",
    # Plotting (optional)
    "plot_conc_time",
    "plot_conc_time_by_subject",
    "plot_lambda_z",
    "plot_residuals",
    "plot_parameter_summary",
    "plot_forest",
]
