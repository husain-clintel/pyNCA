"""NCA analysis functions."""

# Import intervals first (no circular dependencies)
from pynca.analysis.intervals import Interval, create_intervals, auto_select_intervals

# Then import modules that may depend on core
from pynca.analysis.sparse import sparse_auc, sparse_auc_se
from pynca.analysis.steady_state import time_to_steady_state, check_steady_state
from pynca.analysis.superposition import superposition

# Import nca last (depends on core.data)
from pynca.analysis.nca import NCA, run_nca

# Multi-analyte support
from pynca.analysis.multi_analyte import (
    run_multi_analyte_nca,
    calc_metabolite_ratio,
    calc_metabolite_ratios,
    summarize_metabolite_ratios,
    molar_ratio,
    MultiAnalyteResults,
)

__all__ = [
    "NCA",
    "run_nca",
    "Interval",
    "create_intervals",
    "auto_select_intervals",
    "sparse_auc",
    "sparse_auc_se",
    "time_to_steady_state",
    "check_steady_state",
    "superposition",
    # Multi-analyte
    "run_multi_analyte_nca",
    "calc_metabolite_ratio",
    "calc_metabolite_ratios",
    "summarize_metabolite_ratios",
    "molar_ratio",
    "MultiAnalyteResults",
]
