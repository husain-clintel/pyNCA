"""Concentration interpolation and extrapolation functions."""

from pynca.interpolation.interpolate import interpolate_conc, interpolate_conc_linear, interpolate_conc_log
from pynca.interpolation.extrapolate import extrapolate_conc

__all__ = [
    "interpolate_conc",
    "interpolate_conc_linear",
    "interpolate_conc_log",
    "extrapolate_conc",
]
