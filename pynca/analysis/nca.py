"""Main NCA analysis engine."""

from typing import Any, Dict, List, Optional, Union
import numpy as np

from pynca.core.data import NCAData
from pynca.core.results import NCAResults
from pynca.core.options import NCAOptions, get_default_options, options as global_options
from pynca.analysis.intervals import Interval

from pynca.calc.auc import calc_auc_last, calc_auc_all, calc_auc_inf, calc_auc_pct_extrap
from pynca.calc.aumc import calc_aumc_last, calc_aumc_inf
from pynca.calc.cmax import calc_cmax, calc_tmax, calc_cmin, calc_tlast, calc_clast
from pynca.calc.half_life import calc_lambda_z, calc_clast_pred
from pynca.calc.clearance import calc_cl, calc_vz, calc_vss, calc_mrt, calc_mrt_last
from pynca.cleaning.blq import clean_blq


class NCA:
    """
    Main NCA analysis class.

    This class performs non-compartmental analysis on pharmacokinetic data.

    Parameters
    ----------
    data : NCAData
        Combined concentration-dose data
    options : NCAOptions, optional
        Analysis options (uses global options if not provided)

    Examples
    --------
    >>> nca = NCA(data)
    >>> results = nca.calculate()
    >>> print(results.to_dataframe())
    """

    def __init__(
        self,
        data: NCAData,
        options: Optional[NCAOptions] = None,
    ):
        self.data = data
        self.options = options or global_options.copy()
        self._results = None

    def calculate(
        self,
        parameters: Optional[List[str]] = None,
    ) -> NCAResults:
        """
        Run NCA calculations.

        Parameters
        ----------
        parameters : list of str, optional
            Parameters to calculate. If None, calculates all default parameters.

        Returns
        -------
        NCAResults
            Calculated NCA results
        """
        if parameters is None:
            parameters = self.options.default_parameters

        results = NCAResults()

        # Process each subject-interval combination
        for subject, interval, data_dict in self.data.iter_analysis_units():
            interval_results = self._calculate_interval(
                subject=subject,
                interval=interval,
                conc=data_dict["conc"],
                time=data_dict["time"],
                dose=data_dict["dose"],
                duration=data_dict.get("duration", 0.0),
                parameters=parameters,
            )

            results.add_results(interval_results)

        self._results = results
        return results

    def _calculate_interval(
        self,
        subject: Any,
        interval: Interval,
        conc: np.ndarray,
        time: np.ndarray,
        dose: Optional[float],
        duration: float = 0.0,
        parameters: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """Calculate all parameters for a single interval."""
        results = []

        # Clean BLQ values
        conc_clean, time_clean = clean_blq(
            conc,
            time,
            method="standard",
            loq=self.options.loq,
            first_handling=self.options.blq_first_handling,
            middle_handling=self.options.blq_middle_handling,
            last_handling=self.options.blq_last_handling,
        )

        if len(conc_clean) < 2:
            # Not enough data
            return results

        # Filter to interval
        mask = (time_clean >= interval.start) & (time_clean <= interval.end)
        conc_int = conc_clean[mask]
        time_int = time_clean[mask]

        if len(conc_int) < 2:
            return results

        # Calculate all requested parameters
        calc_results = self._calculate_all_parameters(
            conc=conc_int,
            time=time_int,
            dose=dose,
            duration=duration,
            method=self.options.auc_method,
        )

        # Build result list
        params_to_use = parameters or list(calc_results.keys())

        for param in params_to_use:
            if param in calc_results:
                results.append({
                    "subject": subject,
                    "parameter": param,
                    "value": calc_results[param],
                    "interval_start": interval.start,
                    "interval_end": interval.end,
                    "interval_name": interval.name,
                })

        return results

    def _calculate_all_parameters(
        self,
        conc: np.ndarray,
        time: np.ndarray,
        dose: Optional[float],
        duration: float = 0.0,
        method: str = "linear",
    ) -> Dict[str, float]:
        """Calculate all NCA parameters for a concentration-time profile."""
        results = {}

        # Basic parameters
        results["cmax"] = calc_cmax(conc, time)
        results["tmax"] = calc_tmax(conc, time)
        results["cmin"] = calc_cmin(conc, time)
        results["tlast"] = calc_tlast(conc, time, loq=self.options.loq)
        results["clast"] = calc_clast(conc, time, loq=self.options.loq)

        # AUC parameters
        results["auc.last"] = calc_auc_last(conc, time, method=method, loq=self.options.loq)
        results["auc.all"] = calc_auc_all(conc, time, method=method)

        # AUMC parameters
        results["aumc.last"] = calc_aumc_last(conc, time, method=method, loq=self.options.loq)

        # Half-life / lambda_z
        lambda_z_results = calc_lambda_z(
            conc,
            time,
            min_points=self.options.lambda_z_n_min,
            max_points=self.options.lambda_z_n_max,
            selection_method=self.options.lambda_z_selection,
            adj_r_squared_threshold=self.options.adj_r_squared_threshold,
        )

        results["lambda.z"] = lambda_z_results["lambda_z"]
        results["half.life"] = lambda_z_results["half_life"]
        results["r.squared"] = lambda_z_results["r_squared"]
        results["adj.r.squared"] = lambda_z_results["adj_r_squared"]
        results["lambda.z.n.points"] = lambda_z_results["n_points"]

        # Predicted Clast
        if not np.isnan(lambda_z_results["intercept"]) and not np.isnan(results["tlast"]):
            results["clast.pred"] = calc_clast_pred(
                lambda_z_results["intercept"],
                results["lambda.z"],
                results["tlast"],
            )
        else:
            results["clast.pred"] = np.nan

        # AUC to infinity (observed)
        results["auc.inf.obs"] = calc_auc_inf(
            conc, time,
            lambda_z=results["lambda.z"],
            method=method,
            loq=self.options.loq,
            clast_obs=results["clast"],
        )

        # AUC to infinity (predicted)
        results["auc.inf.pred"] = calc_auc_inf(
            conc, time,
            lambda_z=results["lambda.z"],
            method=method,
            loq=self.options.loq,
            clast_obs=results["clast.pred"],
        )

        # Percent extrapolated
        results["auc.pct.extrap.obs"] = calc_auc_pct_extrap(
            results["auc.last"], results["auc.inf.obs"]
        )
        results["auc.pct.extrap.pred"] = calc_auc_pct_extrap(
            results["auc.last"], results["auc.inf.pred"]
        )

        # AUMC to infinity
        results["aumc.inf.obs"] = calc_aumc_inf(
            conc, time,
            lambda_z=results["lambda.z"],
            method=method,
            loq=self.options.loq,
            clast_obs=results["clast"],
            tlast=results["tlast"],
        )

        results["aumc.inf.pred"] = calc_aumc_inf(
            conc, time,
            lambda_z=results["lambda.z"],
            method=method,
            loq=self.options.loq,
            clast_obs=results["clast.pred"],
            tlast=results["tlast"],
        )

        # MRT
        results["mrt.last"] = calc_mrt_last(results["aumc.last"], results["auc.last"])
        results["mrt.inf.obs"] = calc_mrt(results["aumc.inf.obs"], results["auc.inf.obs"])
        results["mrt.inf.pred"] = calc_mrt(results["aumc.inf.pred"], results["auc.inf.pred"])

        # Clearance and volume (require dose)
        if dose is not None and dose > 0:
            # Clearance
            results["cl.obs"] = calc_cl(dose, results["auc.inf.obs"])
            results["cl.pred"] = calc_cl(dose, results["auc.inf.pred"])
            results["cl.last"] = calc_cl(dose, results["auc.last"])

            # Volume of distribution
            results["vz.obs"] = calc_vz(results["cl.obs"], results["lambda.z"])
            results["vz.pred"] = calc_vz(results["cl.pred"], results["lambda.z"])

            # Vss
            results["vss.obs"] = calc_vss(
                dose, results["aumc.inf.obs"], results["auc.inf.obs"]
            )
            results["vss.pred"] = calc_vss(
                dose, results["aumc.inf.pred"], results["auc.inf.pred"]
            )
        else:
            results["cl.obs"] = np.nan
            results["cl.pred"] = np.nan
            results["cl.last"] = np.nan
            results["vz.obs"] = np.nan
            results["vz.pred"] = np.nan
            results["vss.obs"] = np.nan
            results["vss.pred"] = np.nan

        return results

    @property
    def results(self) -> Optional[NCAResults]:
        """Get the last calculated results."""
        return self._results


def run_nca(
    data: NCAData,
    parameters: Optional[List[str]] = None,
    options: Optional[NCAOptions] = None,
) -> NCAResults:
    """
    Run NCA analysis (convenience function).

    Parameters
    ----------
    data : NCAData
        Combined concentration-dose data
    parameters : list of str, optional
        Parameters to calculate
    options : NCAOptions, optional
        Analysis options

    Returns
    -------
    NCAResults
        Calculated NCA results

    Examples
    --------
    >>> results = nca.run_nca(data)
    >>> print(results.to_dataframe())
    """
    nca = NCA(data, options=options)
    return nca.calculate(parameters=parameters)
