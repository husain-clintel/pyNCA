"""Parallel processing utilities for NCA calculations."""

from typing import Any, Callable, Dict, List, Optional, Union
import warnings
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing as mp

import numpy as np
import pandas as pd


def get_n_workers(n_jobs: int = -1) -> int:
    """
    Get number of workers for parallel processing.

    Parameters
    ----------
    n_jobs : int
        Number of parallel jobs.
        -1 = use all available CPUs
        -2 = use all but one CPU
        positive int = use that many workers

    Returns
    -------
    int
        Number of workers to use
    """
    max_workers = mp.cpu_count()

    if n_jobs == -1:
        return max_workers
    elif n_jobs == -2:
        return max(1, max_workers - 1)
    elif n_jobs > 0:
        return min(n_jobs, max_workers)
    else:
        return 1


def run_nca_parallel(
    data: "NCAData",
    n_jobs: int = -1,
    backend: str = "process",
    verbose: bool = False,
    **kwargs,
) -> "NCAResults":
    """
    Run NCA analysis in parallel across subjects.

    Parameters
    ----------
    data : NCAData
        NCA data object
    n_jobs : int
        Number of parallel jobs (-1 = all CPUs)
    backend : str
        Parallelization backend: "process" or "thread"
    verbose : bool
        Print progress information
    **kwargs
        Additional arguments passed to run_nca

    Returns
    -------
    NCAResults
        Combined NCA results from all subjects

    Examples
    --------
    >>> results = run_nca_parallel(data, n_jobs=4)
    """
    from pynca.analysis.nca import run_nca
    from pynca.core.data import NCAData
    from pynca.core.results import NCAResults

    # Get subjects
    subjects = data.conc.subjects

    if len(subjects) <= 1:
        # Single subject - no parallelization needed
        return run_nca(data, **kwargs)

    n_workers = get_n_workers(n_jobs)

    if n_workers == 1 or len(subjects) < 3:
        # Don't parallelize for very small datasets
        return run_nca(data, **kwargs)

    if verbose:
        print(f"Running NCA in parallel with {n_workers} workers for {len(subjects)} subjects")

    # Prepare subject data
    subject_data_list = []
    for subject in subjects:
        # Filter data for this subject
        subj_conc = data.conc.filter(subject=subject)
        subj_dose = data.dose.filter(subject=subject) if data.dose is not None else None

        subject_data_list.append({
            "subject": subject,
            "conc_data": subj_conc.data,
            "conc_col": subj_conc.conc_col,
            "time_col": subj_conc.time_col,
            "subject_col": subj_conc.subject_col,
            "dose_data": subj_dose.data if subj_dose else None,
            "dose_col": data.dose.dose_col if data.dose else None,
            "dose_time_col": data.dose.time_col if data.dose else None,
            "kwargs": kwargs,
        })

    # Select executor
    if backend == "process":
        Executor = ProcessPoolExecutor
    elif backend == "thread":
        Executor = ThreadPoolExecutor
    else:
        raise ValueError(f"Unknown backend: {backend}")

    # Run in parallel
    all_results = []

    try:
        with Executor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(_run_subject_nca, sd): sd["subject"]
                for sd in subject_data_list
            }

            for future in as_completed(futures):
                subject = futures[future]
                try:
                    result = future.result()
                    if result is not None:
                        all_results.append(result)
                    if verbose:
                        print(f"  Completed subject: {subject}")
                except Exception as e:
                    warnings.warn(f"Error processing subject {subject}: {e}")

    except Exception as e:
        warnings.warn(f"Parallel execution failed: {e}. Falling back to sequential.")
        return run_nca(data, **kwargs)

    # Combine results
    if not all_results:
        return NCAResults(data=pd.DataFrame())

    combined_df = pd.concat(all_results, ignore_index=True)
    return NCAResults(data=combined_df)


def _run_subject_nca(subject_data: dict) -> Optional[pd.DataFrame]:
    """
    Run NCA for a single subject (worker function).

    Parameters
    ----------
    subject_data : dict
        Dictionary with subject data and configuration

    Returns
    -------
    pd.DataFrame or None
        Results DataFrame for this subject
    """
    try:
        from pynca.core.concentration import NCAConcentration
        from pynca.core.dose import NCADose
        from pynca.core.data import NCAData
        from pynca.analysis.nca import run_nca

        # Reconstruct objects
        conc = NCAConcentration(
            data=subject_data["conc_data"],
            conc_col=subject_data["conc_col"],
            time_col=subject_data["time_col"],
            subject_col=subject_data["subject_col"],
        )

        dose = None
        if subject_data["dose_data"] is not None:
            dose = NCADose(
                data=subject_data["dose_data"],
                dose_col=subject_data["dose_col"],
                time_col=subject_data["dose_time_col"],
                subject_col=subject_data["subject_col"],
            )

        data = NCAData(conc=conc, dose=dose)
        results = run_nca(data, **subject_data["kwargs"])

        return results.data

    except Exception:
        return None


def parallel_map(
    func: Callable,
    items: List[Any],
    n_jobs: int = -1,
    backend: str = "thread",
    verbose: bool = False,
) -> List[Any]:
    """
    Apply function to items in parallel.

    Parameters
    ----------
    func : callable
        Function to apply to each item
    items : list
        Items to process
    n_jobs : int
        Number of parallel jobs
    backend : str
        Parallelization backend: "process" or "thread"
    verbose : bool
        Print progress

    Returns
    -------
    list
        Results from applying func to each item

    Examples
    --------
    >>> def square(x):
    ...     return x ** 2
    >>> results = parallel_map(square, [1, 2, 3, 4], n_jobs=2)
    >>> print(results)
    [1, 4, 9, 16]
    """
    n_workers = get_n_workers(n_jobs)

    if n_workers == 1 or len(items) < 3:
        return [func(item) for item in items]

    if backend == "process":
        Executor = ProcessPoolExecutor
    elif backend == "thread":
        Executor = ThreadPoolExecutor
    else:
        raise ValueError(f"Unknown backend: {backend}")

    results = [None] * len(items)

    with Executor(max_workers=n_workers) as executor:
        futures = {executor.submit(func, item): i for i, item in enumerate(items)}

        for future in as_completed(futures):
            idx = futures[future]
            try:
                results[idx] = future.result()
                if verbose:
                    print(f"  Completed item {idx + 1}/{len(items)}")
            except Exception as e:
                warnings.warn(f"Error processing item {idx}: {e}")
                results[idx] = None

    return results


def run_multi_analyte_nca_parallel(
    data: "NCAData",
    analytes: Optional[List[str]] = None,
    n_jobs: int = -1,
    verbose: bool = False,
    **kwargs,
) -> Dict[str, "NCAResults"]:
    """
    Run multi-analyte NCA in parallel.

    Parameters
    ----------
    data : NCAData
        NCA data with analyte column
    analytes : list of str, optional
        Analytes to process
    n_jobs : int
        Number of parallel jobs
    verbose : bool
        Print progress
    **kwargs
        Additional arguments for run_nca

    Returns
    -------
    dict
        Dictionary mapping analyte name to NCAResults

    Examples
    --------
    >>> results = run_multi_analyte_nca_parallel(data, n_jobs=4)
    """
    from pynca.analysis.nca import run_nca
    from pynca.core.data import NCAData

    if data.conc.analyte_col is None:
        return {"default": run_nca(data, **kwargs)}

    # Get analytes
    available_analytes = data.conc.data[data.conc.analyte_col].unique().tolist()

    if analytes is None:
        analytes = available_analytes
    else:
        missing = set(analytes) - set(available_analytes)
        if missing:
            raise ValueError(f"Analytes not found: {missing}")

    n_workers = get_n_workers(n_jobs)

    if n_workers == 1 or len(analytes) < 2:
        # Sequential processing
        from pynca.analysis.multi_analyte import run_multi_analyte_nca
        return run_multi_analyte_nca(data, analytes=analytes, **kwargs)

    if verbose:
        print(f"Processing {len(analytes)} analytes with {n_workers} workers")

    # Prepare analyte data
    analyte_configs = []
    for analyte in analytes:
        analyte_conc = data.conc.filter(**{data.conc.analyte_col: analyte})
        analyte_configs.append({
            "analyte": analyte,
            "conc_data": analyte_conc.data,
            "conc_col": analyte_conc.conc_col,
            "time_col": analyte_conc.time_col,
            "subject_col": analyte_conc.subject_col,
            "dose_data": data.dose.data if data.dose else None,
            "dose_col": data.dose.dose_col if data.dose else None,
            "dose_time_col": data.dose.time_col if data.dose else None,
            "kwargs": kwargs,
        })

    # Process in parallel (using threads for simplicity with shared memory)
    results = {}

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(_run_analyte_nca, cfg): cfg["analyte"]
            for cfg in analyte_configs
        }

        for future in as_completed(futures):
            analyte = futures[future]
            try:
                result = future.result()
                if result is not None:
                    results[analyte] = result
                if verbose:
                    print(f"  Completed analyte: {analyte}")
            except Exception as e:
                warnings.warn(f"Error processing analyte {analyte}: {e}")

    return results


def _run_analyte_nca(config: dict) -> Optional["NCAResults"]:
    """
    Run NCA for a single analyte (worker function).
    """
    try:
        from pynca.core.concentration import NCAConcentration
        from pynca.core.dose import NCADose
        from pynca.core.data import NCAData
        from pynca.analysis.nca import run_nca

        conc = NCAConcentration(
            data=config["conc_data"],
            conc_col=config["conc_col"],
            time_col=config["time_col"],
            subject_col=config["subject_col"],
        )

        dose = None
        if config["dose_data"] is not None:
            dose = NCADose(
                data=config["dose_data"],
                dose_col=config["dose_col"],
                time_col=config["dose_time_col"],
                subject_col=config["subject_col"],
            )

        data = NCAData(conc=conc, dose=dose)
        return run_nca(data, **config["kwargs"])

    except Exception:
        return None


class ParallelNCA:
    """
    Convenience class for parallel NCA operations.

    Parameters
    ----------
    n_jobs : int
        Number of parallel jobs (-1 = all CPUs)
    backend : str
        Parallelization backend: "process" or "thread"
    verbose : bool
        Print progress information

    Examples
    --------
    >>> pnca = ParallelNCA(n_jobs=4)
    >>> results = pnca.run(data)
    """

    def __init__(
        self,
        n_jobs: int = -1,
        backend: str = "thread",
        verbose: bool = False,
    ):
        self.n_jobs = n_jobs
        self.backend = backend
        self.verbose = verbose
        self._n_workers = get_n_workers(n_jobs)

    @property
    def n_workers(self) -> int:
        """Number of workers that will be used."""
        return self._n_workers

    def run(self, data: "NCAData", **kwargs) -> "NCAResults":
        """Run parallel NCA."""
        return run_nca_parallel(
            data,
            n_jobs=self.n_jobs,
            backend=self.backend,
            verbose=self.verbose,
            **kwargs,
        )

    def run_multi_analyte(
        self,
        data: "NCAData",
        analytes: Optional[List[str]] = None,
        **kwargs,
    ) -> Dict[str, "NCAResults"]:
        """Run parallel multi-analyte NCA."""
        return run_multi_analyte_nca_parallel(
            data,
            analytes=analytes,
            n_jobs=self.n_jobs,
            verbose=self.verbose,
            **kwargs,
        )

    def __repr__(self) -> str:
        return (
            f"ParallelNCA(n_jobs={self.n_jobs}, "
            f"n_workers={self.n_workers}, backend='{self.backend}')"
        )
