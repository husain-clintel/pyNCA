"""BLQ (Below Limit of Quantification) handling."""

from typing import Optional, Tuple, Literal
import numpy as np

from pynca.utils.helpers import is_blq, find_measurable_range


def clean_blq(
    conc: np.ndarray,
    time: np.ndarray,
    method: str = "standard",
    loq: Optional[float] = None,
    first_handling: str = "zero",
    middle_handling: str = "zero",
    last_handling: str = "drop",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Handle Below Limit of Quantification (BLQ) values.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    method : str
        Overall handling method:
        - "standard": Use first/middle/last handling rules
        - "zero": Replace all BLQ with 0
        - "drop": Remove all BLQ values
        - "loq/2": Replace all BLQ with LOQ/2
    loq : float, optional
        Limit of quantification
    first_handling : str
        How to handle BLQ before first measurable:
        - "zero": Replace with 0
        - "drop": Remove
    middle_handling : str
        How to handle BLQ between measurable values:
        - "zero": Replace with 0
        - "drop": Remove
        - "interpolate": Interpolate from surrounding values
    last_handling : str
        How to handle BLQ after last measurable:
        - "zero": Replace with 0
        - "drop": Remove

    Returns
    -------
    tuple
        (cleaned_conc, cleaned_time)
    """
    conc = np.asarray(conc, dtype=float)
    time = np.asarray(time, dtype=float)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    if method != "standard":
        return _simple_blq_handling(conc, time, method, loq)

    return _standard_blq_handling(
        conc, time, loq, first_handling, middle_handling, last_handling
    )


def _simple_blq_handling(
    conc: np.ndarray,
    time: np.ndarray,
    method: str,
    loq: Optional[float],
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply simple BLQ handling to all values."""
    blq_mask = is_blq(conc, loq)

    if method == "zero":
        result_conc = conc.copy()
        result_conc[blq_mask] = 0.0
        return result_conc, time

    elif method == "drop":
        keep_mask = ~blq_mask
        return conc[keep_mask], time[keep_mask]

    elif method == "loq/2":
        if loq is None:
            raise ValueError("LOQ must be specified for 'loq/2' method")
        result_conc = conc.copy()
        result_conc[blq_mask] = loq / 2
        return result_conc, time

    else:
        raise ValueError(f"Unknown BLQ method: {method}")


def _standard_blq_handling(
    conc: np.ndarray,
    time: np.ndarray,
    loq: Optional[float],
    first_handling: str,
    middle_handling: str,
    last_handling: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply standard BLQ handling with different rules for position."""
    first_idx, last_idx = find_measurable_range(conc, time, loq)

    if first_idx is None:
        # No measurable concentrations
        return np.array([]), np.array([])

    blq_mask = is_blq(conc, loq)
    result_conc = conc.copy()
    keep_mask = np.ones(len(conc), dtype=bool)

    for i in range(len(conc)):
        if not blq_mask[i]:
            continue

        if i < first_idx:
            # Before first measurable
            if first_handling == "zero":
                result_conc[i] = 0.0
            elif first_handling == "drop":
                keep_mask[i] = False
            else:
                raise ValueError(f"Unknown first_handling: {first_handling}")

        elif i > last_idx:
            # After last measurable
            if last_handling == "zero":
                result_conc[i] = 0.0
            elif last_handling == "drop":
                keep_mask[i] = False
            else:
                raise ValueError(f"Unknown last_handling: {last_handling}")

        else:
            # Between measurable values
            if middle_handling == "zero":
                result_conc[i] = 0.0
            elif middle_handling == "drop":
                keep_mask[i] = False
            elif middle_handling == "interpolate":
                # Find surrounding measurable values
                prev_idx = _find_prev_measurable(conc, i, blq_mask)
                next_idx = _find_next_measurable(conc, i, blq_mask)
                if prev_idx is not None and next_idx is not None:
                    result_conc[i] = _linear_interpolate(
                        time[prev_idx],
                        conc[prev_idx],
                        time[next_idx],
                        conc[next_idx],
                        time[i],
                    )
                else:
                    result_conc[i] = 0.0
            else:
                raise ValueError(f"Unknown middle_handling: {middle_handling}")

    return result_conc[keep_mask], time[keep_mask]


def _find_prev_measurable(
    conc: np.ndarray,
    idx: int,
    blq_mask: np.ndarray,
) -> Optional[int]:
    """Find index of previous measurable concentration."""
    for i in range(idx - 1, -1, -1):
        if not blq_mask[i]:
            return i
    return None


def _find_next_measurable(
    conc: np.ndarray,
    idx: int,
    blq_mask: np.ndarray,
) -> Optional[int]:
    """Find index of next measurable concentration."""
    for i in range(idx + 1, len(conc)):
        if not blq_mask[i]:
            return i
    return None


def _linear_interpolate(
    t1: float,
    c1: float,
    t2: float,
    c2: float,
    target_t: float,
) -> float:
    """Linear interpolation between two points."""
    if t2 == t1:
        return (c1 + c2) / 2
    return c1 + (c2 - c1) * (target_t - t1) / (t2 - t1)


def handle_blq(
    conc: np.ndarray,
    time: np.ndarray,
    loq: Optional[float] = None,
    method: str = "standard",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convenience function for BLQ handling with default settings.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    loq : float, optional
        Limit of quantification
    method : str
        Handling method

    Returns
    -------
    tuple
        (cleaned_conc, cleaned_time)
    """
    return clean_blq(conc, time, method=method, loq=loq)


def count_blq(
    conc: np.ndarray,
    loq: Optional[float] = None,
) -> int:
    """
    Count number of BLQ values.

    Parameters
    ----------
    conc : array-like
        Concentration values
    loq : float, optional
        Limit of quantification

    Returns
    -------
    int
        Number of BLQ values
    """
    blq_mask = is_blq(conc, loq)
    return int(np.sum(blq_mask))


def blq_summary(
    conc: np.ndarray,
    time: np.ndarray,
    loq: Optional[float] = None,
) -> dict:
    """
    Get summary of BLQ values.

    Parameters
    ----------
    conc : array-like
        Concentration values
    time : array-like
        Time values
    loq : float, optional
        Limit of quantification

    Returns
    -------
    dict
        Summary with counts and positions
    """
    conc = np.asarray(conc)
    time = np.asarray(time)

    # Sort by time
    sort_idx = np.argsort(time)
    conc = conc[sort_idx]
    time = time[sort_idx]

    blq_mask = is_blq(conc, loq)
    first_idx, last_idx = find_measurable_range(conc, time, loq)

    n_total = len(conc)
    n_blq = int(np.sum(blq_mask))

    n_before = 0
    n_middle = 0
    n_after = 0

    if first_idx is not None:
        for i in range(len(conc)):
            if blq_mask[i]:
                if i < first_idx:
                    n_before += 1
                elif i > last_idx:
                    n_after += 1
                else:
                    n_middle += 1

    return {
        "n_total": n_total,
        "n_blq": n_blq,
        "n_measurable": n_total - n_blq,
        "n_blq_before_first": n_before,
        "n_blq_middle": n_middle,
        "n_blq_after_last": n_after,
        "pct_blq": (n_blq / n_total * 100) if n_total > 0 else 0,
    }
