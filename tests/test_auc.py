"""Tests for AUC calculations."""

import numpy as np
import pytest
from pynca.calc.auc import (
    calc_auc,
    calc_auc_last,
    calc_auc_inf,
    calc_auc_pct_extrap,
)


class TestCalcAUC:
    """Tests for calc_auc function."""

    def test_linear_simple(self):
        """Test linear AUC with simple trapezoid."""
        conc = np.array([0, 10, 10, 0])
        time = np.array([0, 1, 2, 3])
        # Triangle + rectangle + triangle = 5 + 10 + 5 = 20
        auc = calc_auc(conc, time, method="linear")
        assert np.isclose(auc, 20.0)

    def test_linear_triangle(self):
        """Test linear AUC with triangle."""
        conc = np.array([0, 10, 0])
        time = np.array([0, 1, 2])
        # Two triangles = 5 + 5 = 10
        auc = calc_auc(conc, time, method="linear")
        assert np.isclose(auc, 10.0)

    def test_empty_array(self):
        """Test with empty arrays raises error."""
        import pytest
        with pytest.raises(ValueError):
            calc_auc(np.array([]), np.array([])[:0], method="linear")

    def test_single_point(self):
        """Test with single point."""
        auc = calc_auc(np.array([10]), np.array([1]), method="linear")
        assert auc == 0.0

    def test_log_method(self):
        """Test log-linear AUC."""
        conc = np.array([10, 5])
        time = np.array([0, 1])
        # Log trapezoid: (10 - 5) / ln(10/5) = 5 / 0.693 = 7.21
        auc_log = calc_auc(conc, time, method="log")
        auc_linear = calc_auc(conc, time, method="linear")
        # Linear = (10 + 5) / 2 = 7.5
        # Log = 7.21 (slightly lower for declining phase)
        assert np.isclose(auc_linear, 7.5)
        assert np.isclose(auc_log, 7.213, rtol=0.01)

    def test_linear_up_log_down(self):
        """Test linear-up/log-down method."""
        conc = np.array([0, 10, 5])
        time = np.array([0, 1, 2])
        auc = calc_auc(conc, time, method="linear-up/log-down")
        # First segment (up): linear = 5
        # Second segment (down): log
        assert auc > 0

    def test_interval(self):
        """Test AUC over interval."""
        conc = np.array([0, 10, 10, 5])
        time = np.array([0, 1, 2, 3])
        auc_full = calc_auc(conc, time, method="linear")
        auc_partial = calc_auc(conc, time, method="linear", start=1, end=2)
        assert auc_partial == 10.0
        assert auc_full > auc_partial


class TestCalcAUCLast:
    """Tests for calc_auc_last function."""

    def test_basic(self):
        """Test basic AUClast calculation."""
        conc = np.array([0, 10, 8, 6, 4, 2, 0])
        time = np.array([0, 1, 2, 3, 4, 5, 6])
        auc_last = calc_auc_last(conc, time, method="linear")
        # Last measurable is at t=5 (conc=2)
        assert auc_last > 0

    def test_all_blq(self):
        """Test with all BLQ values."""
        conc = np.array([0, 0, 0])
        time = np.array([0, 1, 2])
        auc_last = calc_auc_last(conc, time)
        assert np.isnan(auc_last)


class TestCalcAUCInf:
    """Tests for calc_auc_inf function."""

    def test_basic(self):
        """Test AUCinf calculation."""
        conc = np.array([0, 10, 8, 6, 4, 2])
        time = np.array([0, 1, 2, 3, 4, 5])
        lambda_z = 0.2  # Estimated rate constant
        auc_inf = calc_auc_inf(conc, time, lambda_z=lambda_z)
        auc_last = calc_auc_last(conc, time)
        # AUCinf should be greater than AUClast
        assert auc_inf > auc_last

    def test_invalid_lambda_z(self):
        """Test with invalid lambda_z."""
        conc = np.array([0, 10, 5])
        time = np.array([0, 1, 2])
        assert np.isnan(calc_auc_inf(conc, time, lambda_z=None))
        assert np.isnan(calc_auc_inf(conc, time, lambda_z=np.nan))
        assert np.isnan(calc_auc_inf(conc, time, lambda_z=-0.1))


class TestCalcAUCPctExtrap:
    """Tests for calc_auc_pct_extrap function."""

    def test_basic(self):
        """Test percent extrapolated."""
        auc_last = 80
        auc_inf = 100
        pct = calc_auc_pct_extrap(auc_last, auc_inf)
        assert np.isclose(pct, 20.0)

    def test_no_extrapolation(self):
        """Test with no extrapolation needed."""
        pct = calc_auc_pct_extrap(100, 100)
        assert np.isclose(pct, 0.0)

    def test_invalid_inputs(self):
        """Test with invalid inputs."""
        assert np.isnan(calc_auc_pct_extrap(np.nan, 100))
        assert np.isnan(calc_auc_pct_extrap(100, np.nan))
        assert np.isnan(calc_auc_pct_extrap(100, 0))
