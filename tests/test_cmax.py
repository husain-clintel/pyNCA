"""Tests for Cmax and related parameter calculations."""

import numpy as np
import pytest
from pynca.calc.cmax import (
    calc_cmax,
    calc_tmax,
    calc_cmin,
    calc_tlast,
    calc_clast,
    calc_cav,
    calc_swing,
    calc_ptf,
)


class TestCalcCmax:
    """Tests for Cmax calculation."""

    def test_basic(self):
        """Test basic Cmax calculation."""
        conc = np.array([0, 5, 10, 8, 6, 4])
        time = np.array([0, 1, 2, 3, 4, 5])
        assert calc_cmax(conc, time) == 10

    def test_multiple_maxima(self):
        """Test with multiple equal maxima."""
        conc = np.array([0, 10, 5, 10, 5])
        time = np.array([0, 1, 2, 3, 4])
        assert calc_cmax(conc, time) == 10

    def test_empty(self):
        """Test with empty array raises error."""
        import pytest
        with pytest.raises(ValueError):
            calc_cmax(np.array([]), np.array([]))

    def test_with_nan(self):
        """Test with NaN values."""
        conc = np.array([0, np.nan, 10, 5])
        time = np.array([0, 1, 2, 3])
        assert calc_cmax(conc, time) == 10


class TestCalcTmax:
    """Tests for Tmax calculation."""

    def test_basic(self):
        """Test basic Tmax calculation."""
        conc = np.array([0, 5, 10, 8, 6])
        time = np.array([0, 1, 2, 3, 4])
        assert calc_tmax(conc, time) == 2

    def test_first_max(self):
        """Test first=True with multiple maxima."""
        conc = np.array([0, 10, 5, 10, 5])
        time = np.array([0, 1, 2, 3, 4])
        assert calc_tmax(conc, time, first=True) == 1

    def test_last_max(self):
        """Test first=False with multiple maxima."""
        conc = np.array([0, 10, 5, 10, 5])
        time = np.array([0, 1, 2, 3, 4])
        assert calc_tmax(conc, time, first=False) == 3


class TestCalcCmin:
    """Tests for Cmin calculation."""

    def test_basic(self):
        """Test basic Cmin calculation."""
        conc = np.array([0, 5, 10, 8, 2])
        time = np.array([0, 1, 2, 3, 4])
        # Excluding zeros, min is 2
        assert calc_cmin(conc, time, exclude_zero=True) == 2

    def test_include_zero(self):
        """Test with zero included."""
        conc = np.array([0, 5, 10, 8, 2])
        time = np.array([0, 1, 2, 3, 4])
        assert calc_cmin(conc, time, exclude_zero=False) == 0


class TestCalcTlast:
    """Tests for Tlast calculation."""

    def test_basic(self):
        """Test basic Tlast calculation."""
        conc = np.array([0, 10, 5, 2, 0])
        time = np.array([0, 1, 2, 3, 4])
        # Last measurable (>0) is at t=3
        assert calc_tlast(conc, time) == 3

    def test_all_measurable(self):
        """Test when all points are measurable."""
        conc = np.array([1, 10, 5, 2, 1])
        time = np.array([0, 1, 2, 3, 4])
        assert calc_tlast(conc, time) == 4


class TestCalcClast:
    """Tests for Clast calculation."""

    def test_basic(self):
        """Test basic Clast calculation."""
        conc = np.array([0, 10, 5, 2, 0])
        time = np.array([0, 1, 2, 3, 4])
        assert calc_clast(conc, time) == 2


class TestCalcCav:
    """Tests for Cav calculation."""

    def test_basic(self):
        """Test basic Cav calculation."""
        auc = 100  # concentration*time
        tau = 10   # time
        assert calc_cav(auc, tau) == 10

    def test_invalid_tau(self):
        """Test with invalid tau."""
        assert np.isnan(calc_cav(100, 0))


class TestCalcSwing:
    """Tests for swing calculation."""

    def test_basic(self):
        """Test swing calculation."""
        cmax = 10
        cmin = 2
        # Swing = (10 - 2) / 2 = 4
        assert np.isclose(calc_swing(cmax, cmin), 4.0)

    def test_invalid_cmin(self):
        """Test with invalid Cmin."""
        assert np.isnan(calc_swing(10, 0))


class TestCalcPTF:
    """Tests for PTF calculation."""

    def test_basic(self):
        """Test PTF calculation."""
        cmax = 10
        cmin = 2
        cav = 5
        # PTF = (10 - 2) / 5 * 100 = 160%
        assert np.isclose(calc_ptf(cmax, cmin, cav), 160.0)
