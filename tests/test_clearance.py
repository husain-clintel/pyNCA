"""Tests for clearance and volume calculations."""

import numpy as np
import pytest
from pynca.calc.clearance import (
    calc_cl,
    calc_vz,
    calc_vss,
    calc_mrt,
    calc_mrt_last,
)


class TestCalcCL:
    """Tests for clearance calculation."""

    def test_basic(self):
        """Test basic clearance calculation."""
        dose = 100  # mg
        auc_inf = 50  # mg*h/L
        cl = calc_cl(dose, auc_inf)
        assert np.isclose(cl, 2.0)  # L/h

    def test_with_bioavailability(self):
        """Test clearance with bioavailability factor."""
        dose = 100
        auc_inf = 50
        f = 0.5
        cl = calc_cl(dose, auc_inf, f=f)
        assert np.isclose(cl, 1.0)  # CL * F

    def test_invalid_inputs(self):
        """Test with invalid inputs."""
        assert np.isnan(calc_cl(np.nan, 50))
        assert np.isnan(calc_cl(100, np.nan))
        assert np.isnan(calc_cl(100, 0))
        assert np.isnan(calc_cl(0, 50))


class TestCalcVz:
    """Tests for volume of distribution calculation."""

    def test_basic(self):
        """Test Vz calculation."""
        cl = 2.0  # L/h
        lambda_z = 0.1  # 1/h
        vz = calc_vz(cl, lambda_z)
        assert np.isclose(vz, 20.0)  # L

    def test_invalid_lambda_z(self):
        """Test with invalid lambda_z."""
        assert np.isnan(calc_vz(2.0, 0))
        assert np.isnan(calc_vz(2.0, -0.1))
        assert np.isnan(calc_vz(2.0, np.nan))


class TestCalcVss:
    """Tests for Vss calculation."""

    def test_basic(self):
        """Test Vss calculation."""
        dose = 100
        aumc_inf = 5000
        auc_inf = 50
        vss = calc_vss(dose, aumc_inf, auc_inf)
        # Vss = dose * AUMC / AUC^2 = 100 * 5000 / 2500 = 200
        assert np.isclose(vss, 200.0)

    def test_invalid_inputs(self):
        """Test with invalid inputs."""
        assert np.isnan(calc_vss(100, 5000, 0))
        assert np.isnan(calc_vss(0, 5000, 50))


class TestCalcMRT:
    """Tests for MRT calculation."""

    def test_basic(self):
        """Test MRT calculation."""
        aumc_inf = 500
        auc_inf = 50
        mrt = calc_mrt(aumc_inf, auc_inf)
        assert np.isclose(mrt, 10.0)  # h

    def test_invalid_inputs(self):
        """Test with invalid inputs."""
        assert np.isnan(calc_mrt(np.nan, 50))
        assert np.isnan(calc_mrt(500, 0))


class TestCalcMRTLast:
    """Tests for MRT_last calculation."""

    def test_basic(self):
        """Test MRTlast calculation."""
        aumc_last = 400
        auc_last = 40
        mrt = calc_mrt_last(aumc_last, auc_last)
        assert np.isclose(mrt, 10.0)
