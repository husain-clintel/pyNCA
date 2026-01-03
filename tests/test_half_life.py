"""Tests for half-life calculations."""

import numpy as np
import pytest
from pynca.calc.half_life import calc_lambda_z, calc_half_life


class TestCalcLambdaZ:
    """Tests for lambda_z calculation."""

    def test_monoexponential_decay(self):
        """Test with perfect monoexponential decay."""
        # Generate perfect monoexponential data
        lambda_z_true = 0.1
        t = np.array([0, 1, 2, 4, 6, 8, 10, 12])
        # Cmax at t=1, then exponential decay
        c = np.array([0, 10, 10 * np.exp(-lambda_z_true * 1),
                      10 * np.exp(-lambda_z_true * 3),
                      10 * np.exp(-lambda_z_true * 5),
                      10 * np.exp(-lambda_z_true * 7),
                      10 * np.exp(-lambda_z_true * 9),
                      10 * np.exp(-lambda_z_true * 11)])

        result = calc_lambda_z(c, t, min_points=3)

        # Should recover lambda_z closely
        assert not np.isnan(result["lambda_z"])
        assert np.isclose(result["lambda_z"], lambda_z_true, rtol=0.1)
        assert result["r_squared"] > 0.99

    def test_half_life_calculation(self):
        """Test half-life calculation."""
        lambda_z = 0.1
        t = np.array([0, 1, 2, 4, 6, 8])
        c = np.array([0, 10, 9, 7.4, 6.0, 4.9])

        result = calc_lambda_z(c, t)

        if not np.isnan(result["lambda_z"]):
            expected_half_life = np.log(2) / result["lambda_z"]
            assert np.isclose(result["half_life"], expected_half_life)

    def test_insufficient_points(self):
        """Test with too few points."""
        c = np.array([10, 5])
        t = np.array([0, 1])

        result = calc_lambda_z(c, t, min_points=3)
        assert np.isnan(result["lambda_z"])

    def test_all_same_concentration(self):
        """Test with constant concentration."""
        c = np.array([10, 10, 10, 10])
        t = np.array([0, 1, 2, 3])

        result = calc_lambda_z(c, t)
        # Lambda_z should be ~0 or NaN
        assert result["lambda_z"] <= 0 or np.isnan(result["lambda_z"])

    def test_best_fit_selection(self):
        """Test automatic point selection."""
        # Data with some noise
        np.random.seed(42)
        lambda_z_true = 0.15
        t = np.array([0, 0.5, 1, 2, 4, 6, 8, 10, 12, 24])
        c_true = 10 * np.exp(-lambda_z_true * t)
        c_true[0] = 0  # Pre-dose
        c_true[1] = 5  # Absorption phase
        c = c_true + np.random.normal(0, 0.1, len(c_true))
        c = np.maximum(c, 0.01)  # Keep positive

        result = calc_lambda_z(c, t, selection_method="best_fit")

        assert result["n_points"] >= 3
        assert result["adj_r_squared"] > 0


class TestCalcHalfLife:
    """Tests for calc_half_life wrapper."""

    def test_returns_same_as_lambda_z(self):
        """Test that calc_half_life returns same results as calc_lambda_z."""
        t = np.array([0, 1, 2, 4, 6, 8])
        c = np.array([0, 10, 8, 5.5, 3.8, 2.6])

        result_hl = calc_half_life(c, t)
        result_lz = calc_lambda_z(c, t)

        assert result_hl["lambda_z"] == result_lz["lambda_z"]
        assert result_hl["half_life"] == result_lz["half_life"]
