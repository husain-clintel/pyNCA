"""Tests for multi-analyte NCA functionality."""

import numpy as np
import pandas as pd
import pytest

from pynca.analysis.multi_analyte import (
    calc_metabolite_ratio,
    calc_metabolite_ratios,
    molar_ratio,
    summarize_metabolite_ratios,
)


class TestCalcMetaboliteRatio:
    """Test calc_metabolite_ratio function."""

    def test_metabolite_to_parent(self):
        """Test M/P ratio calculation."""
        ratio = calc_metabolite_ratio(100.0, 50.0, ratio_type="metabolite_to_parent")
        assert ratio == 0.5

    def test_parent_to_metabolite(self):
        """Test P/M ratio calculation."""
        ratio = calc_metabolite_ratio(100.0, 50.0, ratio_type="parent_to_metabolite")
        assert ratio == 2.0

    def test_molar_ratio(self):
        """Test molar ratio calculation."""
        ratio = calc_metabolite_ratio(100.0, 50.0, ratio_type="molar")
        assert ratio == 0.5

    def test_zero_parent(self):
        """Test with zero parent value."""
        ratio = calc_metabolite_ratio(0.0, 50.0, ratio_type="metabolite_to_parent")
        assert np.isnan(ratio)

    def test_zero_metabolite(self):
        """Test with zero metabolite value."""
        ratio = calc_metabolite_ratio(100.0, 0.0, ratio_type="parent_to_metabolite")
        assert np.isnan(ratio)

    def test_nan_values(self):
        """Test with NaN values."""
        ratio = calc_metabolite_ratio(np.nan, 50.0)
        assert np.isnan(ratio)

        ratio = calc_metabolite_ratio(100.0, np.nan)
        assert np.isnan(ratio)

    def test_unknown_ratio_type(self):
        """Test with unknown ratio type."""
        with pytest.raises(ValueError, match="Unknown ratio_type"):
            calc_metabolite_ratio(100.0, 50.0, ratio_type="unknown")


class TestMolarRatio:
    """Test molar_ratio function."""

    def test_basic(self):
        """Test basic molar ratio conversion."""
        # Mass ratio = 0.5, parent MW = 200, metabolite MW = 100
        # Molar ratio = 0.5 * (200/100) = 1.0
        result = molar_ratio(0.5, 200.0, 100.0)
        assert result == 1.0

    def test_different_weights(self):
        """Test with different molecular weights."""
        result = molar_ratio(1.0, 300.0, 150.0)
        assert result == 2.0

    def test_nan_mass_ratio(self):
        """Test with NaN mass ratio."""
        result = molar_ratio(np.nan, 200.0, 100.0)
        assert np.isnan(result)

    def test_invalid_mw(self):
        """Test with invalid molecular weights."""
        result = molar_ratio(0.5, 0.0, 100.0)
        assert np.isnan(result)

        result = molar_ratio(0.5, 200.0, -100.0)
        assert np.isnan(result)


class TestCalcMetaboliteRatios:
    """Test calc_metabolite_ratios function."""

    def test_basic(self):
        """Test basic metabolite ratios calculation."""
        from pynca.core.results import NCAResults

        # Create mock parent results as list of dicts
        parent_results_list = [
            {"subject": 1, "parameter": "cmax", "value": 10.0, "interval_start": 0, "interval_end": 24},
            {"subject": 1, "parameter": "auc.last", "value": 100.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "parameter": "cmax", "value": 12.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "parameter": "auc.last", "value": 120.0, "interval_start": 0, "interval_end": 24},
        ]
        parent_results = NCAResults(results=parent_results_list)

        # Create mock metabolite results
        metabolite_results_list = [
            {"subject": 1, "parameter": "cmax", "value": 5.0, "interval_start": 0, "interval_end": 24},
            {"subject": 1, "parameter": "auc.last", "value": 50.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "parameter": "cmax", "value": 6.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "parameter": "auc.last", "value": 60.0, "interval_start": 0, "interval_end": 24},
        ]
        metabolite_results = NCAResults(results=metabolite_results_list)

        ratios = calc_metabolite_ratios(
            parent_results,
            metabolite_results,
            parameters=["cmax", "auc.last"],
        )

        assert len(ratios) == 4  # 2 subjects x 2 parameters
        assert all(ratios["ratio"] == 0.5)


class TestSummarizeMetaboliteRatios:
    """Test summarize_metabolite_ratios function."""

    def test_basic(self):
        """Test basic ratio summarization."""
        ratios_df = pd.DataFrame({
            "subject": [1, 2, 3, 1, 2, 3],
            "parameter": ["cmax", "cmax", "cmax", "auc.last", "auc.last", "auc.last"],
            "ratio": [0.5, 0.6, 0.55, 0.4, 0.45, 0.42],
        })

        summary = summarize_metabolite_ratios(ratios_df)

        assert len(summary) == 2  # 2 parameters
        assert "cmax_ratio" in summary["parameter"].values
        assert "auc.last_ratio" in summary["parameter"].values
        assert "mean" in summary.columns
        assert "n" in summary.columns
