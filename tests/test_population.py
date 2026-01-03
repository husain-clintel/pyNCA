"""Tests for population NCA summary functions."""

import numpy as np
import pandas as pd
import pytest

from pynca.summary.summarize import (
    summarize_by_group,
    compare_groups,
    bioequivalence_analysis,
    population_pk_summary,
    inter_subject_variability,
    dose_proportionality,
)
from pynca.core.results import NCAResults


def create_results(data_list):
    """Helper to create NCAResults from a list of dicts."""
    return NCAResults(results=data_list)


class TestSummarizeByGroup:
    """Test summarize_by_group function."""

    def test_basic(self):
        """Test basic group summary."""
        data = [
            {"subject": 1, "treatment": "A", "parameter": "cmax", "value": 10.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "treatment": "A", "parameter": "cmax", "value": 12.0, "interval_start": 0, "interval_end": 24},
            {"subject": 3, "treatment": "B", "parameter": "cmax", "value": 11.0, "interval_start": 0, "interval_end": 24},
            {"subject": 4, "treatment": "B", "parameter": "cmax", "value": 13.0, "interval_start": 0, "interval_end": 24},
            {"subject": 1, "treatment": "A", "parameter": "auc.last", "value": 100.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "treatment": "A", "parameter": "auc.last", "value": 120.0, "interval_start": 0, "interval_end": 24},
            {"subject": 3, "treatment": "B", "parameter": "auc.last", "value": 110.0, "interval_start": 0, "interval_end": 24},
            {"subject": 4, "treatment": "B", "parameter": "auc.last", "value": 130.0, "interval_start": 0, "interval_end": 24},
        ]
        results = create_results(data)

        summary = summarize_by_group(results, group_col="treatment")

        assert len(summary) == 4  # 2 groups x 2 parameters
        assert "treatment" in summary.columns
        assert "n" in summary.columns
        assert "mean" in summary.columns

    def test_missing_group_col(self):
        """Test with missing group column."""
        data = [
            {"subject": 1, "parameter": "cmax", "value": 10.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "parameter": "cmax", "value": 12.0, "interval_start": 0, "interval_end": 24},
        ]
        results = create_results(data)

        with pytest.raises(ValueError, match="not found in results"):
            summarize_by_group(results, group_col="treatment")


class TestCompareGroups:
    """Test compare_groups function."""

    def test_basic(self):
        """Test basic group comparison."""
        data = [
            {"subject": 1, "treatment": "test", "parameter": "cmax", "value": 10.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "treatment": "test", "parameter": "cmax", "value": 12.0, "interval_start": 0, "interval_end": 24},
            {"subject": 3, "treatment": "test", "parameter": "cmax", "value": 11.0, "interval_start": 0, "interval_end": 24},
            {"subject": 4, "treatment": "ref", "parameter": "cmax", "value": 9.0, "interval_start": 0, "interval_end": 24},
            {"subject": 5, "treatment": "ref", "parameter": "cmax", "value": 10.0, "interval_start": 0, "interval_end": 24},
            {"subject": 6, "treatment": "ref", "parameter": "cmax", "value": 11.0, "interval_start": 0, "interval_end": 24},
        ]
        results = create_results(data)

        comparison = compare_groups(
            results,
            group_col="treatment",
            test_group="test",
            ref_group="ref",
            parameters=["cmax"],
        )

        assert len(comparison) == 1
        assert "ratio" in comparison.columns
        assert "lower" in comparison.columns
        assert "upper" in comparison.columns

    def test_arithmetic_ratio(self):
        """Test with arithmetic mean ratio."""
        data = [
            {"subject": 1, "treatment": "test", "parameter": "cmax", "value": 12.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "treatment": "test", "parameter": "cmax", "value": 12.0, "interval_start": 0, "interval_end": 24},
            {"subject": 3, "treatment": "ref", "parameter": "cmax", "value": 10.0, "interval_start": 0, "interval_end": 24},
            {"subject": 4, "treatment": "ref", "parameter": "cmax", "value": 10.0, "interval_start": 0, "interval_end": 24},
        ]
        results = create_results(data)

        comparison = compare_groups(
            results,
            group_col="treatment",
            test_group="test",
            ref_group="ref",
            log_transform=False,
        )

        assert comparison.iloc[0]["ratio"] == 1.2


class TestBioequivalenceAnalysis:
    """Test bioequivalence_analysis function."""

    def test_basic(self):
        """Test basic BE analysis."""
        np.random.seed(42)
        values = np.random.lognormal(2, 0.2, 24)
        data = []
        for i in range(12):
            data.append({
                "subject": i + 1,
                "formulation": "test",
                "parameter": "cmax",
                "value": values[i],
                "interval_start": 0,
                "interval_end": 24,
            })
        for i in range(12, 24):
            data.append({
                "subject": i + 1,
                "formulation": "reference",
                "parameter": "cmax",
                "value": values[i],
                "interval_start": 0,
                "interval_end": 24,
            })

        results = create_results(data)

        be_results = bioequivalence_analysis(
            results,
            group_col="formulation",
            test_group="test",
            ref_group="reference",
            parameters=["cmax"],
        )

        assert len(be_results) == 1
        assert "be_conclusion" in be_results.columns
        assert "within_limits" in be_results.columns
        assert "be_lower_limit" in be_results.columns
        assert "be_upper_limit" in be_results.columns


class TestPopulationPKSummary:
    """Test population_pk_summary function."""

    def test_basic(self):
        """Test basic population summary."""
        data = [
            {"subject": 1, "parameter": "cmax", "value": 10.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "parameter": "cmax", "value": 12.0, "interval_start": 0, "interval_end": 24},
            {"subject": 3, "parameter": "cmax", "value": 11.0, "interval_start": 0, "interval_end": 24},
            {"subject": 4, "parameter": "cmax", "value": 13.0, "interval_start": 0, "interval_end": 24},
            {"subject": 5, "parameter": "cmax", "value": 14.0, "interval_start": 0, "interval_end": 24},
        ]
        results = create_results(data)

        summary = population_pk_summary(results)

        assert len(summary) == 1
        assert "mean" in summary.columns
        assert "geomean" in summary.columns
        assert "ci95_lower" in summary.columns
        assert "geo_ci90_lower" in summary.columns


class TestInterSubjectVariability:
    """Test inter_subject_variability function."""

    def test_basic(self):
        """Test basic variability calculation."""
        data = [
            {"subject": 1, "parameter": "cmax", "value": 10.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "parameter": "cmax", "value": 12.0, "interval_start": 0, "interval_end": 24},
            {"subject": 3, "parameter": "cmax", "value": 11.0, "interval_start": 0, "interval_end": 24},
            {"subject": 4, "parameter": "cmax", "value": 13.0, "interval_start": 0, "interval_end": 24},
            {"subject": 5, "parameter": "cmax", "value": 14.0, "interval_start": 0, "interval_end": 24},
        ]
        results = create_results(data)

        isv = inter_subject_variability(results)

        assert len(isv) == 1
        assert "cv_pct" in isv.columns
        assert "geocv_pct" in isv.columns
        assert "variability_category" in isv.columns

    def test_variability_categories(self):
        """Test variability categorization."""
        # Low variability (CV < 30%)
        data_low = [
            {"subject": 1, "parameter": "cmax", "value": 10.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "parameter": "cmax", "value": 10.5, "interval_start": 0, "interval_end": 24},
            {"subject": 3, "parameter": "cmax", "value": 10.2, "interval_start": 0, "interval_end": 24},
            {"subject": 4, "parameter": "cmax", "value": 10.3, "interval_start": 0, "interval_end": 24},
            {"subject": 5, "parameter": "cmax", "value": 10.4, "interval_start": 0, "interval_end": 24},
        ]
        results_low = create_results(data_low)
        isv_low = inter_subject_variability(results_low)
        assert isv_low.iloc[0]["variability_category"] == "Low"

        # High variability (CV > 60%)
        data_high = [
            {"subject": 1, "parameter": "cmax", "value": 5.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "parameter": "cmax", "value": 10.0, "interval_start": 0, "interval_end": 24},
            {"subject": 3, "parameter": "cmax", "value": 15.0, "interval_start": 0, "interval_end": 24},
            {"subject": 4, "parameter": "cmax", "value": 20.0, "interval_start": 0, "interval_end": 24},
            {"subject": 5, "parameter": "cmax", "value": 25.0, "interval_start": 0, "interval_end": 24},
        ]
        results_high = create_results(data_high)
        isv_high = inter_subject_variability(results_high)
        assert isv_high.iloc[0]["variability_category"] == "High"


class TestDoseProportionality:
    """Test dose_proportionality function."""

    def test_power_model(self):
        """Test power model dose proportionality."""
        # Create mock results - dose proportional (beta = 1)
        results_10 = create_results([
            {"subject": 1, "parameter": "auc.last", "value": 100.0, "interval_start": 0, "interval_end": 24},
            {"subject": 2, "parameter": "auc.last", "value": 110.0, "interval_start": 0, "interval_end": 24},
            {"subject": 3, "parameter": "auc.last", "value": 105.0, "interval_start": 0, "interval_end": 24},
        ])
        results_20 = create_results([
            {"subject": 4, "parameter": "auc.last", "value": 200.0, "interval_start": 0, "interval_end": 24},
            {"subject": 5, "parameter": "auc.last", "value": 220.0, "interval_start": 0, "interval_end": 24},
            {"subject": 6, "parameter": "auc.last", "value": 210.0, "interval_start": 0, "interval_end": 24},
        ])
        results_40 = create_results([
            {"subject": 7, "parameter": "auc.last", "value": 400.0, "interval_start": 0, "interval_end": 24},
            {"subject": 8, "parameter": "auc.last", "value": 440.0, "interval_start": 0, "interval_end": 24},
            {"subject": 9, "parameter": "auc.last", "value": 420.0, "interval_start": 0, "interval_end": 24},
        ])

        dp = dose_proportionality(
            {10: results_10, 20: results_20, 40: results_40},
            parameter="auc.last",
            method="power_model",
        )

        assert "beta" in dp
        assert "is_proportional" in dp
        assert "conclusion" in dp
        assert dp["n_doses"] == 3

    def test_insufficient_data(self):
        """Test with insufficient data."""
        results = create_results([
            {"subject": 1, "parameter": "auc.last", "value": 100.0, "interval_start": 0, "interval_end": 24},
        ])

        dp = dose_proportionality(
            {10: results, 20: results},
            parameter="auc.last",
        )

        assert "error" in dp
