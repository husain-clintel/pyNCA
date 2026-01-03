"""Integration tests for complete NCA workflow."""

import os
import numpy as np
import pandas as pd
import pytest

# Get the path to test data
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


class TestFullWorkflow:
    """Integration tests for complete NCA analysis."""

    @pytest.fixture
    def theophylline_data(self):
        """Load theophylline test data."""
        csv_path = os.path.join(TEST_DATA_DIR, "theophylline.csv")
        return pd.read_csv(csv_path)

    def test_import_package(self):
        """Test that package imports correctly."""
        import pynca
        assert hasattr(pynca, "NCAConcentration")
        assert hasattr(pynca, "NCADose")
        assert hasattr(pynca, "NCAData")
        assert hasattr(pynca, "run_nca")

    def test_create_concentration_object(self, theophylline_data):
        """Test creating NCAConcentration object."""
        from pynca import NCAConcentration

        conc = NCAConcentration(
            data=theophylline_data,
            conc_col="conc",
            time_col="Time",
            subject_col="Subject",
        )

        assert conc.n_subjects == 6
        assert len(conc) == len(theophylline_data)

    def test_create_dose_object(self, theophylline_data):
        """Test creating NCADose object."""
        from pynca import NCADose

        dose_df = theophylline_data[["Subject", "Dose"]].drop_duplicates()
        dose_df["Time"] = 0

        dose = NCADose(
            data=dose_df,
            dose_col="Dose",
            time_col="Time",
            subject_col="Subject",
        )

        assert dose.n_subjects == 6

    def test_create_nca_data(self, theophylline_data):
        """Test creating NCAData object."""
        from pynca import NCAConcentration, NCADose, NCAData

        conc = NCAConcentration(
            data=theophylline_data,
            conc_col="conc",
            time_col="Time",
            subject_col="Subject",
        )

        dose_df = theophylline_data[["Subject", "Dose"]].drop_duplicates()
        dose_df["Time"] = 0

        dose = NCADose(
            data=dose_df,
            dose_col="Dose",
            time_col="Time",
            subject_col="Subject",
        )

        data = NCAData(conc=conc, dose=dose)
        assert data.n_subjects == 6

    def test_run_nca_analysis(self, theophylline_data):
        """Test running full NCA analysis."""
        from pynca import NCAConcentration, NCADose, NCAData, run_nca

        conc = NCAConcentration(
            data=theophylline_data,
            conc_col="conc",
            time_col="Time",
            subject_col="Subject",
        )

        dose_df = theophylline_data[["Subject", "Dose"]].drop_duplicates()
        dose_df["Time"] = 0

        dose = NCADose(
            data=dose_df,
            dose_col="Dose",
            time_col="Time",
            subject_col="Subject",
        )

        data = NCAData(conc=conc, dose=dose)
        results = run_nca(data)

        # Check that results were calculated
        assert len(results) > 0
        assert results.n_subjects == 6

        # Check that key parameters exist
        df = results.to_dataframe()
        params = df["parameter"].unique()
        assert "cmax" in params
        assert "tmax" in params
        assert "auc.last" in params

    def test_results_to_wide_format(self, theophylline_data):
        """Test converting results to wide format."""
        from pynca import NCAConcentration, NCADose, NCAData, run_nca

        conc = NCAConcentration(
            data=theophylline_data,
            conc_col="conc",
            time_col="Time",
            subject_col="Subject",
        )

        dose_df = theophylline_data[["Subject", "Dose"]].drop_duplicates()
        dose_df["Time"] = 0

        dose = NCADose(
            data=dose_df,
            dose_col="Dose",
            time_col="Time",
            subject_col="Subject",
        )

        data = NCAData(conc=conc, dose=dose)
        results = run_nca(data)

        df_wide = results.to_dataframe(wide=True)

        # Should have one row per subject
        assert len(df_wide) == 6
        # Should have parameter columns
        assert "cmax" in df_wide.columns

    def test_filter_results(self, theophylline_data):
        """Test filtering results."""
        from pynca import NCAConcentration, NCADose, NCAData, run_nca

        conc = NCAConcentration(
            data=theophylline_data,
            conc_col="conc",
            time_col="Time",
            subject_col="Subject",
        )

        dose_df = theophylline_data[["Subject", "Dose"]].drop_duplicates()
        dose_df["Time"] = 0

        dose = NCADose(
            data=dose_df,
            dose_col="Dose",
            time_col="Time",
            subject_col="Subject",
        )

        data = NCAData(conc=conc, dose=dose)
        results = run_nca(data)

        # Filter by parameter
        filtered = results.filter(parameters=["cmax", "tmax"])
        params = filtered.data["parameter"].unique()
        assert len(params) == 2
        assert "cmax" in params
        assert "tmax" in params

    def test_get_parameter_value(self, theophylline_data):
        """Test getting specific parameter values."""
        from pynca import NCAConcentration, NCADose, NCAData, run_nca

        conc = NCAConcentration(
            data=theophylline_data,
            conc_col="conc",
            time_col="Time",
            subject_col="Subject",
        )

        dose_df = theophylline_data[["Subject", "Dose"]].drop_duplicates()
        dose_df["Time"] = 0

        dose = NCADose(
            data=dose_df,
            dose_col="Dose",
            time_col="Time",
            subject_col="Subject",
        )

        data = NCAData(conc=conc, dose=dose)
        results = run_nca(data)

        cmax_values = results.get_parameter("cmax")
        assert len(cmax_values) == 6
        # All Cmax values should be positive
        assert all(cmax_values > 0)


class TestSingleSubject:
    """Tests for single-subject analysis."""

    def test_single_subject(self):
        """Test NCA on single subject."""
        from pynca import NCAConcentration, NCADose, NCAData, run_nca

        conc_df = pd.DataFrame({
            "time": [0, 0.5, 1, 2, 4, 6, 8, 12, 24],
            "conc": [0, 5, 10, 8, 6, 4, 3, 2, 0.5],
        })

        dose_df = pd.DataFrame({
            "time": [0],
            "dose": [100],
        })

        conc = NCAConcentration(
            data=conc_df,
            conc_col="conc",
            time_col="time",
        )

        dose = NCADose(
            data=dose_df,
            dose_col="dose",
            time_col="time",
        )

        data = NCAData(conc=conc, dose=dose)
        results = run_nca(data)

        # Check key results
        df = results.to_dataframe()

        # Check that parameters were calculated
        assert len(df) > 0

        # Filter to get cmax and tmax
        cmax_row = df[df["parameter"] == "cmax"]
        tmax_row = df[df["parameter"] == "tmax"]

        # Cmax should be 10
        assert cmax_row["value"].iloc[0] == 10
        # Tmax should be 1
        assert tmax_row["value"].iloc[0] == 1
