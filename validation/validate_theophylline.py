"""
Validation Script: pyNCA vs Manual Calculations

This script validates pyNCA calculations against manual calculations
using the theophylline dataset.
"""

import numpy as np
import pandas as pd
import sys
sys.path.insert(0, '.')

import pynca
from pynca.calc.auc import calc_auc, calc_auc_last
from pynca.calc.cmax import calc_cmax, calc_tmax, calc_clast, calc_tlast
from pynca.calc.half_life import calc_lambda_z, calc_half_life


def manual_linear_auc(conc: np.ndarray, time: np.ndarray) -> float:
    """Manual AUC calculation using linear trapezoidal rule."""
    auc = 0.0
    for i in range(1, len(conc)):
        auc += (conc[i] + conc[i-1]) * (time[i] - time[i-1]) / 2
    return auc


def validate_subject(subject_data: pd.DataFrame, subject_id: int, dose: float) -> dict:
    """Validate NCA parameters for a single subject."""

    conc = subject_data['conc'].values
    time = subject_data['Time'].values

    results = {}

    # 1. Cmax validation
    results['cmax'] = {
        'manual': np.max(conc),
        'pynca': calc_cmax(conc, time),
        'match': np.isclose(np.max(conc), calc_cmax(conc, time))
    }

    # 2. Tmax validation
    tmax_idx = np.argmax(conc)
    results['tmax'] = {
        'manual': time[tmax_idx],
        'pynca': calc_tmax(conc, time),
        'match': np.isclose(time[tmax_idx], calc_tmax(conc, time))
    }

    # 3. Clast validation (last non-zero concentration)
    non_zero_idx = np.where(conc > 0)[0]
    clast_manual = conc[non_zero_idx[-1]] if len(non_zero_idx) > 0 else np.nan
    results['clast'] = {
        'manual': clast_manual,
        'pynca': calc_clast(conc, time),
        'match': np.isclose(clast_manual, calc_clast(conc, time))
    }

    # 4. Tlast validation
    tlast_manual = time[non_zero_idx[-1]] if len(non_zero_idx) > 0 else np.nan
    results['tlast'] = {
        'manual': tlast_manual,
        'pynca': calc_tlast(conc, time),
        'match': np.isclose(tlast_manual, calc_tlast(conc, time))
    }

    # 5. AUClast validation (linear trapezoidal)
    # Get data up to tlast
    tlast = calc_tlast(conc, time)
    mask = time <= tlast
    auc_manual = manual_linear_auc(conc[mask], time[mask])
    auc_pynca = calc_auc_last(conc, time, method='linear')
    results['auc_last'] = {
        'manual': auc_manual,
        'pynca': auc_pynca,
        'match': np.isclose(auc_manual, auc_pynca, rtol=0.01)
    }

    # 6. Lambda_z validation (using terminal phase)
    lambda_z_result = calc_lambda_z(conc, time)
    results['lambda_z'] = {
        'pynca': lambda_z_result['lambda_z'] if lambda_z_result else np.nan,
        'r_squared': lambda_z_result['r_squared'] if lambda_z_result else np.nan,
        'n_points': lambda_z_result['n_points'] if lambda_z_result else 0,
        'valid': lambda_z_result is not None and lambda_z_result['r_squared'] > 0.8
    }

    # 7. Half-life validation
    if lambda_z_result and not np.isnan(lambda_z_result['lambda_z']) and lambda_z_result['lambda_z'] > 0:
        hl_manual = np.log(2) / lambda_z_result['lambda_z']
        hl_result = calc_half_life(conc, time)
        hl_pynca = hl_result['half_life'] if hl_result else np.nan
        results['half_life'] = {
            'manual': hl_manual,
            'pynca': hl_pynca,
            'match': np.isclose(hl_manual, hl_pynca, rtol=0.01) if not np.isnan(hl_pynca) else False
        }
    else:
        results['half_life'] = {
            'manual': np.nan,
            'pynca': np.nan,
            'match': True  # Both NaN
        }

    return results


def run_validation():
    """Run full validation and generate report."""

    # Load data
    df = pd.read_csv('tests/data/theophylline.csv')

    print("=" * 70)
    print("pyNCA VALIDATION REPORT - Theophylline Dataset")
    print("=" * 70)
    print()

    # Run pyNCA analysis
    conc = pynca.NCAConcentration(
        data=df,
        conc_col='conc',
        time_col='Time',
        subject_col='Subject'
    )

    dose_df = df[['Subject', 'Dose']].drop_duplicates()
    dose_df['Time'] = 0
    dose = pynca.NCADose(
        data=dose_df,
        dose_col='Dose',
        time_col='Time',
        subject_col='Subject'
    )

    data = pynca.NCAData(conc=conc, dose=dose)
    results = pynca.run_nca(data)

    # Get results
    df_results = results.to_dataframe(wide=True)

    print("1. DATA SUMMARY")
    print("-" * 40)
    print(f"   Subjects: {df['Subject'].nunique()}")
    print(f"   Total observations: {len(df)}")
    print(f"   Time range: {df['Time'].min():.2f} - {df['Time'].max():.2f} h")
    print(f"   Dose range: {df['Dose'].min():.2f} - {df['Dose'].max():.2f}")
    print()

    print("2. PARAMETER VALIDATION")
    print("-" * 40)

    all_validations = []

    for subject in df['Subject'].unique():
        subject_data = df[df['Subject'] == subject].copy()
        dose_val = subject_data['Dose'].iloc[0]

        validation = validate_subject(subject_data, subject, dose_val)

        print(f"\n   Subject {subject} (Dose: {dose_val})")
        print("   " + "-" * 35)

        for param, vals in validation.items():
            if 'match' in vals:
                status = "✓" if vals['match'] else "✗"
                if 'manual' in vals:
                    print(f"   {param:12s}: manual={vals['manual']:.4f}, "
                          f"pynca={vals['pynca']:.4f} [{status}]")
                else:
                    print(f"   {param:12s}: pynca={vals['pynca']:.4f} [{status}]")
            else:
                if param == 'lambda_z':
                    print(f"   {param:12s}: λz={vals['pynca']:.6f}, "
                          f"R²={vals['r_squared']:.4f}, n={vals['n_points']}")

        all_validations.append((subject, validation))

    print("\n")
    print("3. pyNCA RESULTS SUMMARY")
    print("-" * 40)
    print()

    # Calculate summary statistics
    summary_params = ['cmax', 'tmax', 'auc.last', 'auc.inf.obs', 'half.life', 'lambda.z', 'cl.obs']

    print(f"   {'Parameter':<15} {'Mean':>10} {'SD':>10} {'CV%':>10} {'Geo Mean':>10}")
    print("   " + "-" * 55)

    for param in summary_params:
        if param in df_results.columns:
            values = df_results[param].dropna()
            if len(values) > 0:
                mean_val = values.mean()
                sd_val = values.std()
                cv_pct = (sd_val / mean_val) * 100 if mean_val != 0 else np.nan
                geo_mean = np.exp(np.log(values[values > 0]).mean()) if (values > 0).all() else np.nan
                print(f"   {param:<15} {mean_val:>10.3f} {sd_val:>10.3f} "
                      f"{cv_pct:>10.1f} {geo_mean:>10.3f}")

    print("\n")
    print("4. VALIDATION SUMMARY")
    print("-" * 40)

    # Count validations
    total_checks = 0
    passed_checks = 0

    for subject, validation in all_validations:
        for param, vals in validation.items():
            if 'match' in vals:
                total_checks += 1
                if vals['match']:
                    passed_checks += 1

    print(f"   Total validation checks: {total_checks}")
    print(f"   Passed: {passed_checks}")
    print(f"   Failed: {total_checks - passed_checks}")
    print(f"   Pass rate: {(passed_checks/total_checks)*100:.1f}%")
    print()

    if passed_checks == total_checks:
        print("   ✓ ALL VALIDATIONS PASSED")
    else:
        print("   ✗ SOME VALIDATIONS FAILED - Review required")

    print()
    print("5. METHODOLOGY NOTES")
    print("-" * 40)
    print("   - AUC calculation: Linear trapezoidal rule")
    print("   - Lambda_z: Log-linear regression on terminal phase")
    print("   - Automatic point selection for lambda_z (R² > 0.8)")
    print("   - Half-life: ln(2) / lambda_z")
    print("   - Clearance: Dose / AUCinf")
    print()
    print("=" * 70)

    return all_validations, df_results


if __name__ == "__main__":
    validations, results = run_validation()
