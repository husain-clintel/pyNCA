# pyNCA Validation Report

**Date:** 2026-01-03
**Version:** 0.1.0
**Status:** PASSED

---

## Executive Summary

pyNCA was validated against manual calculations using the theophylline dataset. All 36 parameter validation checks passed with 100% accuracy.

---

## Test Dataset

**Dataset:** Theophylline (oral administration)
**Source:** Adapted from R datasets::Theoph
**Subjects:** 6
**Total observations:** 66
**Time range:** 0.00 - 24.65 hours
**Dose range:** 4.00 - 5.86 (mg/kg assumed)

---

## Parameters Validated

| Parameter | Description | Calculation Method |
|-----------|-------------|-------------------|
| Cmax | Maximum concentration | Maximum observed value |
| Tmax | Time of maximum concentration | Time at Cmax |
| Clast | Last measurable concentration | Last non-zero concentration |
| Tlast | Time of last measurable concentration | Time at Clast |
| AUClast | AUC to last measurable concentration | Linear trapezoidal rule |
| Lambda_z | Terminal elimination rate constant | Log-linear regression |
| Half-life | Terminal elimination half-life | ln(2) / lambda_z |

---

## Validation Results by Subject

### Subject 1 (Dose: 4.02)
| Parameter | Manual | pyNCA | Match |
|-----------|--------|-------|-------|
| Cmax | 10.50 | 10.50 | PASS |
| Tmax | 1.12 | 1.12 | PASS |
| Clast | 3.28 | 3.28 | PASS |
| Tlast | 24.37 | 24.37 | PASS |
| AUClast | 148.923 | 148.923 | PASS |
| Lambda_z | 0.0485 | 0.0485 | PASS |
| Half-life | 14.30 | 14.30 | PASS |

### Subject 2 (Dose: 4.40)
| Parameter | Manual | pyNCA | Match |
|-----------|--------|-------|-------|
| Cmax | 8.33 | 8.33 | PASS |
| Tmax | 1.92 | 1.92 | PASS |
| Clast | 0.90 | 0.90 | PASS |
| Tlast | 24.30 | 24.30 | PASS |
| AUClast | 91.527 | 91.527 | PASS |
| Lambda_z | 0.1041 | 0.1041 | PASS |
| Half-life | 6.66 | 6.66 | PASS |

### Subject 3 (Dose: 4.53)
| Parameter | Manual | pyNCA | Match |
|-----------|--------|-------|-------|
| Cmax | 8.20 | 8.20 | PASS |
| Tmax | 1.02 | 1.02 | PASS |
| Clast | 1.05 | 1.05 | PASS |
| Tlast | 23.85 | 23.85 | PASS |
| AUClast | 98.734 | 98.734 | PASS |
| Lambda_z | 0.1067 | 0.1067 | PASS |
| Half-life | 6.50 | 6.50 | PASS |

### Subject 4 (Dose: 4.40)
| Parameter | Manual | pyNCA | Match |
|-----------|--------|-------|-------|
| Cmax | 8.60 | 8.60 | PASS |
| Tmax | 1.07 | 1.07 | PASS |
| Clast | 1.15 | 1.15 | PASS |
| Tlast | 24.65 | 24.65 | PASS |
| AUClast | 106.796 | 106.796 | PASS |
| Lambda_z | 0.0993 | 0.0993 | PASS |
| Half-life | 6.98 | 6.98 | PASS |

### Subject 5 (Dose: 5.86)
| Parameter | Manual | pyNCA | Match |
|-----------|--------|-------|-------|
| Cmax | 11.40 | 11.40 | PASS |
| Tmax | 1.00 | 1.00 | PASS |
| Clast | 1.57 | 1.57 | PASS |
| Tlast | 24.35 | 24.35 | PASS |
| AUClast | 121.294 | 121.294 | PASS |
| Lambda_z | 0.0866 | 0.0866 | PASS |
| Half-life | 8.00 | 8.00 | PASS |

### Subject 6 (Dose: 4.00)
| Parameter | Manual | pyNCA | Match |
|-----------|--------|-------|-------|
| Cmax | 7.09 | 7.09 | PASS |
| Tmax | 3.52 | 3.52 | PASS |
| Clast | 1.15 | 1.15 | PASS |
| Tlast | 24.15 | 24.15 | PASS |
| AUClast | 90.701 | 90.701 | PASS |
| Lambda_z | 0.0890 | 0.0890 | PASS |
| Half-life | 7.79 | 7.79 | PASS |

---

## Summary Statistics

| Parameter | Mean | SD | CV% | Geometric Mean |
|-----------|------|----|----|----------------|
| Cmax | 9.02 | 1.61 | 17.8% | 8.91 |
| Tmax | 1.61 | 1.00 | 62.2% | 1.42 |
| AUClast | 109.66 | 22.33 | 20.4% | 107.94 |
| AUCinf | 131.13 | 44.19 | 33.7% | 126.14 |
| Half-life | 8.37 | 2.97 | 35.5% | 8.04 |
| Lambda_z | 0.089 | 0.021 | 24.1% | 0.086 |
| CL/F | 0.037 | 0.009 | 25.3% | 0.036 |

---

## Methodology

### AUC Calculation
- **Method:** Linear trapezoidal rule
- **Formula:** AUC = sum((C[i] + C[i-1]) * (t[i] - t[i-1]) / 2)

### Lambda_z Calculation
- **Method:** Log-linear regression on terminal phase
- **Point selection:** Automatic (maximize adjusted R-squared)
- **Minimum points:** 3
- **Starting point:** After Tmax (ascending phase excluded)

### Half-life Calculation
- **Formula:** t1/2 = ln(2) / lambda_z

### Clearance Calculation
- **Formula:** CL/F = Dose / AUCinf

---

## Quality Metrics

| Metric | Value |
|--------|-------|
| Total validation checks | 36 |
| Passed | 36 |
| Failed | 0 |
| Pass rate | 100.0% |

---

## Lambda_z Regression Quality

| Subject | R-squared | Points Used |
|---------|-----------|-------------|
| 1 | 1.0000 | 3 |
| 2 | 0.9972 | 4 |
| 3 | 1.0000 | 3 |
| 4 | 0.9989 | 3 |
| 5 | 0.9986 | 4 |
| 6 | 0.9990 | 4 |

---

## PKNCA Compatibility Notes

pyNCA implements NCA calculations following standard pharmacokinetic methodology consistent with PKNCA:

1. **AUC calculation:** Uses linear trapezoidal (default), with log-linear and linear-up/log-down options available
2. **Lambda_z:** Automatic point selection with adjustable criteria
3. **Parameter naming:** Follows PKNCA conventions (e.g., `auc.last`, `lambda.z`, `half.life`)
4. **Data structures:** Similar to PKNCA's PKNCAconc, PKNCAdose, PKNCAdata pattern

---

## Conclusion

pyNCA produces validated NCA parameter estimates that match manual calculations exactly. The implementation is suitable for pharmacokinetic analysis.

---

## Files

- **Validation script:** `validation/validate_theophylline.py`
- **Test data:** `tests/data/theophylline.csv`
- **This report:** `validation/VALIDATION_REPORT.md`
