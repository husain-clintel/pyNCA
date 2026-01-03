# pyNCA

[![Python Version](https://img.shields.io/pypi/pyversions/pynca)](https://pypi.org/project/pynca/)
[![PyPI Version](https://img.shields.io/pypi/v/pynca)](https://pypi.org/project/pynca/)
[![License](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Tests](https://github.com/pyNCA/pynca/actions/workflows/ci.yml/badge.svg)](https://github.com/pyNCA/pynca/actions/workflows/ci.yml)
[![Documentation](https://readthedocs.org/projects/pynca/badge/?version=latest)](https://pynca.readthedocs.io/en/latest/?badge=latest)

A Python package for Pharmacokinetic Non-Compartmental Analysis (NCA).

pyNCA is a Python implementation inspired by the [PKNCA R package](https://github.com/billdenney/pknca), designed to perform all standard NCA calculations for pharmacokinetic data.

## Features

### Core NCA Calculations
- **Peak Parameters**: Cmax, Tmax, Cmin, Clast, Tlast, Cav, Ctrough, Swing, PTF
- **AUC Methods**: Linear, log-linear, linear-up/log-down trapezoidal integration
- **AUC Variants**: AUClast, AUCinf, AUCall, AUCint, partial AUC, dose-normalized AUC
- **AUMC Calculations**: AUMC, AUMClast, AUMCinf
- **Half-Life**: Automatic terminal phase selection, effective half-life, accumulation half-life
- **Clearance/Volume**: CL, Vz, Vss, MRT
- **Bioavailability**: F calculation, accumulation index

### Advanced Features
- **Multi-Analyte Support**: Parent drug and metabolite analysis with ratio calculations
- **Population Summaries**: Group comparisons, bioequivalence analysis, dose proportionality
- **Sparse Sampling**: NCA for destructive/sparse sampling designs (Bailer method)
- **Bootstrap CI**: Percentile and BCa confidence intervals
- **Parallel Processing**: Multi-core support for large datasets
- **Visualization**: Concentration-time plots, lambda_z diagnostics, forest plots

### Data Handling
- **BLQ Handling**: Multiple strategies (zero, drop, LOQ/2, position-aware)
- **Missing Data**: Imputation methods
- **Unit Handling**: Automatic unit assignment and conversion with pint
- **Pandas Integration**: Seamless DataFrame workflows

## Installation

```bash
# Basic installation
pip install pynca

# With plotting support
pip install pynca[plot]

# With all optional dependencies
pip install pynca[all]

# Development installation
pip install pynca[dev]
```

### Docker

```bash
# Build the image
docker build -t pynca .

# Run with docker-compose
docker-compose up jupyter  # Start Jupyter Lab on port 8888
```

## Quick Start

```python
import pynca as nca
import pandas as pd

# Load your data
conc_df = pd.DataFrame({
    'subject': [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
    'time': [0, 1, 2, 4, 8, 0, 1, 2, 4, 8],
    'concentration': [0, 8.5, 6.2, 3.1, 0.8, 0, 9.1, 6.8, 3.5, 0.9],
})

dose_df = pd.DataFrame({
    'subject': [1, 2],
    'time': [0, 0],
    'dose': [100, 100],
})

# Create NCA objects
conc = nca.NCAConcentration(
    data=conc_df,
    conc_col='concentration',
    time_col='time',
    subject_col='subject'
)

dose = nca.NCADose(
    data=dose_df,
    dose_col='dose',
    time_col='time',
    subject_col='subject'
)

# Combine and analyze
data = nca.NCAData(conc=conc, dose=dose)
results = nca.run_nca(data)

# View results
print(results.to_dataframe())
```

## Documentation

Full documentation is available at [https://pynca.readthedocs.io](https://pynca.readthedocs.io)

## Requirements

**Python**: 3.9+

**Required Dependencies**:
- numpy >= 1.21
- pandas >= 1.3
- scipy >= 1.7
- pint >= 0.19

**Optional Dependencies**:
- matplotlib >= 3.5 (for plotting)
- openpyxl >= 3.0 (for Excel export)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Install development dependencies (`pip install -e .[dev]`)
4. Make your changes and add tests
5. Run tests (`pytest tests/`)
6. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
7. Push to the branch (`git push origin feature/AmazingFeature`)
8. Open a Pull Request

## License

AGPL-3.0 (matching PKNCA R package)

## Citation

If you use pyNCA in your research, please cite:

```
pyNCA: Python Non-Compartmental Analysis Package
https://github.com/pyNCA/pynca
```

## Acknowledgments

This package is inspired by and aims to replicate the functionality of the excellent [PKNCA R package](https://github.com/billdenney/pknca) by Bill Denney.
