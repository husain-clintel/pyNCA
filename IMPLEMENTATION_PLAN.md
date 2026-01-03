# pyNCA Implementation Plan

A Python package that replicates PKNCA R package functionality for Pharmacokinetic Non-Compartmental Analysis.

---

## Phase 1: Project Setup & Package Structure

### 1.1 Initialize Package Structure
```
pyNCA/
├── pyproject.toml
├── setup.py
├── README.md
├── LICENSE (AGPL-3 to match PKNCA)
├── pynca/
│   ├── __init__.py
│   ├── core/
│   │   ├── __init__.py
│   │   ├── concentration.py    # PKNCAconc equivalent
│   │   ├── dose.py             # PKNCAdose equivalent
│   │   ├── data.py             # PKNCAdata equivalent
│   │   ├── results.py          # PKNCAresults equivalent
│   │   └── options.py          # PKNCA.options equivalent
│   ├── calc/
│   │   ├── __init__.py
│   │   ├── auc.py              # AUC calculations
│   │   ├── aumc.py             # AUMC calculations
│   │   ├── cmax.py             # Cmax, Cmin, Tmax
│   │   ├── half_life.py        # Half-life calculations
│   │   ├── clearance.py        # CL, Vss, Vz calculations
│   │   ├── bioavailability.py  # F calculations
│   │   └── parameters.py       # Parameter registry
│   ├── cleaning/
│   │   ├── __init__.py
│   │   ├── blq.py              # BLQ handling
│   │   ├── imputation.py       # Data imputation
│   │   └── exclusion.py        # Data exclusion rules
│   ├── interpolation/
│   │   ├── __init__.py
│   │   ├── interpolate.py      # Concentration interpolation
│   │   └── extrapolate.py      # Concentration extrapolation
│   ├── analysis/
│   │   ├── __init__.py
│   │   ├── nca.py              # pk.nca equivalent
│   │   ├── intervals.py        # Interval selection
│   │   ├── sparse.py           # Sparse NCA calculations
│   │   ├── steady_state.py     # Time to steady state
│   │   └── superposition.py    # Superposition calculations
│   ├── summary/
│   │   ├── __init__.py
│   │   ├── summarize.py        # Result summarization
│   │   └── units.py            # Unit handling/conversion
│   └── utils/
│       ├── __init__.py
│       ├── validation.py       # Input validation
│       └── helpers.py          # Utility functions
├── tests/
│   ├── __init__.py
│   ├── test_auc.py
│   ├── test_half_life.py
│   ├── test_clearance.py
│   ├── test_integration.py
│   └── data/                   # Test datasets
│       └── theophylline.csv
└── docs/
    ├── conf.py
    ├── index.rst
    └── tutorials/
```

### 1.2 Dependencies
```toml
[project]
dependencies = [
    "numpy>=1.21",
    "pandas>=1.3",
    "scipy>=1.7",
    "pint>=0.19",           # Unit handling
]

[project.optional-dependencies]
dev = ["pytest", "pytest-cov", "black", "ruff", "mypy"]
docs = ["sphinx", "sphinx-rtd-theme"]
```

---

## Phase 2: Core Data Structures

### 2.1 NCAConcentration Class (`core/concentration.py`)
- Store concentration-time data with subject/group identifiers
- Support formula-like syntax for grouping: `conc ~ time | subject`
- Validate concentration and time data
- Handle multiple analytes
- Methods: `filter()`, `get_subject()`, `to_dataframe()`

### 2.2 NCADose Class (`core/dose.py`)
- Store dosing information (time, amount, route, duration)
- Support multiple dose types: bolus, infusion, extravascular
- Link doses to subjects
- Handle single and multiple dosing scenarios

### 2.3 NCAData Class (`core/data.py`)
- Combine concentration and dose data
- Automatically determine calculation intervals
- Support custom interval specification
- Methods: `set_intervals()`, `get_intervals()`, `merge()`

### 2.4 NCAResults Class (`core/results.py`)
- Store calculated NCA parameters
- Support filtering and subsetting
- Export to DataFrame, CSV, Excel
- Methods: `summary()`, `to_dataframe()`, `filter()`, `exclude()`

### 2.5 Options Class (`core/options.py`)
- Global configuration for calculations
- Business rules customization
- Key options:
  - `auc.method`: "linear", "log-linear", "linear-up/log-down"
  - `lambda_z.n_points`: min points for half-life
  - `lambda_z.method`: regression method
  - `blq.handling`: "zero", "drop", "loq/2"
  - `summary.stats`: default statistics
  - `time.units`, `conc.units`, `dose.units`

---

## Phase 3: Calculation Functions

### 3.1 AUC Calculations (`calc/auc.py`)
```python
def calc_auc(conc, time, method="linear", interval=None):
    """Area Under the Curve calculation"""

def calc_auc_last(conc, time, method="linear"):
    """AUC from first to last measurable concentration"""

def calc_auc_inf(conc, time, lambda_z, clast):
    """AUC extrapolated to infinity"""

def calc_auc_int(conc, time, start, end, method="linear"):
    """AUC over specified interval with interpolation"""

def calc_auc_pct_extrap(auc_last, auc_inf):
    """Percent of AUC extrapolated"""
```

**Integration Methods:**
- Linear trapezoidal
- Log-linear trapezoidal
- Linear-up/log-down (mixed)

### 3.2 AUMC Calculations (`calc/aumc.py`)
```python
def calc_aumc(conc, time, method="linear"):
    """Area Under the Moment Curve"""

def calc_aumc_last(conc, time, method="linear"):
    """AUMC to last measurable concentration"""

def calc_aumc_inf(conc, time, lambda_z, clast, tlast):
    """AUMC extrapolated to infinity"""
```

### 3.3 Peak Parameters (`calc/cmax.py`)
```python
def calc_cmax(conc, time, check_tmax=True):
    """Maximum concentration"""

def calc_cmin(conc, time):
    """Minimum concentration"""

def calc_tmax(conc, time, first=True):
    """Time of maximum concentration"""

def calc_cav(auc, tau):
    """Average concentration over interval"""

def calc_ctrough(conc, time):
    """Trough concentration"""

def calc_swing(cmax, cmin):
    """Swing (Cmax - Cmin) / Cmin"""

def calc_ptf(cmax, cmin, cav):
    """Peak-trough fluctuation"""
```

### 3.4 Half-Life Calculations (`calc/half_life.py`)
```python
def calc_half_life(conc, time, tlast=None, min_points=3,
                   max_points=None, method="ols"):
    """
    Terminal elimination half-life
    Returns: half_life, lambda_z, r_squared, adj_r_squared,
             lambda_z_intercept, n_points, time_range
    """

def select_lambda_z_points(conc, time, method="best_fit"):
    """Automatic selection of terminal phase points"""

def calc_half_life_tobit(conc, time, loq):
    """Half-life using Tobit regression for censored data"""
```

**Point Selection Methods:**
- Best fit (maximize adjusted R²)
- Fixed number of points
- Manual specification

### 3.5 Clearance & Volume (`calc/clearance.py`)
```python
def calc_cl(dose, auc_inf, f=1.0):
    """Total clearance (CL or CL/F)"""

def calc_cl_last(dose, auc_last, f=1.0):
    """Clearance based on AUClast"""

def calc_vz(cl, lambda_z):
    """Volume of distribution (Vz)"""

def calc_vss(dose, aumc_inf, auc_inf, f=1.0):
    """Volume of distribution at steady state"""

def calc_mrt(aumc_inf, auc_inf):
    """Mean residence time"""

def calc_mrt_iv(aumc_inf, auc_inf, duration=0):
    """MRT corrected for IV infusion"""
```

### 3.6 Bioavailability (`calc/bioavailability.py`)
```python
def calc_f(auc_test, auc_ref, dose_test, dose_ref):
    """Absolute or relative bioavailability"""

def calc_f_with_cl(cl_test, cl_ref):
    """Bioavailability from clearance ratio"""
```

### 3.7 Parameter Registry (`calc/parameters.py`)
```python
# Define all available parameters with metadata
PARAMETERS = {
    "cmax": {
        "function": calc_cmax,
        "description": "Maximum concentration",
        "unit_type": "concentration",
        "requires": ["conc", "time"],
    },
    "auc.last": {
        "function": calc_auc_last,
        "description": "AUC to last measurable concentration",
        "unit_type": "concentration*time",
        "requires": ["conc", "time"],
    },
    # ... all parameters
}

def get_parameter_info(param_name):
    """Get parameter metadata"""

def list_parameters():
    """List all available parameters"""
```

---

## Phase 4: Data Cleaning & Handling

### 4.1 BLQ Handling (`cleaning/blq.py`)
```python
def clean_blq(conc, time, method="zero", loq=None):
    """
    Handle Below Limit of Quantification values
    Methods:
    - "zero": Replace with 0
    - "drop": Remove BLQ values
    - "loq/2": Replace with LOQ/2
    - "before_first": Drop BLQ before first measurable
    - "after_last": Keep/drop BLQ after Clast
    """
```

### 4.2 Missing Data (`cleaning/imputation.py`)
```python
def clean_na(conc, time, method="drop"):
    """Handle missing concentration values"""

def impute_conc(conc, time, target_time, method="linear"):
    """Impute concentration at specific time point"""
```

### 4.3 Data Exclusion (`cleaning/exclusion.py`)
```python
def exclude_points(data, criteria):
    """Exclude data points based on criteria"""

def flag_outliers(conc, time, method="mad"):
    """Flag potential outliers"""
```

---

## Phase 5: Interpolation & Extrapolation

### 5.1 Interpolation (`interpolation/interpolate.py`)
```python
def interpolate_conc(conc, time, target_time, method="linear"):
    """
    Interpolate concentration at target time
    Methods: linear, log-linear, dose-aware
    """

def interpolate_conc_dose_aware(conc, time, target_time,
                                 dose_time, dose_amount):
    """Dose-aware concentration interpolation"""
```

### 5.2 Extrapolation (`interpolation/extrapolate.py`)
```python
def extrapolate_conc(conc, time, target_time, lambda_z, clast, tlast):
    """Extrapolate concentration beyond last measurement"""
```

---

## Phase 6: Main Analysis Engine

### 6.1 NCA Analysis (`analysis/nca.py`)
```python
class NCA:
    """Main NCA analysis class"""

    def __init__(self, data: NCAData, options: Options = None):
        self.data = data
        self.options = options or get_default_options()

    def calculate(self, parameters: list = None) -> NCAResults:
        """
        Run NCA calculations for specified parameters
        If parameters=None, calculate all applicable parameters
        """

    def calculate_interval(self, start, end, parameters):
        """Calculate parameters for specific interval"""

    def _calculate_single_profile(self, conc, time, dose, interval):
        """Calculate all parameters for single concentration profile"""
```

### 6.2 Interval Selection (`analysis/intervals.py`)
```python
def auto_select_intervals(conc_data, dose_data):
    """Automatically determine calculation intervals based on dosing"""

def create_intervals(start, end, parameters=None):
    """Create interval specification"""

class Interval:
    """Interval specification class"""
    start: float
    end: float
    parameters: list
    name: str
```

### 6.3 Sparse NCA (`analysis/sparse.py`)
```python
def sparse_auc(conc, time, subject, method="batch"):
    """
    AUC calculation for sparse/destructive sampling
    Methods: Bailer, batch means, bootstrap
    """

def sparse_auc_se(conc, time, subject):
    """Standard error for sparse AUC"""
```

### 6.4 Steady State (`analysis/steady_state.py`)
```python
def time_to_steady_state(conc, time, dose_times, method="monoexponential"):
    """Estimate time to reach steady state"""

def check_steady_state(conc, time, trough_times, method="tost"):
    """Test if steady state has been achieved"""
```

### 6.5 Superposition (`analysis/superposition.py`)
```python
def superposition(conc, time, tau, n_doses=None, method="concentration"):
    """
    Predict concentrations at steady state from single dose data
    using superposition principle
    """
```

---

## Phase 7: Summary & Reporting

### 7.1 Summarization (`summary/summarize.py`)
```python
def summarize_results(results: NCAResults,
                      stats=["n", "mean", "sd", "cv", "median", "min", "max"],
                      by=None):
    """
    Summarize NCA results across subjects
    """

def geometric_mean(values):
    """Calculate geometric mean"""

def geometric_cv(values):
    """Calculate geometric coefficient of variation"""
```

### 7.2 Unit Handling (`summary/units.py`)
```python
def assign_units(results, conc_unit, time_unit, dose_unit):
    """Assign units to calculated parameters"""

def convert_units(value, from_unit, to_unit):
    """Convert between units using pint"""

def derive_parameter_units(param_name, base_units):
    """Derive appropriate units for each parameter"""
```

---

## Phase 8: Testing & Validation

### 8.1 Unit Tests
- Test each calculation function against known values
- Test edge cases (empty data, single point, all zeros)
- Compare results with PKNCA R package output

### 8.2 Integration Tests
- Full workflow tests with theophylline dataset
- Multi-subject analysis
- Multiple dosing scenarios

### 8.3 Validation Dataset
Create validation comparing pyNCA vs PKNCA results:
```python
# tests/test_validation.py
def test_vs_pknca_theophylline():
    """Validate results match PKNCA R package"""
```

---

## Phase 9: Documentation

### 9.1 API Documentation
- Docstrings for all public functions
- Sphinx autodoc generation

### 9.2 Tutorials
1. Quick Start Guide
2. Theophylline Example (matching PKNCA vignette)
3. Custom Intervals
4. Handling BLQ Data
5. Sparse Sampling Analysis
6. Unit Handling

### 9.3 Migration Guide
- Mapping PKNCA functions to pyNCA equivalents

---

## NCA Parameters - Complete List

| Parameter | Description | Function |
|-----------|-------------|----------|
| `cmax` | Maximum concentration | `calc_cmax` |
| `cmin` | Minimum concentration | `calc_cmin` |
| `tmax` | Time of Cmax | `calc_tmax` |
| `tlast` | Time of last measurable conc | `calc_tlast` |
| `clast` | Last measurable concentration | `calc_clast` |
| `auc.last` | AUC to Clast | `calc_auc_last` |
| `auc.inf.obs` | AUC extrapolated (observed Clast) | `calc_auc_inf` |
| `auc.inf.pred` | AUC extrapolated (predicted Clast) | `calc_auc_inf` |
| `auc.all` | AUC including time 0 | `calc_auc_all` |
| `auc.pct.extrap` | % AUC extrapolated | `calc_auc_pct_extrap` |
| `aumc.last` | AUMC to Clast | `calc_aumc_last` |
| `aumc.inf.obs` | AUMC extrapolated | `calc_aumc_inf` |
| `lambda.z` | Terminal rate constant | `calc_half_life` |
| `half.life` | Terminal half-life | `calc_half_life` |
| `r.squared` | R² for lambda.z regression | `calc_half_life` |
| `adj.r.squared` | Adjusted R² | `calc_half_life` |
| `lambda.z.n.points` | Points used for lambda.z | `calc_half_life` |
| `cl.obs` | Clearance (observed) | `calc_cl` |
| `cl.pred` | Clearance (predicted) | `calc_cl` |
| `vz.obs` | Volume of distribution | `calc_vz` |
| `vss.obs` | Vss (observed) | `calc_vss` |
| `mrt.last` | MRT to Clast | `calc_mrt` |
| `mrt.inf.obs` | MRT extrapolated | `calc_mrt` |
| `vss.iv` | Vss for IV dosing | `calc_vss` |
| `cav` | Average concentration | `calc_cav` |
| `ctrough` | Trough concentration | `calc_ctrough` |
| `swing` | Swing | `calc_swing` |
| `ptf` | Peak-trough fluctuation | `calc_ptf` |
| `accumulation.index` | Accumulation index | `calc_accumulation` |
| `f` | Bioavailability | `calc_f` |

---

## Implementation Order

1. **Week 1**: Project setup, core data structures (Phase 1-2)
2. **Week 2**: Basic calculations - AUC, Cmax, Tmax (Phase 3.1-3.3)
3. **Week 3**: Half-life, clearance, volume (Phase 3.4-3.6)
4. **Week 4**: Data cleaning, interpolation (Phase 4-5)
5. **Week 5**: Main analysis engine (Phase 6.1-6.2)
6. **Week 6**: Advanced features - sparse, superposition (Phase 6.3-6.5)
7. **Week 7**: Summary, units, reporting (Phase 7)
8. **Week 8**: Testing, validation, documentation (Phase 8-9)

---

## Usage Example (Target API)

```python
import pynca as nca
import pandas as pd

# Load data
conc_df = pd.read_csv("concentrations.csv")
dose_df = pd.read_csv("doses.csv")

# Create concentration object
conc = nca.NCAConcentration(
    data=conc_df,
    conc_col="concentration",
    time_col="time",
    subject_col="subject"
)

# Create dose object
dose = nca.NCADose(
    data=dose_df,
    dose_col="dose",
    time_col="time",
    subject_col="subject"
)

# Combine into analysis dataset
data = nca.NCAData(conc=conc, dose=dose)

# Set options
nca.options.set(
    auc_method="linear-up/log-down",
    blq_handling="zero"
)

# Run analysis
results = nca.run_nca(data)

# View results
print(results.summary())

# Export
results.to_csv("nca_results.csv")
```

---

## References

- [PKNCA CRAN](https://cran.r-project.org/web/packages/PKNCA/index.html)
- [PKNCA GitHub](https://github.com/billdenney/pknca)
- [PKNCA Documentation](https://humanpred.github.io/pknca/)
