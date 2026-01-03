Getting Started
===============

This guide will walk you through the basics of using pyNCA for non-compartmental
pharmacokinetic analysis.

Installation
------------

Install pyNCA using pip:

.. code-block:: bash

   pip install pynca

For development:

.. code-block:: bash

   pip install pynca[dev]

Dependencies
------------

pyNCA requires:

* Python >= 3.9
* numpy >= 1.21
* pandas >= 1.3
* scipy >= 1.7
* pint >= 0.19

Basic Concepts
--------------

pyNCA is built around four core data structures:

1. **NCAConcentration**: Stores concentration-time data
2. **NCADose**: Stores dosing information
3. **NCAData**: Combines concentration and dose data with analysis intervals
4. **NCAResults**: Stores and exports NCA results

Basic Workflow
--------------

Step 1: Prepare Your Data
^^^^^^^^^^^^^^^^^^^^^^^^^

Your data should be in a pandas DataFrame format:

.. code-block:: python

   import pandas as pd
   import pynca as nca

   # Concentration data
   conc_df = pd.DataFrame({
       'subject': [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
       'time': [0, 1, 2, 4, 8, 0, 1, 2, 4, 8],
       'conc': [0, 10, 8, 6, 3, 0, 12, 9, 7, 4]
   })

   # Dose data
   dose_df = pd.DataFrame({
       'subject': [1, 2],
       'time': [0, 0],
       'dose': [100, 100]
   })

Step 2: Create NCA Objects
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Create concentration object
   conc = nca.NCAConcentration(
       data=conc_df,
       conc_col='conc',
       time_col='time',
       subject_col='subject'
   )

   # Create dose object
   dose = nca.NCADose(
       data=dose_df,
       dose_col='dose',
       time_col='time',
       subject_col='subject'
   )

   # Combine into NCAData
   data = nca.NCAData(conc=conc, dose=dose)

Step 3: Run NCA Analysis
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Run analysis
   results = nca.run_nca(data)

   # View results as DataFrame
   df = results.to_dataframe()
   print(df)

   # Wide format (one row per subject)
   df_wide = results.to_dataframe(wide=True)
   print(df_wide)

Step 4: Export Results
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Export to CSV
   results.to_csv("nca_results.csv")

   # Export to Excel
   results.to_excel("nca_results.xlsx")

   # Get summary statistics
   summary = results.summary()
   print(summary)

Calculated Parameters
---------------------

pyNCA calculates the following parameters by default:

**Exposure Parameters**

* ``cmax`` - Maximum concentration
* ``tmax`` - Time of maximum concentration
* ``clast`` - Last measurable concentration
* ``tlast`` - Time of last measurable concentration
* ``auc.last`` - AUC to last measurable concentration
* ``auc.inf.obs`` - AUC extrapolated to infinity (observed)
* ``auc.pct.extrap.obs`` - Percent of AUC extrapolated
* ``aumc.last`` - AUMC to last measurable concentration
* ``aumc.inf.obs`` - AUMC extrapolated to infinity

**Terminal Phase Parameters**

* ``lambda.z`` - Terminal elimination rate constant
* ``half.life`` - Terminal half-life

**Clearance and Volume**

* ``cl.obs`` - Clearance (observed)
* ``vz.obs`` - Volume of distribution (terminal phase)
* ``vss.obs`` - Volume of distribution at steady state
* ``mrt.last`` - Mean residence time

Configuration Options
---------------------

pyNCA provides global options for customizing calculations:

.. code-block:: python

   from pynca import options

   # Change AUC method
   options.auc_method = "linear-up/log-down"

   # Change minimum points for lambda_z
   options.lambda_z_min_points = 4

   # Change BLQ handling
   options.blq_handling = "drop"

Available AUC Methods
^^^^^^^^^^^^^^^^^^^^^

* ``"linear"`` - Linear trapezoidal rule (default)
* ``"log"`` - Log-linear trapezoidal rule
* ``"linear-up/log-down"`` - Linear for ascending, log for descending

Available BLQ Handling Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``"zero"`` - Replace with zero
* ``"drop"`` - Remove BLQ values
* ``"loq/2"`` - Replace with LOQ/2
* ``"position"`` - Position-aware handling (leading BLQ to zero, trailing BLQ dropped)

Next Steps
----------

* See the :doc:`tutorials/theophylline` tutorial for a complete example
* Explore the :doc:`api/core` reference for detailed API documentation
