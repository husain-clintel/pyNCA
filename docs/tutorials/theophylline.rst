Theophylline Example
====================

This tutorial demonstrates a complete NCA analysis using the theophylline dataset,
similar to the PKNCA vignette. Theophylline is a bronchodilator used to treat
respiratory diseases.

Dataset Description
-------------------

The theophylline dataset contains concentration-time data from 6 subjects
who received oral doses of theophylline.

.. code-block:: python

   import pandas as pd
   import pynca as nca

   # Load the theophylline dataset
   df = pd.read_csv("tests/data/theophylline.csv")
   print(df.head(12))

Output::

       Subject    Wt  Dose   Time   conc
    0        1  79.6  4.02   0.00   0.74
    1        1  79.6  4.02   0.25   2.84
    2        1  79.6  4.02   0.57   6.57
    3        1  79.6  4.02   1.12  10.50
    4        1  79.6  4.02   2.02   9.66
    5        1  79.6  4.02   3.82   8.58
    6        1  79.6  4.02   5.10   8.36
    7        1  79.6  4.02   7.03   7.47
    8        1  79.6  4.02   9.05   6.89
    9        1  79.6  4.02  12.12   5.94
   10        1  79.6  4.02  24.37   3.28
   11        2  72.4  4.40   0.00   0.00

Setting Up the Analysis
-----------------------

Step 1: Create Concentration Object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   conc = nca.NCAConcentration(
       data=df,
       conc_col='conc',
       time_col='Time',
       subject_col='Subject'
   )

   print(f"Number of subjects: {conc.n_subjects}")
   print(f"Subjects: {conc.subjects}")

Output::

   Number of subjects: 6
   Subjects: [1, 2, 3, 4, 5, 6]

Step 2: Create Dose Object
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Extract unique dose information
   dose_df = df[['Subject', 'Dose']].drop_duplicates()
   dose_df['Time'] = 0  # Dose at time 0

   dose = nca.NCADose(
       data=dose_df,
       dose_col='Dose',
       time_col='Time',
       subject_col='Subject'
   )

   print(f"Dose information:\n{dose.data}")

Step 3: Combine Data
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   data = nca.NCAData(conc=conc, dose=dose)
   print(data)

Output::

   NCAData(n_subjects=6, n_intervals=6, has_dose=True)

Running the Analysis
--------------------

.. code-block:: python

   results = nca.run_nca(data)
   print(f"Parameters calculated: {len(results)}")

Viewing Results
---------------

Long Format
^^^^^^^^^^^

.. code-block:: python

   df_results = results.to_dataframe()
   print(df_results.head(20))

Output::

       subject  interval_start  interval_end parameter      value
    0        1             0.0         24.37      cmax      10.50
    1        1             0.0         24.37      tmax       1.12
    2        1             0.0         24.37     tlast      24.37
    3        1             0.0         24.37     clast       3.28
    4        1             0.0         24.37  auc.last    148.92
    ...

Wide Format
^^^^^^^^^^^

.. code-block:: python

   df_wide = results.to_dataframe(wide=True)
   print(df_wide[['subject', 'cmax', 'tmax', 'auc.last', 'half.life']])

Output::

      subject   cmax  tmax   auc.last  half.life
   0        1  10.50  1.12    148.923      14.30
   1        2   8.33  1.92     91.527       6.66
   2        3   8.20  1.02     98.734       6.50
   3        4   8.60  1.07    106.796       6.98
   4        5  11.40  1.00    121.294       8.00
   5        6   7.09  3.52     90.701       7.79

Summary Statistics
------------------

.. code-block:: python

   summary = results.summary()
   print(summary)

This produces summary statistics (mean, SD, CV%, geometric mean) for each parameter.

Filtering Results
-----------------

Filter by Parameter
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Get only Cmax and AUClast
   filtered = results.filter(parameters=['cmax', 'auc.last'])
   print(filtered.to_dataframe())

Filter by Subject
^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Get results for specific subjects
   filtered = results.filter(subjects=[1, 2])
   print(filtered.to_dataframe(wide=True))

Get Specific Parameter Values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Get Cmax values as array
   cmax_values = results.get_parameter('cmax')
   print(f"Cmax values: {cmax_values}")
   print(f"Mean Cmax: {cmax_values.mean():.2f}")

Output::

   Cmax values: [10.5   8.33  8.2   8.6  11.4   7.09]
   Mean Cmax: 9.02

Exporting Results
-----------------

.. code-block:: python

   # Export to CSV
   results.to_csv("theophylline_results.csv")

   # Export to Excel (requires openpyxl)
   results.to_excel("theophylline_results.xlsx")

Individual Parameter Calculations
---------------------------------

You can also calculate individual parameters directly:

.. code-block:: python

   import numpy as np

   # Get data for Subject 1
   subj1 = df[df['Subject'] == 1]
   conc_array = subj1['conc'].values
   time_array = subj1['Time'].values

   # Calculate individual parameters
   cmax = nca.calc_cmax(conc_array, time_array)
   tmax = nca.calc_tmax(conc_array, time_array)
   auc_last = nca.calc_auc_last(conc_array, time_array)

   # Lambda_z returns a dictionary with regression details
   lambda_z_result = nca.calc_lambda_z(conc_array, time_array)

   print(f"Cmax: {cmax}")
   print(f"Tmax: {tmax}")
   print(f"AUClast: {auc_last:.2f}")
   print(f"Lambda_z: {lambda_z_result['lambda_z']:.4f}")
   print(f"R-squared: {lambda_z_result['r_squared']:.4f}")
   print(f"Points used: {lambda_z_result['n_points']}")

Output::

   Cmax: 10.5
   Tmax: 1.12
   AUClast: 148.92
   Lambda_z: 0.0485
   R-squared: 1.0000
   Points used: 3

Customizing the Analysis
------------------------

Using Different AUC Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from pynca import options

   # Use linear-up/log-down method
   options.auc_method = "linear-up/log-down"

   results_luld = nca.run_nca(data)

Specifying Calculation Intervals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Create data with custom interval
   data_custom = nca.NCAData(conc=conc, dose=dose)

   # Add specific interval (e.g., 0-12 hours)
   data_custom.add_interval(start=0, end=12)

   results_partial = nca.run_nca(data_custom)

Comparison with PKNCA
---------------------

pyNCA follows PKNCA naming conventions and calculation methods:

+----------------+----------------+
| PKNCA          | pyNCA          |
+================+================+
| PKNCAconc      | NCAConcentration|
+----------------+----------------+
| PKNCAdose      | NCADose        |
+----------------+----------------+
| PKNCAdata      | NCAData        |
+----------------+----------------+
| PKNCAresults   | NCAResults     |
+----------------+----------------+
| pk.nca()       | run_nca()      |
+----------------+----------------+

Parameter names also follow PKNCA conventions (e.g., ``auc.last``, ``lambda.z``, ``half.life``).
