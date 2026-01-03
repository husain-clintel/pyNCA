pyNCA Documentation
===================

**pyNCA** is a Python package for Pharmacokinetic Non-Compartmental Analysis (NCA),
designed to replicate the functionality of the `PKNCA R package <https://github.com/billdenney/pknca>`_.

.. note::

   This project is under active development.

Features
--------

* **Comprehensive NCA calculations**: AUC, AUMC, Cmax, Tmax, half-life, clearance, and more
* **Multiple AUC methods**: Linear trapezoidal, log-linear, linear-up/log-down
* **Automatic lambda_z selection**: Optimal terminal phase point selection
* **BLQ handling**: Multiple strategies for below-limit-of-quantification data
* **Sparse NCA**: Bailer's method for destructive sampling designs
* **Steady-state tools**: Superposition and accumulation calculations
* **Unit support**: Physical unit handling with pint

Installation
------------

.. code-block:: bash

   pip install pynca

Quick Start
-----------

.. code-block:: python

   import pynca as nca
   import pandas as pd

   # Load concentration data
   conc_df = pd.read_csv("concentrations.csv")
   dose_df = pd.read_csv("doses.csv")

   # Create NCA objects
   conc = nca.NCAConcentration(
       data=conc_df,
       conc_col="concentration",
       time_col="time",
       subject_col="subject"
   )

   dose = nca.NCADose(
       data=dose_df,
       dose_col="dose",
       time_col="time",
       subject_col="subject"
   )

   # Run analysis
   data = nca.NCAData(conc=conc, dose=dose)
   results = nca.run_nca(data)

   # View results
   print(results.to_dataframe(wide=True))

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   getting_started
   tutorials/theophylline

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/core
   api/calc
   api/analysis

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
