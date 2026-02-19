PyWake
======

| `Documentation <https://topfarm.pages.windenergy.dtu.dk/PyWake/>`_
| `GitHub <https://github.com/DTUWindEnergy/PyWake>`_
| `Examples <https://topfarm.pages.windenergy.dtu.dk/PyWake/notebooks.html>`_

PyWake is an open-source wind farm simulation tool developed by DTU Wind Energy.
It provides fast and flexible wake modeling capabilities for wind farm flow and power estimation,
supporting both time series and probabilistic (Weibull) wind resource inputs.

Configuration Options
---------------------

PyWake-specific options are configured in the ``analysis.yaml`` file:

**Wake Deficit Model:**

.. code-block:: yaml

   wind_deficit_model:
     name: Bastankhah2014  # Options: Jensen, Bastankhah2014, SuperGaussian, TurbOPark, FUGA
     wake_expansion_coefficient:
       k_a: 0.0            # TI-dependent coefficient: k = k_a * TI + k_b
       k_b: 0.04           # Base wake expansion
       free_stream_ti: false
     ceps: 0.2             # Epsilon parameter (Bastankhah models)
     use_effective_ws: true

**Deflection Model:**

.. code-block:: yaml

   deflection_model:
     name: Jimenez         # Options: None, Jimenez
     beta: 0.1             # Jimenez deflection coefficient

**Turbulence Model:**

.. code-block:: yaml

   turbulence_model:
     name: CrespoHernandez  # Options: None, STF2005, STF2017, CrespoHernandez
     c1: 1.0                # STF model coefficients (if applicable)
     c2: 1.0

**Superposition Model:**

.. code-block:: yaml

   superposition_model:
     ws_superposition: Linear   # Options: Linear, Squared
     ti_superposition: Squared  # Options: Linear, Squared

**Rotor Averaging:**

.. code-block:: yaml

   rotor_averaging:
     name: center           # Options: center, avg_deficit
     grid: grid
     n_x_grid_points: 4
     n_y_grid_points: 4

**Blockage Model:**

.. code-block:: yaml

   blockage_model:
     name: None             # Options: None, SelfSimilarityDeficit2020, FUGA
     ss_alpha: 0.888        # SelfSimilarity parameter

API Reference
-------------

Usage
~~~~~

**Command Line:**

.. code-block:: console

   wifa_pywake path/to/system.yaml

**Python:**

.. code-block:: python

   from wifa.pywake_api import run_pywake

   # Run PyWake simulation
   aep = run_pywake("path/to/system.yaml")

   # Or pass a pre-loaded dictionary
   aep = run_pywake(system_dict, output_dir="results")

.. py:function:: wifa.pywake_api.run_pywake(yamlFile, output_dir="output")

   Run a PyWake simulation from WindIO input.

   :param yamlFile: Path to WindIO YAML file or pre-loaded dictionary
   :type yamlFile: str or dict
   :param output_dir: Output directory (overridden by YAML if specified)
   :type output_dir: str
   :returns: Annual Energy Production in GWh
   :rtype: float

   **Example:**

   .. code-block:: python

      from wifa.pywake_api import run_pywake

      aep = run_pywake("system.yaml")
      print(f"AEP: {aep:.2f} GWh")

Supported Input Formats
~~~~~~~~~~~~~~~~~~~~~~~

- **Time series**: NetCDF with ``time`` dimension for wind speed/direction
- **Weibull distribution**: NetCDF with sector probabilities and Weibull parameters
- **Turbine-specific data**: Per-turbine wind conditions with ``wind_turbine`` dimension

Outputs
~~~~~~~

- ``turbine_data.nc``: Power and effective wind speed per turbine
- ``FarmFlow.nc``: Flow field (wind speed, TI) on specified grid
- ``output.yaml``: Summary of simulation configuration


See also
~~~~~~~~
See also :doc:`pywake_ellipsys` for the PyWakeEllipSys CFD coupling.
