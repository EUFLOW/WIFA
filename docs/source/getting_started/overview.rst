Overview
========

What is WIFA?
-------------

WIFA (Wind Farm API) is a multi-fidelity wind farm simulation framework that integrates different flow modeling tools through a unified Python interface. It allows researchers and engineers to run wind farm simulations using various models—from fast engineering wake models to high-fidelity CFD—all with the same input format.

**Key Features:**

- Unified interface for multiple flow models
- Standardized input/output using WindIO schema
- Support for time series and probabilistic wind resources
- Flexible model configuration through YAML files
- CLI tools for easy execution

Architecture Overview
---------------------

WIFA integrates four flow modeling tools, each suited for different use cases:

.. image:: ../../img/wifa_diagram.png
   :align: center
   :width: 80%

**Engineering Wake Models (Fast)**

- **PyWake**: DTU's wind farm modeling framework, optimized for rapid AEP calculations
- **foxes**: Fraunhofer IWES's vectorized wind farm modeling framework, designed for large-scale scenarios

**Atmospheric Perturbation Model (Medium Fidelity)**

- **wayve**: KU Leuven's atmospheric perturbation model capturing wind farm-atmosphere interactions and gravity wave effects

**CFD (High Fidelity)**

- **code_saturne**: EDF's open-source CFD solver with atmospheric boundary layer capabilities

Choosing the Right Tool
-----------------------

+----------------+------------------+-------------------+----------------------------+
| Tool           | Speed            | Physics           | Best For                   |
+================+==================+===================+============================+
| PyWake         | Very fast        | Wake deficits     | AEP estimation, layout opt |
+----------------+------------------+-------------------+----------------------------+
| foxes          | Very fast        | Wake deficits     | Large farms, time series   |
+----------------+------------------+-------------------+----------------------------+
| wayve          | Medium           | Gravity waves,    | Wind farm blockage,        |
|                |                  | atmospheric       | atmospheric coupling       |
|                |                  | coupling          |                            |
+----------------+------------------+-------------------+----------------------------+
| code_saturne   | Slow (HPC)       | Full RANS CFD     | Detailed flow fields,      |
|                |                  |                   | complex terrain            |
+----------------+------------------+-------------------+----------------------------+

Quickstart Example
------------------

**1. Prepare Input Files**

WIFA uses WindIO-formatted YAML files to describe:

- Wind energy system configuration (``system.yaml``)
- Site and energy resource (NetCDF wind data)
- Wind farm layout and turbine specifications
- Analysis options and model parameters

Example ``system.yaml``:

.. code-block:: yaml

   name: My Wind Farm Study
   site: !include ../site/site.yaml
   wind_farm: !include ../wind_farm/farm.yaml
   attributes:
     flow_model:
       name: foxes  # or pywake, wayve, codesaturne
     analysis: !include analysis.yaml
     model_outputs_specification:
       output_folder: "results"
       turbine_outputs:
         turbine_nc_filename: 'turbine_data.nc'
         output_variables: ['power', 'rotor_effective_velocity']

**2. Run a Simulation**

Using the CLI:

.. code-block:: console

   wifa path/to/system.yaml

Using Python:

.. code-block:: python

   from wifa.main_api import run_api

   # Run simulation
   run_api("path/to/system.yaml")

**3. Access Results**

Results are saved to the specified output folder as NetCDF files:

.. code-block:: python

   import xarray as xr

   # Load turbine data
   turbine_data = xr.open_dataset("results/turbine_data.nc")
   print(turbine_data.power)

   # Load flow field (if requested)
   flow_field = xr.open_dataset("results/flow_field.nc")
   print(flow_field.wind_speed)

Main API
--------

The ``flow_model.name`` field in your YAML selects the tool:

- ``pywake`` → ``wifa.pywake_api.run_pywake()``
- ``foxes`` → ``wifa.foxes_api.run_foxes()``
- ``wayve`` → ``wifa.wayve_api.run_wayve()``
- ``codesaturne`` → ``wifa.cs_api.run_code_saturne()``

.. py:function:: wifa.main_api.run_api(yaml_input)

   Main entry point for WIFA simulations.

   :param yaml_input: Path to WindIO-formatted YAML file
   :type yaml_input: str
   :raises ValueError: If flow_model.name is invalid

   **Example:**

   .. code-block:: python

      from wifa.main_api import run_api

      run_api("examples/cases/windio_4turbines/wind_energy_system/system.yaml")

.. py:function:: wifa.main_api.run()

   CLI entry point. Parses command line arguments and calls ``run_api()``.

   This function is registered as the ``wifa`` console script.

CLI Tools
---------

WIFA provides five CLI commands, installed automatically with the package:

.. code-block:: console

   wifa <system.yaml>            # Auto-selects tool from flow_model.name
   wifa_pywake <system.yaml>     # Run directly with PyWake
   wifa_foxes <system.yaml>      # Run directly with foxes
   wifa_wayve <system.yaml>      # Run directly with wayve
   wifa_saturne <system.yaml>    # Run directly with code_saturne

All commands return exit code ``0`` on success and ``1`` on error.
See each tool's page for detailed options and examples.

Next Steps
----------

- See :doc:`installation` for detailed installation instructions
- Explore the :doc:`examples` for complete worked cases
- Read :doc:`/windio` for WindIO schema details
