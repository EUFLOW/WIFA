code_saturne
============

.. image:: ../../img/Logo_code_saturne.png
    :align: left
    :width: 150

`code_saturne <https://www.code-saturne.org/cms/web/>`_ is a free open-source finite volume CFD solver for the Navier-Stokes equations, developped primarily by EDF. It solves the Navier-Stokes equations with scalar transport for 2D, 2D-axisymmetric, and 3D flows, whether steady or unsteady, laminar or turbulent, incompressible, dilatable, or weakly compressible, isothermal or not.

code_saturne contains modules dedicated to specific physics, namely for atmospheric flows. A detailed description of the modelling possibilities, including the atmospheric module, can be found in `code_saturne's online documentation <https://www.code-saturne.org/cms/web/documentation/v80/>`_. This page focuses on the modelling assumptions and added source files allowing to handle wind farm flow modelling in WIFA.

This API can be launched with version 8.0 of code_saturne and requires access to a high performance computer cluster.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   cs_equations
   cs_turbines
   cs_mesh
   cs_inflow

Configuration Options
---------------------

code_saturne-specific options are configured in the ``analysis.yaml`` file:

**Run Type:**

.. code-block:: yaml

   run_type: "simulate"    # Options: simulate, postprocess

**HPC Configuration:**

.. code-block:: yaml

   HPC_config:
     # Simulation run settings
     run_node_number: 1
     run_ntasks_per_node: 48
     run_wall_time_hours: 6
     run_partition: ""
     #
     # Meshing settings
     mesh_node_number: 1
     mesh_ntasks_per_node: 48
     mesh_wall_time_hours: 1
     mesh_partition: ""
     #
     wckey: ""              # Cluster accounting key

**Mesh Configuration:**

.. code-block:: yaml

   mesh:
     remesh: True           # Generate new mesh or use existing
     mesh_file_name: MESH/mesh.med

Workflow
~~~~~~~~

The code_saturne workflow consists of several automated steps:

1. **Initialization**: Parse WindIO input, create case directory structure
2. **Meshing**: Generate computational mesh using salome (if ``remesh: True``)
3. **Precursor** (optional): Run precursor simulation for inflow profiles
4. **Simulation**: Submit main wind farm simulation to HPC queue
5. **Postprocessing**: Extract turbine power and flow field data

Input Requirements
~~~~~~~~~~~~~~~~~~

code_saturne supports two types of inflow specification:

**Hub-height time series:**

- Wind speed and direction at hub height
- Stability via LMO (Monin-Obukhov length)
- Generates analytical profiles internally

**Vertical profiles:**

- Full atmospheric profiles (u, v, temperature, TKE, epsilon)
- Requires precursor simulation for turbulent inflow
- Supports capping inversion for stable stratification

Simulation Modes
~~~~~~~~~~~~~~~~

**Neutral boundary layer:**

- Uses log-law profiles
- No precursor required
- Energy equation disabled

**Stable stratification with capping inversion:**

- Generates Nieuwstadt profiles or uses provided profiles
- Precursor simulation for turbulent inflow
- Damping layer for gravity waves
- Energy equation enabled

API Reference
-------------

.. note::

   code_saturne simulations require access to an HPC cluster and pre-installed code_saturne (v8.0) and salome executables. See :doc:`/getting_started/installation` for setup instructions.

Usage
~~~~~

**Command Line:**

.. code-block:: console

   wifa_saturne path/to/system.yaml

**Python:**

.. code-block:: python

   from wifa.cs_api.cs_modules.csLaunch.cs_run_function import run_code_saturne

   # Run code_saturne simulation
   run_code_saturne("path/to/system.yaml")

   # Test mode (reduced iterations)
   run_code_saturne("system.yaml", test_mode=True)

   # Postprocess existing results
   run_code_saturne("system.yaml", postprocess_only=True)

.. py:function:: wifa.cs_api.cs_modules.csLaunch.cs_run_function.run_code_saturne(windio_input, test_mode=False, output_dir=None, postprocess_only=False)

   Run a code_saturne simulation from WindIO input.

   :param windio_input: Path to WindIO YAML file
   :type windio_input: str
   :param test_mode: Enable test mode (reduced iterations for validation)
   :type test_mode: bool
   :param output_dir: Custom output directory
   :type output_dir: str, optional
   :param postprocess_only: Only run postprocessing on existing results
   :type postprocess_only: bool

   **Example:**

   .. code-block:: python

      from wifa.cs_api.cs_modules.csLaunch.cs_run_function import run_code_saturne

      # Production run
      run_code_saturne("system.yaml")

      # Quick validation
      run_code_saturne("system.yaml", test_mode=True)

.. py:function:: wifa.cs_api.cs_modules.csLaunch.cs_run_function.initialize_cs_case_from_windio(windio_input, output_dir)

   Initialize a code_saturne case from WindIO input without running.

   :param windio_input: Path to WindIO YAML file
   :type windio_input: str
   :param output_dir: Output directory for the case
   :type output_dir: str
   :returns: Configured CS_study object
   :rtype: CS_study

Outputs
~~~~~~~

Results are stored in the case directory:

- ``Farm/RESU/<case_name>/power.txt``: Turbine power time series
- ``Farm/RESU/<case_name>/turbine_data.nc``: Per-turbine NetCDF output
- ``Farm/RESU/<case_name>/flow_field.nc``: Flow field data (if requested)
