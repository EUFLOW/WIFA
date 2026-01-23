Examples
========

All examples are in ``examples/cases/``. Run any example with:

.. code-block:: console

   wifa examples/cases/<case_name>/wind_energy_system/system.yaml

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Case
     - Description
   * - ``windio_4turbines``
     - 4-turbine farm, time series input, flow field output
   * - ``simple_wind_rose``
     - Weibull wind rose with sector probabilities
   * - ``heterogeneous_wind_rose_at_turbines``
     - Per-turbine heterogeneous wind rose
   * - ``heterogeneous_wind_rose_map``
     - Gridded heterogeneous wind rose with spatial interpolation
   * - ``turbine_specific_speeds_timeseries``
     - Per-turbine wind speeds and directions (e.g. SCADA)
   * - ``timeseries_with_operating_flag``
     - Time series with turbine on/off operating flags
   * - ``open_source_scada``
     - Real SCADA data format
   * - ``windio_4turbines_multipleTurbines``
     - Multiple turbine types in one farm
   * - ``windio_4turbines_ABL``
     - Atmospheric boundary layer with capping inversion
   * - ``windio_4turbines_ABL_stable``
     - Stable atmospheric stratification
   * - ``windio_4turbines_profiles_stable``
     - Full vertical profile input for stable conditions
   * - ``AWAKEN``
     - AWAKEN field campaign layout
   * - ``KUL_LES``
     - KU Leuven LES comparison case

Example Structure
-----------------

Each example follows this directory structure:

.. code-block:: text

   examples/cases/<case_name>/
   ├── wind_energy_system/
   │   ├── system.yaml          # Main entry point
   │   └── analysis.yaml        # Model configuration
   ├── plant_energy_site/
   │   └── site.yaml            # Site boundaries, resource reference
   ├── plant_energy_resource/
   │   └── resource.nc          # Wind resource NetCDF
   └── plant_wind_farm/
       ├── farm.yaml            # Farm layout
       └── turbine.yaml         # Turbine specification

Creating Your Own Case
----------------------

1. Copy an existing example closest to your use case
2. Modify the wind resource NetCDF with your data
3. Update turbine coordinates in ``farm.yaml``
4. Adjust turbine specifications in ``turbine.yaml``
5. Configure model parameters in ``analysis.yaml``
6. Run with your preferred tool

See :doc:`/windio` for detailed schema documentation.
