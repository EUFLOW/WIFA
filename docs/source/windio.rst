WindIO
======

WIFA uses the WindIO schema to standardize inputs across all flow models. This ensures that the same input files work with PyWake, foxes, WAYVE, and code_saturne.

Input Structure
---------------

A WIFA simulation requires a hierarchical set of YAML files:

.. code-block:: text

   wind_energy_system/
   ├── system.yaml          # Main entry point
   └── analysis.yaml        # Model-specific configuration

   plant_energy_site/
   └── site.yaml            # Site boundaries and resource reference

   plant_energy_resource/
   └── resource.nc          # Wind resource data (NetCDF)

   plant_wind_farm/
   ├── farm.yaml            # Turbine layout
   └── turbine.yaml         # Turbine specifications

Wind Energy System (system.yaml)
--------------------------------

The main entry point that references all other components:

.. code-block:: yaml

   name: My Wind Farm Study
   site: !include ../plant_energy_site/site.yaml
   wind_farm: !include ../plant_wind_farm/farm.yaml
   attributes:
     flow_model:
       name: foxes           # Required: pywake, foxes, wayve, codesaturne
     analysis: !include analysis.yaml
     model_outputs_specification:
       output_folder: "results"
       run_configuration:
         times_run:
           all_occurences: false
           subset: [0, 2, 5]  # Specific time indices to run
       turbine_outputs:
         turbine_nc_filename: 'turbine_data.nc'
         output_variables: ['power', 'rotor_effective_velocity']
       flow_field:
         report: true
         flow_nc_filename: flow_field.nc
         output_variables: ['wind_speed', 'wind_direction']
         z_planes:
           z_sampling: "hub_heights"  # or "plane_list"
           z_list: [100, 119, 150]    # if plane_list
           xy_sampling: "grid"
           x_bounds: [-1000, 5000]
           y_bounds: [-1000, 1000]
           dx: 150
           dy: 150

Site Configuration (site.yaml)
------------------------------

Defines site boundaries and references the wind resource:

.. code-block:: yaml

   name: My Site
   boundaries:
     polygons:
       - x: [0, 5000, 5000, 0]      # Polygon vertices (meters)
         y: [0, 0, 2000, 2000]
   energy_resource: !include ../plant_energy_resource/resource.yaml

Wind Farm Layout (farm.yaml)
----------------------------

**Single Turbine Type:**

.. code-block:: yaml

   name: My Wind Farm
   layouts:
     - coordinates:
         x: [0, 1000, 2000, 3000]   # Turbine x positions (meters)
         y: [0, 0, 0, 0]            # Turbine y positions (meters)
   turbines: !include turbine.yaml

**Multiple Turbine Types:**

.. code-block:: yaml

   name: My Wind Farm
   layouts:
     - coordinates:
         x: [0, 1000, 2000, 3000]
         y: [0, 0, 0, 0]
       turbine_types: [1, 2, 2, 3]  # Index into turbine_types dict
   turbine_types:
     1: !include turbine_type1.yaml
     2: !include turbine_type2.yaml
     3: !include turbine_type3.yaml

Turbine Definition (turbine.yaml)
---------------------------------

.. code-block:: yaml

   name: IEA 15MW Reference
   hub_height: 150.0              # meters
   rotor_diameter: 240.0          # meters
   performance:
     # Option A: Power curve
     power_curve:
       power_values: [0, 500000, ..., 15000000]    # Watts
       power_wind_speeds: [3, 4, ..., 25]          # m/s
     # Option B: Cp curve (alternative to power_curve)
     Cp_curve:
       Cp_values: [0.0, 0.2, 0.45, ...]            # dimensionless
       Cp_wind_speeds: [3, 4, 5, ...]              # m/s
     # Required with either option:
     Ct_curve:
       Ct_values: [0.0, 0.8, 0.75, ...]            # dimensionless
       Ct_wind_speeds: [3, 4, 5, ...]              # m/s

Wind Resource Formats
---------------------

Wind resource data is can be provided as NetCDF files linked with the ``!include`` directive. WIFA supports several formats:

Time Series
~~~~~~~~~~~

NetCDF structure for time-varying wind conditions:

.. code-block:: text

   Dimensions:
     time: N                       # Number of time steps

   Coordinates:
     time: float64 (time,)         # Time index or timestamp

   Data Variables:
     wind_speed: float64 (time,)           # Hub-height wind speed [m/s]
     wind_direction: float64 (time,)       # Wind direction [degrees, 0=N, 90=E]
     turbulence_intensity: float64 (time,) # TI as fraction [0-1]

**Creating time series in Python:**

.. code-block:: python

   import xarray as xr
   import numpy as np

   n_times = 1000
   ds = xr.Dataset(
       coords={"time": np.arange(n_times, dtype=np.float64)},
       data_vars={
           "wind_speed": ("time", np.random.weibull(2, n_times) * 8 + 4),
           "wind_direction": ("time", np.random.uniform(0, 360, n_times)),
           "turbulence_intensity": ("time", np.full(n_times, 0.06)),
       },
   )
   ds.to_netcdf("timeseries_resource.nc")

Turbine-Specific Time Series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For per-turbine wind conditions (e.g., SCADA data):

.. code-block:: text

   Dimensions:
     time: N                       # Number of time steps
     wind_turbine: M               # Number of turbines

   Coordinates:
     time: float64 (time,)
     wind_turbine: int64 (wind_turbine,)   # Turbine indices [0, 1, 2, ...]

   Data Variables:
     wind_speed: float64 (time, wind_turbine)
     wind_direction: float64 (time, wind_turbine)
     turbulence_intensity: float64 (time, wind_turbine)
     operating: float64 (time, wind_turbine)  # 1.0=on, 0.0=off

Weibull Wind Rose
~~~~~~~~~~~~~~~~~

For probabilistic AEP calculations:

.. code-block:: text

   Dimensions:
     wind_direction: D             # Number of direction sectors

   Coordinates:
     wind_direction: float64 (wind_direction,)  # Sector centers [degrees]

   Data Variables:
     sector_probability: float64 (wind_direction,)  # Sum = 1.0
     weibull_a: float64 (wind_direction,)           # Scale parameter [m/s]
     weibull_k: float64 (wind_direction,)           # Shape parameter [-]
     turbulence_intensity: float64 (wind_direction,)

Heterogeneous Wind Rose (per turbine)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For spatially varying wind resources:

.. code-block:: text

   Dimensions:
     wind_turbine: M               # Number of turbines
     wind_direction: D             # Number of direction sectors

   Coordinates:
     wind_turbine: int64 (wind_turbine,)
     wind_direction: float64 (wind_direction,)

   Data Variables:
     sector_probability: float64 (wind_turbine, wind_direction)
     weibull_a: float64 (wind_turbine, wind_direction)
     weibull_k: float64 (wind_turbine, wind_direction)
     turbulence_intensity: float64 (wind_turbine, wind_direction)
     x: float64 (wind_turbine,)    # Turbine x position [m]
     y: float64 (wind_turbine,)    # Turbine y position [m]
     height: float64 (wind_turbine,)  # Reference height [m]

Variable Conventions
--------------------

+---------------------------+-------------+---------------+----------------------------------+
| Variable                  | Units       | Valid Range   | Notes                            |
+===========================+=============+===============+==================================+
| ``wind_speed``            | m/s         | >= 0          | At hub height                    |
+---------------------------+-------------+---------------+----------------------------------+
| ``wind_direction``        | degrees     | [0, 360)      | Meteorological: 0=N, 90=E        |
+---------------------------+-------------+---------------+----------------------------------+
| ``turbulence_intensity``  | fraction    | [0, 1]        | 0.06 = 6% TI                     |
+---------------------------+-------------+---------------+----------------------------------+
| ``sector_probability``    | fraction    | [0, 1]        | Sum = 1.0                        |
+---------------------------+-------------+---------------+----------------------------------+
| ``weibull_a``             | m/s         | > 0           | Scale parameter                  |
+---------------------------+-------------+---------------+----------------------------------+
| ``weibull_k``             | dimensionless| > 0          | Shape (2.0 = Rayleigh)           |
+---------------------------+-------------+---------------+----------------------------------+
| ``operating``             | flag        | 0 or 1        | 1 = on, 0 = off                  |
+---------------------------+-------------+---------------+----------------------------------+

Outputs
-------

WIFA outputs results as NetCDF files in the specified output folder:

Turbine Data (turbine_data.nc)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Per-turbine results for each simulation state:

.. code-block:: text

   Dimensions:
     turbine: M                    # Number of turbines
     time/wd/ws: N                 # Simulation states

   Data Variables:
     power: float64 (states, turbine)              # Power [W]
     rotor_effective_velocity: float64 (states, turbine)  # [m/s]

Flow Field (flow_field.nc)
~~~~~~~~~~~~~~~~~~~~~~~~~~

Spatial flow data on the requested grid:

.. code-block:: text

   Dimensions:
     x: Nx
     y: Ny
     z: Nz
     time/wd/ws: N

   Data Variables:
     wind_speed: float64 (states, x, y, z)         # [m/s]
     wind_direction: float64 (states, x, y, z)     # [degrees]
     turbulence_intensity: float64 (states, x, y, z)  # [fraction]

See the `windIO documentation <https://windio.readthedocs.io/>`_ for the full schema specification.
