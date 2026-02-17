wayve
=====

*WAYVE* is an open-source code Python framework for the atmospheric perturbation model (APM).

The APM is an atmospheric perturbation model developed to simulate the interaction between wind farms and the atmosphere. It focuses on the generation and feedback effects of gravity waves. This code implements the model described in the following papers:

- Allaerts D., Meyers J. **Sensitivity and feedback of wind-farm-induced gravity waves**. Journal of Fluid Mechanics, 2019. `DOI <https://doi.org/10.1017/jfm.2018.969>`_
- Devesse K., Lanzilao L., & Meyers J. **A meso-micro atmospheric perturbation model for wind farm blockage**. Preprint, 2023. `DOI <https://doi.org/10.48550/arXiv.2310.18748>`_
- Devesse K., Stipa S., Brinkerhoff J., Allaerts D., & Meyers J. **Comparing methods for coupling wake models to an atmospheric perturbation model in WAYVE**. Journal of Physics: Conference Series, 2024. `DOI <https://iopscience.iop.org/article/10.1088/1742-6596/2767/9/092079/meta>`_

Configuration Options
---------------------

WAYVE-specific options are configured in the ``analysis.yaml`` file:

**APM Grid Settings:**

.. code-block:: yaml

   apm_grid:
     Lx: 1.e6              # Grid size in x-direction [m]
     Ly: 1.e6              # Grid size in y-direction [m]
     dx: 500               # Grid spacing [m]
     L_filter: 1.e3        # Spatial filter length [m]

**Layer Description:**

.. code-block:: yaml

   layers_description:
     farm_layer_height: 238.       # Lower layer height [m]
     number_of_fa_layers: 1        # Number of free atmosphere layers (1 or more)

**Wake Model Coupling:**

.. code-block:: yaml

   wm_coupling:
     method: "PB"          # Options: PB (Pressure-Based), VM (Velocity Matching), US (Upstream)
     settings:             # Method-specific settings
       alpha: 1.0          # VM method parameter
       distance: 500       # US method parameter [m]
     subgrid:
       include_subgrid: false
       D_to_dx: 0.5

**Additional APM Terms:**

.. code-block:: yaml

   APM_additional_terms:
     momentum_entrainment:
       mfp_type: "constant_flux"   # Options: None, constant_flux
       apm_mfp_settings:
         a_mfp: 0.120              # Momentum flux coefficient
         d_mfp: 27.80              # Distance parameter
     apm_disp_stresses:
       ds_type: subgrid            # Dispersive stresses type

**Wake Model Settings:**

.. code-block:: yaml

   wake_tool: wayve        # Options: wayve, foxes
   wind_deficit_model:
     name: Bastankhah2014
     wake_expansion_coefficient:
       k_a: 0.0
       k_b: 0.04
     ceps: 0.2

Input Requirements
~~~~~~~~~~~~~~~~~~

WAYVE requires atmospheric profile data for proper ABL initialization:

**Hub-height time series (minimal):**

- Wind speed and direction at hub height
- Optional: LMO (Monin-Obukhov length), friction velocity, roughness height

**Vertical profiles (recommended for stratified flows):**

- Wind speed and direction profiles
- Potential temperature profile
- Turbulence intensity profile
- Optional: TKE and dissipation profiles

**Capping inversion data (for gravity waves):**

- ABL height
- Inversion strength (dtheta)
- Inversion depth (dH)
- Free atmosphere lapse rate

API Reference
-------------

Usage
~~~~~

**Command Line:**

.. code-block:: console

   wifa_wayve path/to/system.yaml

**Python:**

.. code-block:: python

   from wifa.wayve_api import run_wayve

   # Run WAYVE simulation
   run_wayve("path/to/system.yaml", output_dir="results")

.. py:function:: wifa.wayve_api.run_wayve(yamlFile, output_dir="output", debug_mode=False)

   Run a WAYVE (APM) simulation from WindIO input.

   :param yamlFile: Path to WindIO YAML file or pre-loaded dictionary
   :type yamlFile: str or dict
   :param output_dir: Output directory for results
   :type output_dir: str
   :param debug_mode: Enable debug mode (wake model only, no APM solve)
   :type debug_mode: bool

   **Example:**

   .. code-block:: python

      from wifa.wayve_api import run_wayve

      run_wayve("system.yaml", output_dir="wayve_results")

Outputs
~~~~~~~

- ``turbine_data.nc``: Power and rotor effective velocity per turbine
- ``flow_field.nc``: Wind speed and direction on specified grid (if requested)
