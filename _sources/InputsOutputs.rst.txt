Inputs and Outputs - WindIO
---------------------------

Inputs
~~~~~~~~~~~~~~~~~~~~~~~
The inputs to WIFA can be specified using the windIO schema. Few example can be found in the `WIFA examples folder <https://github.com/EUFLOW/WIFA/-/tree/main/examples/cases?ref_type=heads>`_.

Tool choice should be specified in the "wind_energy_system" yaml file. `This one <https://github.com/EUFLOW/WIFA/-/blob/main/examples/cases/windio_4turbines/wind_energy_system/system.yaml?ref_type=heads>`_ for example when running the  `windio_4turbines example case <https://github.com/EUFLOW/WIFA/-/tree/main/examples/cases/windio_4turbines?ref_type=heads>`_

Depending on the tool you choose, you may specify additional computation options. This can be specified through "analysis" yaml. `This one <https://github.com/EUFLOW/WIFA/-/blob/main/examples/cases/windio_4turbines/wind_energy_system/analysis.yaml?ref_type=heads>`_ for example when running the  `windio_4turbines example case <https://github.com/EUFLOW/WIFA/-/tree/main/examples/cases/windio_4turbines?ref_type=heads>`_. If any options of the other tools are kept in this file, they are simply ignored.

An example of this is the wake expansion model parameters for the engineering tools Pywake and FOXES:

.. code-block:: console

  #pywake and foxes
  wind_deficit_model:
    name: Bastankhah2014
    wake_expansion_coefficient: # k = ka*ti + kb
      k_a: 0.0
      k_b: 0.04
      free_stream_ti: false
    ceps: 0.2
    use_effective_ws: true

And optional slurm parameters when running the CFD tool code_saturne:

.. code-block:: console

  #code_saturne
  HPC_config:
    run_node_number: 1
    run_ntasks_per_node: 48
    run_wall_time_hours: 6
    run_partition: ""
    #
    mesh_node_number: 1
    mesh_ntasks_per_node: 48
    mesh_wall_time_hours: 1
    run_partition: ""
    #
    wckey: ""

Farm layout, turbine specifications, can be described as explained in the windIO documentation.

The wind resource for which the power should be computed can be described in a netcdf file, and specified in the "energy_resource" yaml file. `Here <https://github.com/EUFLOW/WIFA/-/tree/main/examples/cases/windio_4turbines/plant_energy_resource?ref_type=heads>`_ for example when running the `windio_4turbines example case <https://github.com/EUFLOW/WIFA/-/tree/main/examples/cases/windio_4turbines?ref_type=heads>`_. The minimal requirement in this netcdf is to have a "wind_velocity" and "wind_direction" variables, with dimension "time", to give a hub height information. The "time" dimension does not necessarily need to be physical time, but can also be a list of realizations. Other variables (LMO for stability, z0, capping_inversion_strength, etc.) can also be specified.

To specify vertical profiles, one has to add a dimension "height" to the "wind_speed" and "wind_direction" variables. The "potential_temperature" can also be specified as a variable in this latter case to account for thermal stratification, which is of particular interest when using the APM and CFD tools.

Outputs
~~~~~~~~~~~~~~~~~~~~~~~
The output of the WIFA is a results folder containing two netcdf files:
 - turbine data : with default power and rotor_effective_velocity for each turbine and at each time
 - flow data : with default wind_speed and wind_direction at x and y coordinates, and at z=hub height.

The options of the api will evolve to allow additionnal outputs.
