API
===

All functions accept a WindIO YAML file path (or pre-loaded dict).

.. code-block:: python

   from wifa import run_api, run_foxes, run_pywake, run_wayve, run_code_saturne

wifa.main_api
-------------

.. py:function:: run_api(yaml_input)

   Run simulation using the tool specified by ``flow_model.name`` in the YAML.

wifa.foxes_api
--------------

.. py:function:: run_foxes(input_yaml, input_dir=None, output_dir=None, engine="default", n_procs=None, chunksize_states=None, chunksize_points=None, verbosity=1, **kwargs)

   Run a foxes simulation. Returns ``(farm_results, point_results, outputs)``.

wifa.pywake_api
---------------

.. py:function:: run_pywake(yamlFile, output_dir="output")

   Run a PyWake simulation. Returns AEP in GWh.

wifa.wayve_api
--------------

.. py:function:: run_wayve(yamlFile, output_dir="output", debug_mode=False)

   Run a WAYVE simulation. Outputs written to NetCDF files.

wifa.cs_api
-----------

.. py:function:: run_code_saturne(windio_input, test_mode=False, output_dir=None, postprocess_only=False)

   Run a code_saturne simulation. Requires HPC with code_saturne v8.0 installed.

.. py:function:: initialize_cs_case_from_windio(windio_input, output_dir)

   Initialize a code_saturne case without running. Returns ``CS_study`` object.
