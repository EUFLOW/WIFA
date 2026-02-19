CLI
===

All commands accept a WindIO YAML file as input. Exit code ``0`` on success, ``1`` on error.

.. code-block:: text

   wifa <input.yaml>
   wifa_pywake <input.yaml>
   wifa_foxes <input.yaml> [options]
   wifa_wayve <input.yaml>
   wifa_saturne <input.yaml>

wifa_foxes options
------------------

.. code-block:: text

   -o, --output_dir         Output directory
   -e, --engine             Computation engine (default, process, thread, dask)
   -n, --n_procs            Number of processes
   -c, --chunksize_states   Chunk size for states dimension
   -C, --chunksize_points   Chunk size for points dimension (default: 5000)
   -it, --iterative         Use iterative algorithm
   -v, --verbosity          Verbosity level, 0 = silent (default: 1)
