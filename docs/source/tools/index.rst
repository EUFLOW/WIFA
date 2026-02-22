Tools
=====

WIFA integrates four flow modeling tools, each suited for different use cases and fidelity levels.

.. list-table::
   :header-rows: 1
   :widths: 15 20 65

   * - Tool
     - Speed
     - Use Case
   * - :doc:`PyWake <pywake>`
     - Fast
     - AEP estimation, layout optimization
   * - :doc:`foxes <foxes>`
     - Fast
     - Large farms, long time series
   * - :doc:`wayve <wayve>`
     - Medium
     - Gravity waves, farm blockage
   * - :doc:`code_saturne <codesaturne>`
     - Slow (HPC)
     - Detailed flow analysis

.. toctree::
   :hidden:

   pywake
   foxes
   wayve
   codesaturne
   pywake_ellipsys
