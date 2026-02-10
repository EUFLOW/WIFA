.. WIFA documentation master file

=====================================
WIFA Documentation
=====================================

**WIFA** (Wind Farm API) is a multi-fidelity wind farm simulation framework that integrates multiple flow modeling tools through a unified Python interface.

.. image:: ../img/wifa_diagram.png
   :align: center
   :width: 80%
   :alt: WIFA Architecture Diagram

----

Supported Tools
---------------

.. list-table::
   :header-rows: 1
   :widths: 15 20 15 50

   * - Tool
     - Type
     - Speed
     - Use Case
   * - **PyWake**
     - Engineering wake model
     - Fast
     - AEP estimation, layout optimization
   * - **foxes**
     - Engineering wake model
     - Fast
     - Large farms, long time series
   * - **wayve**
     - Atmospheric perturbation
     - Medium
     - Gravity waves, farm blockage
   * - **code_saturne**
     - CFD (RANS)
     - Slow (HPC)
     - Detailed flow analysis

All tools use the common **WindIO** schema, enabling seamless comparison across fidelities.

----

Quick Example
-------------

.. tab-set::

    .. tab-item:: CLI

        .. code-block:: console

            # Run with automatic tool selection from YAML
            wifa system.yaml

            # Run with specific tool
            wifa_foxes system.yaml
            wifa_pywake system.yaml

    .. tab-item:: Python

        .. code-block:: python

            from wifa.main_api import run_api

            # Run simulation (tool selected from YAML)
            run_api("path/to/system.yaml")

            # Or use tool-specific API
            from wifa.foxes_api import run_foxes
            results = run_foxes("system.yaml", engine="process", n_procs=4)

.. toctree::
   :hidden:

   getting_started/index
   windio
   tools/index
   API
   CLI
   references

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
