.. WIFA documentation master file, created by
   sphinx-quickstart on Thu Apr 11 11:03:12 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. image:: ../img/Logo_FLOW.png
    :align: center
    :width: 150

=====================================
Welcome to WIFA's documentation!
=====================================

WIFA integrates different tools (PyWake, foxes, wayve, code_saturne) of different fidelities to model the atmospheric flow of a wind farm, through a Python framework. This includes engineering wake models, an Atmospheric Perturbation Model (APM), and 3D CFD, all of which allow to estimate the power production of wind farms. These tools are interfaced using a common input and output schema: the WindIO ontology.

The objective of WIFA is to make the use of these models harmonized and easier, although they are initially different both in formulation and use. In particular, this would allow chaining to other tools (e.g. load models).

The presented API is a culmination of close collaboration within the FLOW project. Each supported tool is a product of the institutions represented within the project.

License
-------


Contents
--------
  .. toctree::
      :maxdepth: 2

      installation

  .. toctree::
      :maxdepth: 2

      examples

  .. toctree::
      :maxdepth: 3

      API

  .. toctree::
      :maxdepth: 3

      references


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
