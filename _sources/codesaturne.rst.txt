code_saturne
-------------------------------------
.. image:: ../img/Logo_code_saturne.png
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
