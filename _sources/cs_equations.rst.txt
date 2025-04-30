Governing Equations for 3D atmospheric flow
-------------------------------------------

Incompressible Reynolds Averaged Navier Stokes (RANS) mass and momentum conservation equations with anelastic approximation are solved numerically for wind farm simulations

.. math::
    \dfrac{\partial \rho u_i}{\partial x_i}

.. math::
    \dfrac{\partial \rho u_i}{\partial t} + \dfrac{\partial \rho u_i u_j}{\partial x_j} = -\dfrac{\partial p}{\partial x_i} + \dfrac{\partial}{\partial x_i}\left[ \mu^t \left(\dfrac{\partial u_i}{\partial x_j} + \dfrac{\partial u_j}{\partial x_i} - \dfrac{2}{3}\dfrac{\partial u_k}{\partial x_k}\delta_{ij}\right)\right] + \rho g_i

where :math:`\boldsymbol{u}=(u,v,w)^T` is the velocity, :math:`\boldsymbol{x}=(x,y,z)^T` is the position vector, :math:`\rho` is the air density, :math:`p` is the pressure, :math:`\boldsymbol{\delta}` is the Kronecker symbol and :math:`\boldsymbol{g}` is the acceleration of gravity. The turbulent viscosity :math:`\mu^t` is modeled using a :math:`k-\varepsilon` turbulence model with linear production.

Depending on the input parameters and modeling choices, a transport equation for the potential temperature :math:`\theta` is solved

.. math::
   \dfrac{\partial \rho \theta}{\partial t} + \dfrac{\partial \rho u_i \theta}{\partial x_i} = \dfrac{\partial }{\partial x_i}\left[\dfrac{\mu^t}{\sigma_t}\dfrac{\partial \theta}{\partial x_i}\right]

with :math:`\sigma=0.7` the thermal Schmidt number.

More details about the theoretical framework and numerical methods can be found in `code_saturne's online documentation <https://www.code-saturne.org/cms/web/documentation/v80/>`_.


