Turbines representation
-----------------------
The presence of wind turbines is accounted for using the Actuator Disk Model (ADM), i.e. a uniform and steady body force (thrust) locally applied for each turbine.

The thrust and power coefficients are in general provided as a function of the free stream velocity :math:`u_0`. In code_saturne and for CFD codes in general, the definition of the latter is unclear in the presence of wake effects. Hence, it is more convenient to compute the thrust force and power as a function of the velocity at disk position :math:`u_D`. In this case, the thrust force formulation is modified to account for this new velocity.

Using the induction coefficient definition :math:`a=1-u_D/u_0`, the axial thrust force can be written:

.. math::
   F_z = \dfrac{1}{2} \rho A u_{D}^2\frac{c_t(u_{0})}{(1-a)^2} = \dfrac{1}{2} \rho A u_{D}^2 c_t^*

where :math:`\rho` is the density, :math:`A` the actuator disk area and :math:`c_t` and :math:`c_t^*=c_t/(1-a)^2` thrust coefficients. Similarly, the Power can be reformulated as:

.. math::
   P(u_0) = \dfrac{1}{2} \rho A u_{D}^3\frac{c_p(u_{0})}{(1-a)^3} = \dfrac{1}{2} \rho A u_{D}^3c_p^*

with :math:`c_p` and :math:`c_p^*=c_p/(1-a)^3` the power coefficients.
