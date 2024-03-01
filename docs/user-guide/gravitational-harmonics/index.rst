##########################
Gravitational Harmonics
##########################

Here, we show how to use Swiftest's Gravitational Harmonics capabilities. This is based on ``/spherical_harmonics_cb`` 
in ``swiftest/examples``. Swiftest uses `SHTOOLS <https://shtools.github.io/SHTOOLS/>`__ to compute the gravitational 
harmonics coefficients for a given body and calculate it's associated acceleration kick. The additional acceleration 
kick is based on the gravitaional potential described `here <https://sseh.uchicago.edu/doc/Weiczorek_2015.pdf>`__.

..math::

    U(r) = \frac{GM}{r} \sum_{l=0}^{\infty} \sum_{m=-l}^{l} \left( \frac{R_0}{r} \right)^l C_{lm} Y_{lm} (\theta, \phi) \label{grav_pot}

* Gravitational potential:math:`U` at a point:math:`\Vec{r}`;:math:`\theta` is the polar angle;:math:`\phi` is the azimuthal angle;:math:`R_0` is the central body 
radius;:math:`G` is the gravitational constant;:math:`Y_{lm}` is the spherical harmonic function at degree:math:`l` and order:math:`m`;:math:`C_{lm}` is the corresponding coefficient.


Gravitational Harmonics coefficients
=====================================

Swiftest adopts the :math:`4\pi` or geodesy normalization for the gravitational harmonics coefficients.
The coefficients can be computed in a number of ways: 

- Using the axes measurements of the body ( :func:`clm_from_ellipsoid <swiftest.shgrav.clm_from_ellipsoid>`)
  
- Using a surface relief grid ( :func:`clm_from_relief <swiftest.shgrav.clm_from_relief>`). 
  
   - This function is still in development and may not work as expected.
  
- Manually entering the coefficients when adding the central body ( :func:`add_body <swiftest.Simulation.add_body>`)



Computing coefficients from axes measurements
===============================================

.. .. toctree::
..    :maxdepth: 2
..    :hidden: