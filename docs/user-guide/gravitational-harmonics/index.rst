##########################
Gravitational Harmonics
##########################

Here, we show how to use Swiftest's Gravitational Harmonics capabilities. This is based on ``/spherical_harmonics_cb`` 
in ``swiftest/examples``. Swiftest uses `SHTOOLS <https://shtools.github.io/SHTOOLS/>`__ to compute the gravitational 
harmonics coefficients for a given body and calculate it's associated acceleration kick. The conventions and formulae used 
to calculate the additional kick are described `here <https://sseh.uchicago.edu/doc/Weiczorek_2015.pdf>`__. The gravitational
potential is given by the following equation:

.. math::

    U(r) = \frac{GM}{r} \sum_{l=0}^{\infty} \sum_{m=-l}^{l} \left( \frac{R_0}{r} \right)^l C_{lm} Y_{lm} (\theta, \phi) \label{grav_pot}

Gravitational potential :math:`U` at a point :math:`\vec{r}`; :math:`\theta` is the polar angle; :math:`\phi` is the azimuthal angle; 
:math:`R_0` is the central body radius; :math:`G` is the gravitational constant; :math:`Y_{lm}` is the spherical harmonic function at 
degree :math:`l` and order :math:`m`; :math:`C_{lm}` is the corresponding coefficient.
 
Gravitational Harmonics coefficients
=====================================

Swiftest adopts the  :math:`4\pi` or geodesy normalization for the gravitational harmonics coefficients described 
in `Weiczorek et al. (2015) <https://sseh.uchicago.edu/doc/Weiczorek_2015.pdf>`__.

The coefficients can be computed in a number of ways: 

- Using the axes measurements of the body. (:func:`clm_from_ellipsoid <swiftest.shgrav.clm_from_ellipsoid>`)
- Using a surface relief grid (:func:`clm_from_relief <swiftest.shgrav.clm_from_relief>`). 
   - Note: This function is still in development and may not work as expected.
- Manually entering the coefficients when adding the central body. (:func:`add_body <swiftest.Simulation.add_body>`)
..    - Ensure to correctly normalize the coefficients if manually entering them. 


Computing coefficients from axes measurements
===============================================

Given the axes measurements of a body, the gravitational harmonics coefficients can be computed in a straightforward 
manner. Let's start with setting up the simulation object with units of `km`, `days`, and `kg`. ::
    
    import swiftest

    sim = swiftest.Simulation(DU2M = 1e3, TU = 'd', MU = 'kg', integrator = 'symba')
    sim.clean() 

Define the central body parameters. Here we use Chariklo as an example body and refer to Jacobi Ellipsoid model from 
`Leiva et al. (2017) <https://iopscience.iop.org/article/10.3847/1538-3881/aa8956>`__ for the axes measurements. ::

    cb_mass = 6.1e18
    cb_radius = 123
    cb_a = 157 
    cb_b = 139 
    cb_c = 86 
    cb_volume = 4.0 / 3 * np.pi * cb_radius**3 
    cb_density = cb_mass / cb_volume 
    cb_T_rotation = 7.004 # hours
    cb_T_rotation/= 24.0 # converting to julian days (TU)
    cb_rot = [[0, 0, 360.0 / cb_T_rotation]] # degrees/TU

Once the central body parameters are defined, we can compute the gravitational harmonics coefficients (:math:`C_{lm}`). 
The output coefficients are already correctly normalized. ::

    c_lm, cb_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, c = cb_c, lmax = 6, lref_radius = True)

Add the central body to the simulation. ::

    sim.add_body(name = 'Chariklo', mass = cb_mass, rot = cb_rot, radius = cb_radius, c_lm = c_lm)

Set the parameters for the simulation and run. ::

    sim.set_parameter(tstart=0.0, tstop=10.0, dt=0.01, istep_out=10, dump_cadence=0, compute_conservation_values=True, mtiny=mtiny)

    # Run the simulation
    sim.run()

Setting a reference radius for the coefficients
==================================================

The coefficients can be computed with respect to a reference radius. This is useful when the user wants to explicitly set the reference radius. ::

    c_lm, cb_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, c = cb_c, lmax = 6, lref_radius = True, ref_radius = cb_radius)




.. .. toctree::
..    :maxdepth: 2
..    :hidden: