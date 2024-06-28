###################################################
Using SHTOOLS to model non-spherical central bodies
###################################################

.. rubric:: by Kaustub Anand

Here, we show how to use Swiftest's Gravitational Harmonics capabilities. This is based on ``spherical_harmonics_cb`` 
in ``swiftest/examples``. Swiftest uses `SHTOOLS <https://shtools.github.io/SHTOOLS/>`__ to compute the gravitational 
harmonics coefficients for a given body and calculate it's associated acceleration kick. The conventions and formulae used here 
are described in `Weiczorek et al. (2015) <https://sseh.uchicago.edu/doc/Weiczorek_2015.pdf>`__. The gravitational
potential is given by the following equation:

.. math::

    U(r) = \frac{GM}{r} \sum_{l=0}^{\infty} \sum_{m=-l}^{l} \left( \frac{R_0}{r} \right)^l C_{lm} Y_{lm} (\theta, \phi) \label{grav_pot}

Gravitational potential :math:`U` at a point :math:`\vec{r}`; :math:`\theta` is the polar angle; :math:`\phi` is the azimuthal angle; 
:math:`R_0` is the central body radius; :math:`G` is the gravitational constant; :math:`Y_{lm}` is the spherical harmonic function at 
degree :math:`l` and order :math:`m`; :math:`C_{lm}` is the corresponding coefficient.

Set up a Simulation
====================

Let's start with setting up the simulation object with units of `km`, `days`, and `kg`. 

.. code-block:: python
    
    import swiftest

    sim = swiftest.Simulation(DU2M = 1e3, TU = 'd', MU = 'kg', integrator = 'symba')
 
Gravitational Harmonics Coefficients
=====================================

Swiftest adopts the  :math:`4\pi` or geodesy normalization for the gravitational harmonics coefficients described 
in `Weiczorek et al. (2015) <https://sseh.uchicago.edu/doc/Weiczorek_2015.pdf>`__. 

The coefficients can be computed in a number of ways: 

- Using the axes measurements of the body. (:func:`clm_from_ellipsoid <swiftest.shgrav.clm_from_ellipsoid>`)

- Using a surface relief grid (:func:`clm_from_relief <swiftest.shgrav.clm_from_relief>`). *Note: This function is still in development and may not work as expected.*

- Manually entering the coefficients when adding the central body. (:func:`add_body <swiftest.Simulation.add_body>`)

Computing Coefficients from Axes Measurements
===============================================

Given the axes measurements of a body, the gravitational harmonics coefficients can be computed in a straightforward 
manner. Here we use Chariklo as an example body and refer to Jacobi Ellipsoid model from 
`Leiva et al. (2017) <https://iopscience.iop.org/article/10.3847/1538-3881/aa8956>`__ for the axes measurements. 

.. code-block:: python

    # Define the central body parameters. 
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
The output coefficients are already correctly normalized. 

.. code-block:: python

    c_lm, cb_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, c = cb_c, lref_radius = True)

*Note: Here we set the reference radius flag to* `True` *and ask the function to return the reference radius. More in the 
additional capabilities section below. The maximum degree is set to 6 by default to ensure computational efficiency.*

Add the central body to the simulation along with the coefficients. 

.. code-block:: python

    sim.add_body(name = 'Chariklo', mass = cb_mass, rot = cb_rot, radius = cb_radius, c_lm = c_lm)

If the :math:`J_{2}` and :math:`J_{4}` terms are passed as well, Swiftest ignores them and uses the :math:`C_{lm}` terms instead.
Now the user can set up the rest of the simulation as usual. 

.. code-block:: python

    # add other bodies and set simulation parameters
    .
    .
    .

    # run the simulation
    sim.run()

Manually Adding the Coefficients
================================

If the user already has the coefficients, they can be added to the central body directly. Ensure that they are correctly normalized and 
the right shape. The dimensions of ``c_lm`` is ``[sign, l, m]`` where: 

- ``sign`` indicates coefficients of positive (``[0]``) and negative (``[1]``) ``m`` and is of length 2. 
- The dimension ``l`` corresponds to the degree of the Spherical Harmonic and is of length :math:`l_{max} + 1`.
- The dimension ``m`` corresponds to the order of the Spherical Harmonic and is of length :math:`l_{max} + 1`.

:math:`l_{max}` is the highest order of the coefficients. 

.. code-block:: python

    c_lm = ..... # defined by the user
    sim.add_body(name = 'Body', mass = cb_mass, rot = cb_rot, radius = cb_radius, c_lm = c_lm)

Additional Capabilities of Swiftest's Coefficient Generator Functions
===========================================================================================

The output from :func:`~swiftest.shgrav.clm_from_ellipsoid` and :func:`~swiftest.shgrav.clm_from_relief`
can be customised to the user's needs. Here we show some of the additional capabilities of these functions.

Setting a Reference Radius for the Coefficients
-------------------------------------------------

The coefficients are computed with respect to a reference radius. `SHTOOLS <https://shtools.github.io/SHTOOLS/>`__ calculates it's own radius from 
the axes passed, but there are different ways to calculate the reference radius for non-spherical bodies in the literature. As a result, Swiftest allows 
the user to explicitly set a reference radius (``ref_radius``) which scales the coefficients accordingly. This is particularly useful when a 
specific radius is desired.

This is done by setting ``lref_radius = True`` and passing a ``ref_radius``. Here we pass the Central Body radius (``cb_radius``) manually set earlier as 
the reference. 

.. code-block:: python

    c_lm, ref_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, c = cb_c, lref_radius = True, ref_radius = cb_radius)

When ``lref_radius == True``, it tells the function to return the reference radius used to calculate the 
coefficients and look for any reference radius (``ref_radius``) passed. If no reference radius is passed, the function returns the radius calculated
internally. 

.. code-block:: python
        
    c_lm, ref_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, c = cb_c, lref_radius = True)

By default, ``lref_radius`` is set to ``False``. In this case, the function only returns the coefficients. 

.. code-block:: python

    c_lm = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, c = cb_c)

We recommend extracting the ``ref_radius`` from the function output and using it when adding the central body to the simulation.

Combinations of Principal Axes
-------------------------------

The user can pass any combinations of the principal axes (``a``, ``b``, and ``c``) with ``a`` being the only required one. This is particularly 
useful for cases like oblate spheroids (:math:`a = b \neq c`). For example, the following statements are equivalent: 

.. code-block:: python
    
    c_lm, ref_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, c = cb_c, lref_radius = True)

    c_lm, ref_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_a, c = cb_c, lref_radius = True)

    c_lm, ref_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, c = cb_c, lref_radius = True)

For bodies with :math:`a \neq b = c`, the following statements are equivalent: 

.. code-block:: python
    
    c_lm, ref_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, c = cb_c, lref_radius = True)

    c_lm, ref_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, c = cb_b, lref_radius = True)

    c_lm, ref_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, lref_radius = True)


Setting the Maximum Degree :math:`l`
-------------------------------------

The user can set the maximum degree :math:`l` for the coefficients. 

.. code-block:: python

    lmax = 4
    c_lm, ref_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, c = cb_c, lmax = lmax, lref_radius = True)

``lmax`` is by currently capped to 6 to ensure computational efficiency. This is derived from Jean's law by setting the 
characteristic wavelength (:math:`\lambda`) of a harmonic degree (:math:`l`) to the radius (:math:`R`) of the body.

.. math:: 

    \lambda = \frac{2\pi R}{\sqrt{l(l+1)}} 

    \lambda = R \Rightarrow l = 6
