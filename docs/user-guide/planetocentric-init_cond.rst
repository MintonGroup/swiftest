#################################
Planetocentric Initial Conditions
#################################

.. rubric:: by David A. Minton

The central body is treated differently than other massive bodies in a Swiftest simulation, owing to its use of symplectic integrators and the kinds of systems that Swiftest is designed to simulate. 
 When adding planetary bodies using the :func:`add_solar_system_body <swiftest.Simulation.add_solar_system_body>`, the largest body in the dataset is used to define the central body of the system. Swiftest uses the JPL/Horizons ephemerides service to fetch initial conditions of a variety of solar system bodies, and by default solar system bodies are initialized in a Sun-centered, ecliptic reference frame. When a body other than the Sun is set to the central body, the initial conditions of all bodies are translated into the central body's frame, but the z-axis of the simulation coordinate frame remains aligned with the ecliptic north pole. When simulating systems with a body other than the Sun, it is more common to align the coordinate system to the central body's rotation pole. This is done by setting the ``align_to_central_body_rotation`` argument to ``True`` when adding bodies to the simulation.

Here we will demonstrate how this works by setting up a simulation of the two moons of Mars, Phobos and Deimos, demonstrating how the central body frame shift is done and how to use the ``align_to_central_body_rotation`` argument.

Mars-Centered Ecliptic Frame
============================

By default, the initial state vectors of any bodies added via the Horizons ephemerides service are translated into the central body's frame (i.e. the largest body in the system). We can demonstrate this by adding Mars, Phobos, and Deimos to a simulation and printing their inclinations and rotation vectors.

.. doctest::
    :options: +ELLIPSIS, +NORMALIZE_WHITESPACE

    >>> import swiftest
    >>> sim = swiftest.Simulation()
    >>> sim.add_solar_system_body(['Mars', 'Phobos', 'Deimos'])
    Fetching ephemerides data for Mars from JPL/Horizons
    Found ephemerides data for Mars Barycenter (4) from JPL/Horizons
    Physical properties found for Mars (499)
    Fetching ephemerides data for Phobos from JPL/Horizons
    Found ephemerides data for Phobos (401) (Phobos) from JPL/Horizons
    Physical properties found for Phobos (401)
    Fetching ephemerides data for Deimos from JPL/Horizons
    Found ephemerides data for Deimos (402) (Deimos) from JPL/Horizons
    Physical properties found for Deimos (402)

    >>> sim.data.sel(time=0,name=['Phobos','Deimos']).inc.values
    array([26.55406151, 24.06273241])

    >>> sim.data.sel(time=0).rot.values
    array([[ 57176.72129971,  -7162.31013346, 114478.65911179],
       [183617.63838933, -15311.68519899, 368719.81695279],
       [ 42015.96468119,  -6035.89004401,  95062.95557518]])

The inclinations of Phobos and Deimos are 26.55 and 24.06 degrees, respectively, which seems rather high. Now let's initialize the simulation with the ``align_to_central_body_rotation`` argument set to ``True`` and see how the inclinations and rotation vectors change.


.. doctest::
    :options: +ELLIPSIS, +NORMALIZE_WHITESPACE

    >>> import swiftest
    >>> sim = swiftest.Simulation()
    >>> sim.add_solar_system_body(['Mars', 'Phobos', 'Deimos'],align_to_central_body_rotation=True)
    Fetching ephemerides data for Mars from JPL/Horizons
    Found ephemerides data for Mars Barycenter (4) from JPL/Horizons
    Physical properties found for Mars (499)
    Fetching ephemerides data for Phobos from JPL/Horizons
    Found ephemerides data for Phobos (401) (Phobos) from JPL/Horizons
    Physical properties found for Phobos (401)
    Fetching ephemerides data for Deimos from JPL/Horizons
    Found ephemerides data for Deimos (402) (Deimos) from JPL/Horizons
    Physical properties found for Deimos (402)

    >>> sim.data.sel(time=0,name=['Phobos','Deimos']).inc.values
    array([1.0753048, 2.6925768])

    >>> sim.data.sel(time=0).rot.values
    array([[-8.71125394e-12,  9.09494702e-13,  1.28163331e+05],
        [-3.81633732e+02,  7.73720290e+03,  4.12121558e+05],
        [-4.89032977e+03, -1.60117170e+02,  1.03994220e+05]])

Now we can see that the inclinations of Phobos and Deimos are 1.08 and 2.69 degrees, respectively, which is much more reasonable. The rotation vectors have also changed, as expected, with the x and y components of the rotation vector of Mars being nearly zero (well within floating point precision). This is because the coordinate frame of the system has been rotated to align the z-axis with the central body's north pole.


Mixing Central Body Frames
==========================

When adding bodies to a simulation, it is possible to mix central body frames. In order to prevent confusion, we demonstrate some impractical examples to illustrate the behavior of the ``align_to_central_body_rotation`` argument. There are two distinct scenarios that can occur when passing this argument. The first is when a new body is added to a simulation that is larger than any previous ones. When this occurs, the large body is set to be the central body, and *all* bodies in the simulation are translated into the new frame and then rotated, if ``align_to_central_body_rotation`` is set to ``True``. The second scenario is when one or more new bodies are added to a simulation that are all smaller than the current central body. In this case, the new body or bodies are translated into the central body's frame, and *only the new bodies are rotated* if ``align_to_central_body_rotation`` is set to ``True``. This is demonstrated in the following examples.

.. doctest::
    :options: +ELLIPSIS, +NORMALIZE_WHITESPACE

    >>> import swiftest
    >>> sim = swiftest.Simulation()
    >>> sim.add_solar_system_body(['Phobos', 'Deimos'])
    Fetching ephemerides data for Phobos from JPL/Horizons
    Found ephemerides data for Phobos (401) (Phobos) from JPL/Horizons
    Physical properties found for Phobos (401)
    Fetching ephemerides data for Deimos from JPL/Horizons
    Found ephemerides data for Deimos (402) (Deimos) from JPL/Horizons
    Physical properties found for Deimos (402) 

    >>> sim.data.name.values
    array(['Phobos', 'Deimos'], dtype='<U32')

    >>> sim.data.sel(time=0).inc.values
    array([      nan, 28.318414])

This sets up a simulation with Phobos and Deimos, with Phobos as the central body. The inclination of Phobos is ``nan`` because it is the central body. Now let's add Mars to the simulation, with ``align_to_central_body_rotation`` set to ``True`` and see how the inclinations change:

    >>> sim.add_solar_system_body('Mars',align_to_central_body_rotation=True)
    Fetching ephemerides data for Mars from JPL/Horizons
    Found ephemerides data for Mars Barycenter (4) from JPL/Horizons
    Physical properties found for Phobos (401)
    Physical properties found for Mars (499)

    >>> sim.data.name.values
    array(['Mars', 'Deimos', 'Phobos'], dtype='<U32')

    >>> sim.data.sel(time=0).inc.values
    array([      nan, 2.6925768, 1.0753048])

Because Mars is now the most massive body in the system, it has replaced Phobos as the central body. Because the central body has changed, the ``align_to_central_body_rotation=True`` argument rotates all bodies in the system to align with Mars's rotation vector. 

Now Let's see what happens when we add a new body to the simulation that is smaller than the current central body and set ``align_to_central_body_rotation`` to ``True``:

    >>> sim = swiftest.Simulation()
    >>> sim.add_solar_system_body(['Mars', 'Phobos'])
    Fetching ephemerides data for Mars from JPL/Horizons
    Found ephemerides data for Mars Barycenter (4) from JPL/Horizons
    Physical properties found for Mars (499)
    Fetching ephemerides data for Phobos from JPL/Horizons
    Found ephemerides data for Phobos (401) (Phobos) from JPL/Horizons
    Physical properties found for Phobos (401)

    >>> sim.data.name.values
    array(['Mars', 'Phobos'], dtype='<U32')

    >>> sim.data.inc.values
    array([[        nan, 26.55406151]])

We have not aligned the pole of Mars when the simulation was initialized, so the inclination of Phobos is its value relative to the ecliptic. Now we will add Deimos and set ``align_to_central_body_rotation`` to ``True``:

    >>> sim.add_solar_system_body('Deimos',align_to_central_body_rotation=True)
    Fetching ephemerides data for Deimos from JPL/Horizons
    Found ephemerides data for Deimos (402) (Deimos) from JPL/Horizons
    Physical properties found for Mars (499)
    Physical properties found for Deimos (402)
    Physical properties found for Mars (499)


    >>> sim.data.name.values
    array(['Mars', 'Phobos', 'Deimos'], dtype='<U32')

    >>> sim.data.sel(time=0).inc.values
    array([        nan, 26.55406151,  2.6925768 ])


We can see that *only* the inclination of Deimos was rotated. Of course, this leads to an inconsistent set of initial conditions for this system. This is why it is important to be careful when using the ``align_to_central_body_rotation`` argument.

.. note:: 
    The ability of Swiftest to shift the central body frame and rotate the system can be convenient for some scenarios, however it is better to minimize the amount of rotations and shifts when generating a set
    of initial conditions. We recommend that the first body added to a simulation should be the central body (it does not matter if it is first in a list of bodies or added separately), and that the ``align_to_central_body_rotation`` argument is either always set to ``True`` or always set to ``False`` for all bodies added to the system, including the central body. 