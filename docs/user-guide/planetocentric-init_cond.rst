#################################
Planetocentric Initial Conditions
#################################

.. rubric:: by David A. Minton

The central body is treated differently than other massive bodies in a Swiftest simulation, owing to its use of symplectic integrators and the kinds of systems that Swiftest is designed to simulate. 
 When adding planetary bodies using the :func:`add_solar_system_body <swiftest.Simulation.add_solar_system_body>`, the largest body in the dataset is used to define the central body of the system. Swiftest uses the JPL/Horizons ephemerides service to fetch initial conditions of a variety of solar system bodies, and by default solar system bodies are initialized in a Sun-centered, ecliptic reference frame. When a body other than the Sun is used as the central body,

Here we will demonstrate how this works by setting up a simulation of the two moons of Mars, Phobos and Deimos, demonstrating how the central body frame shift is done and how to use the ``align_to_central_body_rotation`` argument.

Mars-Centered Ecliptic Frame
============================

By default, the initial state vectors of any bodies added via the Horizons ephemerides service are translated into the central body's frame (i.e. the largest body in the system). We can demonstrate this by adding Phobos and Deimos to a simulation of Mars and then plotting their inclination.

.. code-block:: python

    import swiftest
    sim = swiftest.Simulation()
    sim.add_solar_system_body(['Mars', 'Phobos', 'Deimos'])
    print(sim.init_cond.inc.values)

.. testoutput::

    Fetching ephemerides data for Mars from JPL/Horizons
    Found ephemerides data for Mars Barycenter (4) from JPL/Horizons
    Physical properties found for Mars (499)
    Fetching ephemerides data for Phobos from JPL/Horizons
    Found ephemerides data for Phobos (401) (Phobos) from JPL/Horizons
    Physical properties found for Phobos (401)
    Fetching ephemerides data for Deimos from JPL/Horizons
    Found ephemerides data for Deimos (402) (Deimos) from JPL/Horizons
    Physical properties found for Deimos (402)

    array([        nan, 26.55406151, 24.06273241])
