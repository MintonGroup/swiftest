#################################
Planetocentric Initial Conditions
#################################

.. rubric:: by David A. Minton

The central body is treated differently than other massive bodies in a Swiftest simulation, owing to its use of symplectic integrators and the kinds of systems that Swiftest is designed to simulate. 
 When adding planetary bodies using the :func:`add_solar_system_body <swiftest.Simulation.add_solar_system_body>`, the largest body in
the dataset is used to define the central body of the system. Swiftest uses the JPL/Horizons ephemerides service to fetch initial conditions of a variety of solar system bodies, and by default solar system bodies are initialized in a Sun-centered, ecliptic reference frame. When a body other than the Sun is used as the central body, 
the initial state vectors and orbital elements of all other bodies are shifted into that body's frame, with one exception: the z-axis of the simulation's coordinate frame is aligned to the ecliptic pole of the solar system unless the user specifies otherwise by passing the ``align_to_central_body_rotation=True`` as an argument.

Here we will demonstrate how this works by setting up a simulation of the two moons of Mars, Phobos and Deimos, demonstrating how the central body frame shift is done and how to use the ``align_to_central_body_rotation`` argument.

