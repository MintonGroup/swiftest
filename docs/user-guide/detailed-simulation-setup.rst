#####################
Detailed Simulation
#####################

.. rubric:: by Kaustub Anand

Here, we will walk you through the basic features of Swiftest and using them in Python. 
This is based on ``/Basic_Simulation`` in ``swiftest/examples``.

Start with importing Swiftest and other packages we will use in this tutorial. ::
    
    import swiftest
    import numpy as np 

Initial Simulation Setup 
===========================

Create a Swiftest Simulation object and clean the simulation directory of any previous Swiftest objects, if any.
Outputs are stored in the ``/simdata`` directory by default. ::

   sim = swiftest.Simulation()
   sim.clean()

An optional argument can be passed to specify the simulation directory ::

   simdir = '/path/to/simdir'
   sim = swiftest.Simulation(simdir=simdir)
   sim.clean()

Now that we have a simulation object set up (with default parameters), we can add bodies to the simulation. 
The biggest body in the simulation is taken as the central body.

Solar System Bodies
=========================

We can add solar system bodies to the simulation using the :func:`add_solar_system_body <swiftest.Simulation.add_solar_system_body>` method. 
This method uses JPL Horizons to extract the parameters. ::
   
   # Add the modern planets and the Sun using the JPL Horizons Database.
   sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])

We can add other small bodies too. ::

   # Add in some main belt asteroids
   sim.add_solar_system_body(name=["Ceres","Vesta","Pallas","Hygiea"],id_type="smallbody")

   # Add in some big KBOs
   sim.add_solar_system_body(name=["Pluto","Eris","Haumea","Quaoar"])

   # Add in some Centaurs
   sim.add_solar_system_body(name=["Chiron","Chariklo"])

User Defined Bodies
=========================

For completeness, let's also add some bodies with user defined parameters using :func:`sim.add_body <swiftest.Simulation.add_body>`.
We will randomize the initial conditions and therefore import the `numpy.random <https://numpy.org/doc/stable/reference/random/index.html#module-numpy.random>`__ module.::

   from numpy.random import default_rng
   rng = default_rng(seed=123)

Starting with **massive bodies:** ::

   npl = 5 # number of massive bodies
   density_pl  = 3000.0 / (sim.param['MU2KG'] / sim.param['DU2M'] ** 3)
   name_pl     = ["SemiBody_01", "SemiBody_02", "SemiBody_03", "SemiBody_04", "SemiBody_05"]

   M_pl        = np.array([6e20, 8e20, 1e21, 3e21, 5e21]) * sim.KG2MU # mass in simulation units
   R_pl        = np.full(npl, (3 * M_pl/ (4 * np.pi * density_pl)) ** (1.0 / 3.0)) # radius
   Ip_pl       = np.full((npl,3),0.4,) # moment of inertia
   rot_pl      = np.zeros((npl,3)) # initial rotation vector in degrees/TU
   mtiny       = 1.1 * np.max(M_pl) # threshold mass for semi-interacting bodies in SyMBA.

Depending on the simulation parameters, we can add bodies with Orbital Elements or Cartesian Coordinates.

Orbital Elements
-------------------

Initialize orbital elements and then add the bodies. ::
   
   a_pl        = rng.uniform(0.3, 1.5, npl) # semi-major axis
   e_pl        = rng.uniform(0.0, 0.2, npl) # eccentricity
   inc_pl      = rng.uniform(0.0, 10, npl) # inclination (degrees)
   capom_pl    = rng.uniform(0.0, 360.0, npl) # longitude of the ascending node
   omega_pl    = rng.uniform(0.0, 360.0, npl) # argument of periapsis
   capm_pl     = rng.uniform(0.0, 360.0, npl) # mean anomaly

   sim.add_body(name=name_pl, a=a_pl, e=e_pl, inc=inc_pl, capom=capom_pl, omega=omega_pl, capm=capm_pl, mass=M_pl, radius=R_pl,  Ip=Ip_pl, rot=rot_pl)

Cartesian Coordinates
----------------------

The process is similar for adding bodies with cartesian coordinates. However, the parameter `init_cond_format` must be set to `XV` before adding the bodies.
The process of setting parameters is explained in the next section. 
Start by defining the position and velocity vectors. Here we define the orbital velocities and scale them by a random value. ::
   
   # position vectors
   rh_pl = rng.uniform(-5, 5, (npl,3))
   rh_pl_mag = np.linalg.norm(rh_pl, axis=1) # magnitudes of the position vector

   # General velocity vectors

      # define the magnitudes
   velocity_scale = rng.uniform(0.5, 1.5, npl) # scale the orbital velocity
   vh_pl_mag = velocity_scale * np.sqrt(sim.GU * M_pl / rh_pl_mag) # magnitude of the velocity vector

      # initialize the vectors using the position vectors
   vx = rh_pl.T[0] * vh_pl_mag / rh_pl_mag
   vy = rh_pl.T[1] * vh_pl_mag / rh_pl_mag
   vz = rh_pl.T[2] * vh_pl_mag / rh_pl_mag
   
      # rotate the velocity vectors to the XY plane for orbital motion
   vh_pl = np.array([vx, vy, vz]).T
   vh_pl = np.cross(vh_pl, np.array([0,0,1])) # velocity vectors

   sim.add_body(name=name_pl, rh=rh_pl, vh=vh_pl, mass=M_pl, radius=R_pl,  Ip=Ip_pl, rot=rot_pl)

The process is similar for **test particles**. They only need the orbital elements or the cartesian coordinates. 
Here is an example with orbital elements: ::

    # Add 10 user-defined test particles.
    ntp = 10

    name_tp     = ["TestParticle_01", "TestParticle_02", "TestParticle_03", "TestParticle_04", "TestParticle_05", "TestParticle_06", "TestParticle_07", "TestParticle_08", "TestParticle_09", "TestParticle_10"]
    a_tp        = rng.uniform(0.3, 1.5, ntp)
    e_tp        = rng.uniform(0.0, 0.2, ntp)
    inc_tp      = rng.uniform(0.0, 10, ntp)
    capom_tp    = rng.uniform(0.0, 360.0, ntp)
    omega_tp    = rng.uniform(0.0, 360.0, ntp)
    capm_tp     = rng.uniform(0.0, 360.0, ntp)

    sim.add_body(name=name_tp, a=a_tp, e=e_tp, inc=inc_tp, capom=capom_tp, omega=omega_tp, capm=capm_tp)


Customising Simulation Parameters
==================================

Now that we have added the bodies, we can set the simulation parameters. ``tstop`` and ``dt`` need to be set before running the simulation.
This can be done in multiple ways:

- When creating the initial Swiftest simulation object ::
    
    sim = swiftest.Simulation(simdir = simdir, integrator = 'symba', init_cond_format = 'EL', tstart=0.0, tstop=1.0e6, dt=0.01, 
                                istep_out=100, dump_cadence=0, compute_conservation_values=True, mtiny=mtiny)
    
- :func:`sim.set_parameter <swiftest.Simulation.set_parameter>`: Set individual parameters in the simulation. The user can set one or multiple at a time. ::

    sim.set_parameter(tstart=0.0, tstop=1.0e6, dt=0.01, istep_out=100, dump_cadence=0, compute_conservation_values=True, mtiny=mtiny)
    sim.set_parameter(rmin = 0.05)

We now set up the simulation parameters. Here we have a simulation starting from `0.0 y` and running for `1 My = 1e6 years` 
with time steps of `0.01 years`. The timestep should be less than or equal to 1/10 of the orbital period of the innermost body. 

The user can then write the parameters to the `param.in` file by using :func:`write_param <swiftest.Simulation.write_param>`.
To see the parameters of the simulation, use :func:`sim.get_parameter <swiftest.Simulation.get_parameter>`.

Running the Simulation
========================

Once everything is set up, we can save the simulation object and then run it: ::

    sim.save()
    sim.run()

.. .. toctree::
..    :maxdepth: 2
..    :hidden:
