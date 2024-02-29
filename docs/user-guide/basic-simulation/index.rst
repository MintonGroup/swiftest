#################
Basic Simulation
#################

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

We can add solar system bodies to the simulation using the ``add_solar_system_body`` method. 
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

For completeness, let's also add some bodies with user defined parameters using ``sim.add_body()``.
We will randomize the initial conditions and therefore import the ``numpy.random`` module.::

   from numpy.random import default_rng
   rng = default_rng(seed=123)

Starting with massive bodies: ::

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

The process is similar for adding bodies with Cartesian coordinates. Start by defining the position and velocity vectors.
Here we define the orbital velocities and scale them by a random value. ::
   
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



Customising Simulation Parameters
==================================

The user can also specify simulation parameters when creating the Swiftest object as defined in <>, but we will come back to this later.


    


.. .. toctree::
..    :maxdepth: 2
..    :hidden:
