#########################
Detailed simulation setup
#########################

.. rubric:: by Kaustub Anand and David A. Minton

Here, we will walk you through the basic features of Swiftest and using them in Python. 
This is based on ``swiftest/examples/Basic_Simulation``.

Start with importing Swiftest and other packages we will use in this tutorial. 

.. ipython:: python
   :okwarning:

   import swiftest
   import numpy as np 

Initial Simulation Setup 
===========================

Create a Swiftest Simulation object and clean the simulation directory of any previous Swiftest objects, if any.
Outputs are stored in the ``./simdata`` directory by default. 

.. ipython:: python
   :okwarning:

   sim = swiftest.Simulation()

.. ipython:: python
   :okwarning:
   :suppress:

    import os
    import tempfile
    tmpdir=tempfile.TemporaryDirectory()
    os.environ['OMP_NUM_THREADS'] = '1'
    sim = swiftest.Simulation(simdir=tmpdir.name)


An optional argument can be passed to specify the simulation directory 

.. code-block:: python

   sim = swiftest.Simulation(simdir='/path/to/simdata')

The argument to `simdir` can either be a relative path or an absolute path.  Now that we have a simulation object set up (with default parameters), we can add bodies to the simulation.  The most massive body in the simulation is taken as the central body.

Solar System Bodies
=========================

We can add solar system bodies to the simulation using the :meth:`add_solar_system_body <swiftest.Simulation.add_solar_system_body>` method. 
This method uses JPL Horizons to extract the parameters of a particular body given a name.

.. ipython:: python
   :okwarning:
   
   # Add the modern planets and the Sun using the JPL Horizons Database.
   sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])

We can add other small bodies too. 

.. ipython:: python
   :okwarning:

   # Add in some main belt asteroids
   sim.add_solar_system_body(name=["Ceres","Vesta","Pallas","Hygiea"],id_type="smallbody")

   # Add in some big KBOs
   sim.add_solar_system_body(name=["Pluto","Eris","Haumea","Quaoar"])

   # Add in some Centaurs
   sim.add_solar_system_body(name=["Chiron","Chariklo"])

.. note::
   The :meth:`add_solar_system_body <swiftest.Simulation.add_solar_system_body>` method is designed for ease of use. If you pass the `name` argument alone, it will make a "best guess" as to which body to retrieve, using the `astroquery.jplhorizons` module. If you want more control over which body to pass to JPL Horizons, you can supply the optional argument `ephemeris_id` in addition to `name`. The string argument passed to `name` is then used internally by Swiftest to identify the body, but the query to JPL Horizons is made with `ephemeris_id`. Therefore the following two calls are equivalent

.. code-block:: python

      sim.add_solar_system_body(name="Sun")
      sim.add_solar_system_body(name="Sun", ephemeris_id="0")

.. note::
   The arguments `name="Earth"` and `name="Pluto"` are handled as special cases in :meth:`add_solar_system_body <swiftest.Simulation.add_solar_system_body>` due to their unusually massive satellites. When "Earth" (or "Pluto") is requested, then the mass the Moon (Charon) is added to the body mass and the initial conditions are set to the Earth-Moon (Pluto-Charon) barycenter. If you wish to instead request the planet directly, you should pass `ephemeris_id` instead. 

User-defined bodies
=========================

You can add a user-defined body with arbitrary initial conditions using using :meth:`sim.add_body <swiftest.Simulation.add_body>`. This method contains a number of optional arguments, and different combinations of arguments can result in different kinds of bodies. 

- id: This is a unique, positive integer id for the body. Usually you would not pass this argument, as an id will be automatically assigned in the order in which it was added, and the central body is always assigned to be id 0.

- name: This is a unique, string name for the body. If you do not pass this, then the name will be set to "Body{id}"

- rh,vh: These are the position and velocity vectors of the body in Cartesian coordinates. These are used to set the initial conditions for the body when the Simulation is set to `init_cond_format="XV"` (or the equivalent `param["IN_FORM"] = "XV"`).

- a,e,inc,capom,omega,capm: These are used to set the initial osculating orbital elements for the body when the Simulation is set to `init_cond_format="EL"` (or the equivalent `param["IN_FORM"] = "EL"`). 


We will randomize the initial conditions and therefore import the `numpy.random <https://numpy.org/doc/stable/reference/random/index.html#module-numpy.random>`__ module.

.. ipython:: python
   :okwarning:

   from numpy.random import default_rng
   rng = default_rng(seed=123)

Starting with **massive bodies:** 

.. ipython:: python
   :okwarning:

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

Initialize orbital elements and then add the bodies.

.. ipython:: python
   :okwarning:
   
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
Start by defining the position and velocity vectors. Here we define the orbital velocities and scale them by a random value. 

.. code-block:: python

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
Here is an example with orbital elements:

.. ipython:: python
   :okwarning:

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

- When creating the initial Swiftest simulation object

.. code-block:: python

    sim = swiftest.Simulation(integrator = 'symba', tstart=0.0, tstop=1.0e3, dt=0.01, 
                                tstep_out=1.0, dump_cadence=0, compute_conservation_values=True, mtiny=mtiny)
    
- :meth:`sim.set_parameter <swiftest.Simulation.set_parameter>`: Set individual parameters in the simulation. The user can set one or multiple at a time.

.. ipython:: python
   :okwarning:

    sim.set_parameter(tstart=0.0, tstop=1.0e3, dt=0.01, tstep_out=1.0, dump_cadence=0, compute_conservation_values=True, mtiny=mtiny)
    sim.set_parameter(rmin = 0.05)

We now set up the simulation parameters. Here we have a simulation starting from `0.0 y` and running for `1 My = 1e6 years` 
with time steps of `0.01 years`. The timestep should be less than or equal to 1/20 of the orbital period of the innermost body. 

The user can then write the parameters to the `param.in` file by using :meth:`write_param <swiftest.Simulation.write_param>`.
To see the parameters of the simulation, use :meth:`sim.get_parameter <swiftest.Simulation.get_parameter>`.

Running the Simulation
========================

Once everything is set up, we can save the simulation object and then run it.

.. ipython:: python
   :okwarning:

    sim.run()

Once this is finished, you should be able to access the output data stored in the :attr:`~swiftest.Simulation.data` attribute.

.. ipython:: python
  :suppress:

  # Import xarray and set its output to show more lines
  import xarray as xr
  xr.set_options(display_max_rows=50)

.. ipython:: python
   :okwarning:

    sim.data

Or, say, plot the eccentricity history of just the test particles:

.. ipython:: python
   :okwarning:

   @savefig detailed_simulation_e_vs_t_tp.png width=800px
   sim.data['e'].where(sim.data.particle_type == 'Test Particle',drop=True).plot(x='time',hue='name');


Modifying and Removing Bodies
==============================

Modifying the properties of initial conditions bodies and removing them is easily done with :meth:`~swiftest.Simulation.modify_body` and :meth:`~swiftest.Simulation.remove_body`. Any property (other than ``name``, which is the unique identifier) can be modified.


.. code-block:: python

   sim.modify_body(name="TestParticle_01", a=1.0, e=0.1, inc=0.0, capom=0.0, omega=0.0, capm=0.0)

Removing bodies is also straightforward:

.. code-block:: python

   sim.remove_body(name="TestParticle_02")

You can also alter the central body. For instance, if you wanted to use a set of coefficients from the `SHTOOLS library <https://shtools.github.io/SHTOOLS/>`, you could do the following.

.. code-block:: python

   import swiftest
   import pyshtools as pysh

   sim = swiftest.Simulation()
   c_lm = pysh.datasets.Mars.GMM3(lmax = 6).coeffs
   sim.add_solar_system_body(["Mars","Phobos","Deimos"])
   sim.modify_body(name="Mars", c_lm=c_lm)

.. toctree::
   :maxdepth: 2
   :hidden:
