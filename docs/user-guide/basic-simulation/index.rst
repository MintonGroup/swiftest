#################
Basic Simulation
#################

Here, we will walk you through the basic features of Swiftest and using them in Python. 
This is based on ``/Basic_Simulation`` in ``swiftest/examples``.

Start with importing Swiftest. ::
    
    import swiftest

Initial Simulation Setup 
===========================

Create a Swiftest :class:`Simulation <swiftest.Simulation>` object.
Outputs are stored in the ``./simdata`` directory by default. ::

   sim = swiftest.Simulation()

Now that we have a simulation object set up (with default parameters), we can add bodies to the simulation. 
The biggest body in the simulation is taken as the central body. Swiftest sets a Simulation object up with a set of default parameters, 
including a default unit system of AU, y, and solar masses.

Solar System Bodies
=========================

We can add solar system bodies to the simulation using the :func:`add_solar_system_body <swiftest.Simulation.add_solar_system_body>` method. 
This method uses JPL Horizons to extract the parameters. ::
   
   # Add the modern planets and the Sun using the JPL Horizons Database.
   sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])

Add other small bodies too: ::

   # Add in some main belt asteroids
   sim.add_solar_system_body(name=["Ceres","Vesta","Pallas","Hygiea"],id_type="smallbody")

   # Add in some big KBOs and Centaurs
   sim.add_solar_system_body(name=["Pluto","Eris","Haumea","Quaoar", "Chiron","Chariklo"])

Running the Simulation
========================

We now can some simulation parameters using the :func:`set_parameter <swiftest.Simulation.set_parameter>` method. 
Here we have a simulation that runs for 1 My a step size of 0.01 y. We will also save the system every 1000 y and wait until the end of the simulation to write the simulation data to file using the ``dump_cadence=0`` argument ::

    sim.set_parameter(tstop=1.0e6, tstep_out=1e3, dt=0.01, dump_cadence=0)

Once everything is set up, we call the :func:`run <swiftest.Simulation.run>` method to integrate the system forward in time::

    sim.run()

Swiftest is relatively flexible with arguments. You can pass the parameters in when initializing the simulation object, or even later when running.
So the following are all equivalent::

    sim = swiftest.Simulation(tstop=1.0e6, tstep_out=1e3, dt=0.01, dump_cadence=0)
    sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])
    sim.run()

    sim = swiftest.Simulation()
    sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])
    sim.set_parameter(tstop=1.0e6, tstep_out=1e3, dt=0.01, dump_cadence=0)
    sim.run()

    sim = swiftest.Simulation()
    sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])
    sim.run(tstop=1.0e6, tstep_out=1e3, dt=0.01, dump_cadence=0)


Analayzing Simulation Output
=============================

Once a simulation has been run, its output data is stored in the ``./simdata`` directory. The main data is stored in a file with a 
default name of ``data.nc``, which is a netCDF file. It is read in and stored as an `Xarray Dataset <https://docs.xarray.dev/en/stable/>`__ object in the ``sim.data`` attribute.


.. .. toctree::
..    :maxdepth: 2
..    :hidden:
