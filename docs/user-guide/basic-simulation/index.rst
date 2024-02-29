#################
Basic Simulation
#################

Here, we will walk you through the basic features of Swiftest and using them in Python. 
This is based on ``/Basic_Simulation`` in ``swiftest/examples``.

Start with importing Swiftest and other packages we will use in this tutorial. ::
    
    import swiftest

Initial Simulation Setup 
===========================

Create a Swiftest Simulation object.
Outputs are stored in the ``/simdata`` directory by default. ::

   sim = swiftest.Simulation()

Now that we have a simulation object set up (with default parameters), we can add bodies to the simulation. 
The biggest body in the simulation is taken as the central body.

Solar System Bodies
=========================

We can add solar system bodies to the simulation using the ``add_solar_system_body`` method. 
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

We now set up the simulation parameters. Here we have a simulation starting from :math: `0.0 y` and running for :math: `1 My = 1e6 years` 
with time steps of :math: `0.01 years`. ::

    sim.set_parameter(tstart=0.0, tstop=1.0e6, dt=0.01)

Once everything is set up, we can save the simulation object and then run it: ::

    sim.save()
    sim.run()

.. .. toctree::
..    :maxdepth: 2
..    :hidden:
