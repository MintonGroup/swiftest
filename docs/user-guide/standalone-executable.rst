###############################
Using the standalone executable
###############################

.. rubric:: by David A. Minton

For many users, it may be desirable to separate the generation of initial conditions from the actual simulation. In such cases, 
users may find it useful to use the standalone executable that is built with Swiftest. The standalone executable is a command-line 
tool that can be used to run simulations from a set of initial conditions and a configuration file. The configuration file is a 
simple text file that contains the parameters for the simulation, such as the total time of the simulation, the time step size, and 
the output cadence. By default this file is called `param.in`. Running the executable is similar to using the older 
`Swift <https://www.boulder.swri.edu/~hal/swift.html>`_ and `Swifter <https://www.boulder.swri.edu/swifter/>`_ packages.


Create a basic solar system simulation
=======================================

First we will create a basic solar system simulation

.. ipython:: python
  :okwarning:

  import swiftest
  sim = swiftest.Simulation(tstop=1.0, dt=0.01, tstep_out=0.1)
  sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])

Instead of running it, we will simply save the current state of the system to a file.

.. ipython:: python
  :okwarning:

  sim.save()

Outputs are stored in the ``./simdata`` directory by default. There should be two files created: `param.in` and `init_cond.nc`. 
These files contain the parameters for the simulation and the initial conditions of the simulation. Now to run the simulation 
from the terminal, simply navigated to the `simdata` directory and run the executable.

.. ipython:: python
  :okwarning:
  
  %%bash
  cd simdata
  swiftest whm param.in 

This will exectute the simulation using the basic Wisdom-Holman symplectic map integrator and report the progress to the terminal. 
For a more streamlined output that shows a progress bar and saves the standard terminal output to a file called `swiftest.log` you can instead run

.. code-block:: bash

  $ swiftest whm param.in progress

Once the simulation is done, the results can be read in from Python

.. ipython:: python
  :okwarning:

  import swiftest
  sim = swiftest.Simulation(read_data=True)

Passing the argument `read_data=True` informs the Simulation class to read in a pre-existing data file rather than starting a new simulation. 

