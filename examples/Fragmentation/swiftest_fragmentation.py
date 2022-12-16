"""
 Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest.
 If not, see: https://www.gnu.org/licenses.
"""
#!/usr/bin/env python3
"""
Generates and runs a set of Swiftest input files from initial conditions with the SyMBA integrator. All simulation 
outputs for the disruption case are stored in the /disruption subdirectory. All simulation outputs for the hit and run 
case are stored in the /hitandrun subdirectory. All simulation outputs for the super-catastrophic disruption case are 
stored in the /supercat subdirectory.

Input
------
None.

Output
------
Three subdirectories:
disruption/
hitandrun/
supercat/

Each subdirectory contains:
data.nc             : A NetCDF file containing the simulation output.
init_cond.nc        : A NetCDF file containing the initial conditions for the simulation.
collision_000001.nc : A NetCDF file containing the data for the collision.
encounter_000001.nc : A NetCDF file containing the data for the close encounter.
dump_bin1.nc        : A NetCDF file containing the necessary inputs to restart a simulation from t!=0.
dump_bin2.nc        : A NetCDF file containing the necessary inputs to restart a simulation from t!=0.
dump_param1.in      : An ASCII file containing the necessary parameters to restart a simulation.
dump_param2.in      : An ASCII file containing the necessary parameters to restart a simulation.
fraggle.log         : An ASCII file containing the information of any collisional events that occured.
param.in            : An ASCII file containing the parameters for the simulation.
swiftest.log        : An ASCII file containing the information on the status of the simulation as it runs.

"""
import swiftest
import numpy as np
from numpy.random import default_rng

run_arguments = {"tstart": 0.0, "tstop":1e-5, "dt": 1e-5, "istep_out": 1, "fragmentation":True, "encounter_save":"both", "compute_conservation_values":True,
                 "minimum_fragment_gmass":1.0e-11, "gmtiny":1.0e-11, "output_format":"XVEL", "init_cond_format":"XV"}

# Initialize the simulation object as a variable with arguments.
sim_disruption = swiftest.Simulation(simdir="disruption", **run_arguments)
# Add the Sun using the JPL Horizons Database.
sim_disruption.add_solar_system_body(["Sun"])
# Add a user-defined target body.
sim_disruption.add_body(name="Target", rh=[1.0, -1.807993e-05, 0.0], vh=[-2.562596e-04, 6.280005, 0.0], Gmass=1e-7, radius=7e-6)
# Add a user-defined projectile body.
sim_disruption.add_body(name="Projectile", rh=[1.0, 1.807993e-05, 0.0], vh=[-2.562596e-04, -6.280005, 0.0], Gmass=7e-10, radius=3.25e-6)
# Display the run configuration parameters.
sim_disruption.get_parameter()
# Write the parameters to the param.in
sim_disruption.write_param()
# Run the simulation.
sim_disruption.run()

# Do the same as above for the hit and run case.
sim_hitandrun = swiftest.Simulation(simdir="hitandrun", **run_arguments)
sim_hitandrun.add_solar_system_body(["Sun"])
sim_hitandrun.add_body(name="Target", rh=[1.0, -4.2e-05, 0.0], vh=[0.0, 6.28, 0.0], Gmass=1e-7, radius=7e-6, rot=[0.0, 0.0, 6.0e4])
sim_hitandrun.add_body(name="Projectile", rh=[1.0, 4.2e-05, 0.0], vh=[-1.5, -6.28, 0.0], Gmass=7e-10, radius=3.25e-6, rot=[0.0, 0.0, 1.0e5])
sim_hitandrun.get_parameter()
sim_hitandrun.write_param()
sim_hitandrun.run()

# Do the same as above for the super-catastrophic disruption case.
sim_supercat = swiftest.Simulation(simdir="supercat", **run_arguments)
sim_supercat.add_solar_system_body(["Sun"])
sim_supercat.add_body(name="Target", rh=[1.0, -4.2e-05, 0.0], vh=[0.0, 6.28, 0.0], Gmass=1e-7, radius=7e-6, rot=[0.0, 0.0, -6.0e4])
sim_supercat.add_body(name="Projectile", rh=[1.0, 4.2e-05, 0.0], vh=[1.0, -6.28, 0.0], Gmass=1e-8, radius=3.25e-6, rot=[0.0, 0.0, 1.0e5])
sim_supercat.get_parameter()
sim_supercat.write_param()
sim_supercat.run()
