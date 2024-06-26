#!/usr/bin/env python3

"""
 Copyright 2024 - The Minton Group at Purdue University
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest. 
 If not, see: https://www.gnu.org/licenses. 
"""

"""
Generates and runs a set of Swiftest input files from initial conditions with the SyMBA integrator. All simulation 
outputs are stored in the /simdata subdirectory.

Input
------
None.

Output
------
collisions.log   : An ASCII file containing the information of any collisional events that occured.
collisions.nc    : A NetCDF file containing the collision output.
data.nc          : A NetCDF file containing the simulation output.
init_cond.nc     : A NetCDF file containing the initial conditions for the simulation.
param.00...0.in  : A series of parameter input files containing the parameters for the simulation at every output stage.
param.in         : An ASCII file containing the parameters for the simulation.
param.restart.in : An ASCII file containing the parameters for the simulation at the last output. 
swiftest.log     : An ASCII file containing the information on the status of the simulation as it runs.
"""

import swiftest
import numpy as np

# Initialize the simulation object as a variable. Arguments may be defined here or through the sim.run() method.
sim = swiftest.Simulation(simdir='massive_sundiver',compute_conservation_values=True, rotation=True)

# Add the modern planets and the Sun using the JPL Horizons Database.
sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune","Pluto"])

density  = 3000.0 * sim.KG2MU / sim.M2DU**3

# Make a body with a periapsis inside the Sun's radius
q = 0.01 * swiftest.RSun * sim.M2DU
a = 0.1
e = 1.0 - q / a
M = 2e0 * swiftest.MEarth * sim.KG2MU
R = (3 * M  / (4 * np.pi * density)) ** (1.0 / 3.0)
rot = 4 * sim.init_cond.sel(name="Earth")['rot']
sim.add_body(name="Sundiver", a=a, e=e, inc=0.0, capom=0.0, omega=0.0, capm=180.0, mass=M, radius=R, Ip=[0.4,0.4,0.4], rot=rot)
sim.get_parameter()

# Run the simulation. Arguments may be defined here or thorugh the swiftest.Simulation() method.
sim.run(tstart=0.0, tstop=5e-2, dt=0.0001, istep_out=1, dump_cadence=0, integrator="rmvs")

sim = swiftest.Simulation(simdir='massless_sundiver',compute_conservation_values=False, integrator="rmvs")

sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune","Pluto"])
sim.add_body(name="Sundiver", a=a, e=e, inc=0.0, capom=0.0, omega=0.0, capm=180.0)
sim.get_parameter()

# Run the simulation. Arguments may be defined here or thorugh the swiftest.Simulation() method.
sim.run(tstart=0.0, tstop=5e-2, dt=0.0001, istep_out=1, dump_cadence=0)


