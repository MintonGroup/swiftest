#!/usr/bin/env python3

"""
 Copyright 2023 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
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
param.in         : An ASCII file containing the inital parameters for the simulation.
param.restart.in : An ASCII file containing the parameters for the simulation at the last output. 
swiftest.log     : An ASCII file containing the information on the status of the simulation as it runs.
"""

import swiftest
import numpy as np
from numpy.random import default_rng

# Initialize the simulation object as a variable. Arguments may be defined here or through the sim.run() method.
sim = swiftest.Simulation()
rng = default_rng(seed=123)

# Add the modern planets and the Sun using the JPL Horizons Database.
sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])

# Add in some main belt asteroids
sim.add_solar_system_body(name=["Ceres","Vesta","Pallas","Hygiea"],id_type="smallbody")

# Add in some big KBOs
sim.add_solar_system_body(name=["Pluto","Eris","Sedna","Haumea","Makemake","Quaoar","Orcus","Gonggong","Salacia"])

# Add 5 user-defined massive bodies.
npl         = 5
density_pl  = 3000.0 / (sim.param['MU2KG'] / sim.param['DU2M'] ** 3)

name_pl     = ["SemiBody_01", "SemiBody_02", "SemiBody_03", "SemiBody_04", "SemiBody_05"]
a_pl        = rng.uniform(0.3, 1.5, npl)
e_pl        = rng.uniform(0.0, 0.2, npl)
inc_pl      = rng.uniform(0.0, 10, npl)
capom_pl    = rng.uniform(0.0, 360.0, npl)
omega_pl    = rng.uniform(0.0, 360.0, npl)
capm_pl     = rng.uniform(0.0, 360.0, npl)
M_pl        = np.array([6e20, 8e20, 1e21, 3e21, 5e21]) * sim.KG2MU
R_pl        = np.full(npl, (3 * M_pl/ (4 * np.pi * density_pl)) ** (1.0 / 3.0))
Ip_pl       = np.full((npl,3),0.4,)
rot_pl      = np.zeros((npl,3))
mtiny       = 1.1 * np.max(M_pl)

sim.add_body(name=name_pl, a=a_pl, e=e_pl, inc=inc_pl, capom=capom_pl, omega=omega_pl, capm=capm_pl, mass=M_pl, radius=R_pl,  Ip=Ip_pl, rot=rot_pl)

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
sim.set_parameter(tstart=0.0, tstop=1.0e3, dt=0.01, istep_out=100, dump_cadence=0, compute_conservation_values=True, mtiny=mtiny)
# Display the run configuration parameters.
sim.get_parameter()
sim.save()

# Run the simulation. Arguments may be defined here or thorugh the swiftest.Simulation() method.
sim.run()
