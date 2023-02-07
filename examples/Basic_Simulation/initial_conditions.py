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

#!/usr/bin/env python3
"""
Generates and runs a set of Swiftest input files from initial conditions with the SyMBA integrator. All simulation 
outputs are stored in the /simdata subdirectory.

Input
------
None.

Output
------
data.nc         : A NetCDF file containing the simulation output.
dump_bin1.nc   : A NetCDF file containing the necessary inputs to restart a simulation from t!=0.
dump_bin2.nc   : A NetCDF file containing the necessary inputs to restart a simulation from t!=0.
dump_param1.in : An ASCII file containing the necessary parameters to restart a simulation.
dump_param2.in : An ASCII file containing the necessary parameters to restart a simulation.
fraggle.log    : An ASCII file containing the information of any collisional events that occured.
init_cond.nc   : A NetCDF file containing the initial conditions for the simulation.
param.in       : An ASCII file containing the parameters for the simulation.
swiftest.log   : An ASCII file containing the information on the status of the simulation as it runs.
"""

import swiftest
import numpy as np
from numpy.random import default_rng

# Initialize the simulation object as a variable. Arguments may be defined here or through the sim.run() method.
#sim = swiftest.Simulation(fragmentation=True, minimum_fragment_mass = 2.5e-11, mtiny=2.5e-8)
sim = swiftest.Simulation()

# Add the modern planets and the Sun using the JPL Horizons Database.
sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune","Pluto"])

# Add 5 user-defined massive bodies.
npl         = 5
density_pl  = 3000.0 / (sim.param['MU2KG'] / sim.param['DU2M'] ** 3)

name_pl     = ["MassiveBody_01", "MassiveBody_02", "MassiveBody_03", "MassiveBody_04", "MassiveBody_05"]
a_pl        = default_rng().uniform(0.3, 1.5, npl)
e_pl        = default_rng().uniform(0.0, 0.3, npl)
inc_pl      = default_rng().uniform(0.0, 90, npl)
capom_pl    = default_rng().uniform(0.0, 360.0, npl)
omega_pl    = default_rng().uniform(0.0, 360.0, npl)
capm_pl     = default_rng().uniform(0.0, 360.0, npl)
GM_pl       = (np.array([6e23, 8e23, 1e24, 3e24, 5e24]) / sim.param['MU2KG']) * sim.GU
R_pl        = np.full(npl, (3 * (GM_pl / sim.GU) / (4 * np.pi * density_pl)) ** (1.0 / 3.0))
Rh_pl       = a_pl * ((GM_pl) / (3 * sim.GU)) ** (1.0 / 3.0)
Ip_pl      = np.full((npl,3),0.4,)
rot_pl     = np.zeros((npl,3))

sim.add_body(name=name_pl, a=a_pl, e=e_pl, inc=inc_pl, capom=capom_pl, omega=omega_pl, capm=capm_pl, Gmass=GM_pl, radius=R_pl, rhill=Rh_pl, Ip=Ip_pl, rot=rot_pl)

# Add 10 user-defined test particles.
ntp = 10

name_tp     = ["TestParticle_01", "TestParticle_02", "TestParticle_03", "TestParticle_04", "TestParticle_05", "TestParticle_06", "TestParticle_07", "TestParticle_08", "TestParticle_09", "TestParticle_10"]
a_tp        = default_rng().uniform(0.3, 1.5, ntp)
e_tp        = default_rng().uniform(0.0, 0.3, ntp)
inc_tp      = default_rng().uniform(0.0, 90, ntp)
capom_tp    = default_rng().uniform(0.0, 360.0, ntp)
omega_tp    = default_rng().uniform(0.0, 360.0, ntp)
capm_tp     = default_rng().uniform(0.0, 360.0, ntp)

sim.add_body(name=name_tp, a=a_tp, e=e_tp, inc=inc_tp, capom=capom_tp, omega=omega_tp, capm=capm_tp)
# Display the run configuration parameters.
sim.write_param()
sim.get_parameter()

# Run the simulation. Arguments may be defined here or thorugh the swiftest.Simulation() method.
sim.run(tstart=0.0, tstop=1.0e3, dt=0.01, istep_out=100, dump_cadence=10)
