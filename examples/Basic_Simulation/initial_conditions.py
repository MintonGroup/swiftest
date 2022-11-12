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
Generates a set of Swiftest input files from initial conditions.

Returns
-------
param.in : ASCII text file
    Swiftest parameter input file.
pl.in    : ASCII text file
    Swiftest massive body input file.
tp.in    : ASCII text file
    Swiftest test particle input file.
cb.in    : ASCII text file
    Swiftest central body input file.
"""

import swiftest
import numpy as np
from numpy.random import default_rng

# Initialize the simulation object as a variable
sim = swiftest.Simulation(tstart=0.0, tstop=10.0, dt=0.005, tstep_out=1.0, fragmentation=True, minimum_fragment_mass = 2.5e-11, mtiny=2.5e-8)
sim.get_parameter()

# Add the modern planets and the Sun using the JPL Horizons Database
sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune","Pluto"])

# Add 5 user-defined massive bodies
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
Ip1_pl      = [0.4, 0.4, 0.4, 0.4, 0.4]
Ip2_pl      = [0.4, 0.4, 0.4, 0.4, 0.4]
Ip3_pl      = [0.4, 0.4, 0.4, 0.4, 0.4]
rotx_pl     = [0.0, 0.0, 0.0, 0.0, 0.0]
roty_pl     = [0.0, 0.0, 0.0, 0.0, 0.0]
rotz_pl     = [0.0, 0.0, 0.0, 0.0, 0.0]

sim.add_body(name_pl, a_pl, e_pl, inc_pl, capom_pl, omega_pl, capm_pl, GMpl=GM_pl, Rpl=R_pl, rhill=Rh_pl, Ip1=Ip1_pl, Ip2=Ip2_pl, Ip3=Ip3_pl, rotx=rotx_pl, roty=roty_pl, rotz=rotz_pl)

# Add 10 user-defined test particles
ntp = 10

name_tp     = ["TestParticle_01", "TestParticle_02", "TestParticle_03", "TestParticle_04", "TestParticle_05", "TestParticle_06", "TestParticle_07", "TestParticle_08", "TestParticle_09", "TestParticle_10"]
a_tp        = default_rng().uniform(0.3, 1.5, ntp)
e_tp        = default_rng().uniform(0.0, 0.3, ntp)
inc_tp      = default_rng().uniform(0.0, 90, ntp)
capom_tp    = default_rng().uniform(0.0, 360.0, ntp)
omega_tp    = default_rng().uniform(0.0, 360.0, ntp)
capm_tp     = default_rng().uniform(0.0, 360.0, ntp)

sim.add_body(name_tp, a_tp, e_tp, inc_tp, capom_tp, omega_tp, capm_tp)

# Save everything to a set of initial conditions files
sim.save('param.in')
