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
Generates and runs two sets of Swiftest input files from initial conditions with the WHM integrator. All simulation 
outputs for the general relativity run are stored in the /gr subdirectory while all simulation outputs for the run 
without general reelativity are stored in the /nogr subdirectory.

Input
------
None.

Output
------
whm_gr_mercury_precession.png : Portable Network Graphic file depicting the precession of Mercury's perihelion over time
                                with data sourced from the JPL Horizons database, Swiftest run with general relativity, 
                                and Swiftest run without general relativity.
                                
Two subdirectories:
gr/
nogr/

Each subdirecotry contains:
collisions.log               : A NetCDF file containing the collision output.
data.nc                      : A NetCDF file containing the simulation output.
init_cond.nc                 : A NetCDF file containing the initial conditions for the simulation.
param.00...0.in              : A series of parameter input files containing the parameters for the simulation at every output stage.
param.in                     : An ASCII file containing the parameters for the simulation.
param.restart.in             : An ASCII file containing the parameters for the simulation at the last output. 
swiftest.log                 : An ASCII file containing the information on the status of the simulation as it runs.
"""

import swiftest
from astroquery.jplhorizons import Horizons
import datetime
import numpy as np
import matplotlib.pyplot as plt

# Initialize the simulation object as a variable. Define the directory in which the output will be placed.
sim_gr = swiftest.Simulation(simdir="gr")
sim_gr.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])

# Initialize the simulation object as a variable. Define the directory in which the output will be placed.
sim_nogr = swiftest.Simulation(simdir="nogr")
sim_nogr.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])

# Define a set of arguments that apply to both runs. For a list of possible arguments, see the User Manual.
run_args = {"tstop":1000.0, "dt":0.005, "tstep_out":10.0, "dump_cadence": 0,"integrator":"whm"}

# Run both simulations.
sim_gr.run(**run_args,general_relativity=True)
sim_nogr.run(**run_args,general_relativity=False)

# Get the start and end date of the simulation so we can compare with the real solar system.
start_date = sim_gr.ephemeris_date
tstop_d = sim_gr.param['TSTOP'] * sim_gr.param['TU2S'] / swiftest.JD2S

stop_date = (datetime.datetime.fromisoformat(start_date) + datetime.timedelta(days=tstop_d)).isoformat()

#Get the ephemerides of Mercury for the same timeframe as the simulation.
obj = Horizons(id='1', location='@sun',
               epochs={'start':start_date, 'stop':stop_date,
                       'step':'10y'})
el = obj.elements()
t = (el['datetime_jd']-el['datetime_jd'][0]) / 365.25
varpi_obs = el['w'] + el['Omega']

varpisim_gr= sim_gr.data['varpi'].sel(name="Mercury")
varpisim_nogr= sim_nogr.data['varpi'].sel(name="Mercury")
tsim = sim_gr.data['time']

dvarpi_gr = np.diff(varpisim_gr)  * 3600 * 100 / run_args['tstep_out']
dvarpi_nogr = np.diff(varpisim_nogr)  * 3600 * 100 / run_args['tstep_out']
dvarpi_obs = np.diff(varpi_obs) / np.diff(t) * 3600 * 100

# Plot of the data and save the output plot.
fig, ax = plt.subplots()

ax.plot(t, varpi_obs, label="JPL Horizons",linewidth=2.5)
ax.plot(tsim, varpisim_gr, label="Swiftest WHM GR",linewidth=1.5)
ax.plot(tsim, varpisim_nogr, label="Swiftest WHM No GR",linewidth=1.5)
ax.set_xlabel('Time (y)')
ax.set_ylabel('Mercury $\\varpi$ (deg)')
ax.legend()
plt.savefig("whm_gr_mercury_precession.png",dpi=300)

# Print the data to the terminal.
print('Mean precession rate for Mercury long. peri. (arcsec/100 y)')
print(f'JPL Horizons         : {np.mean(dvarpi_obs)}')
print(f'Swiftest No GR       : {np.mean(dvarpi_nogr)}')
print(f'Swiftest GR          : {np.mean(dvarpi_gr)}')
print(f'Obs - Swiftest GR    : {np.mean(dvarpi_obs - dvarpi_gr)}')
print(f'Obs - Swiftest No GR : {np.mean(dvarpi_obs - dvarpi_nogr)}')
