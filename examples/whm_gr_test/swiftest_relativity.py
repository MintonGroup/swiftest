#!/usr/bin/env python
import swiftest
from astroquery.jplhorizons import Horizons
import datetime
import numpy as np
import matplotlib.pyplot as plt

sim_gr = swiftest.Simulation(simdir="gr")
sim_gr.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])

sim_nogr = swiftest.Simulation(simdir="nogr")
sim_nogr.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])

tstep_out = 10.0
sim_gr.run(tstop=1000.0, dt=0.005, tstep_out=tstep_out, integrator="whm",general_relativity=True)
sim_nogr.run(tstop=1000.0, dt=0.005, tstep_out=tstep_out, integrator="whm",general_relativity=False)

# Get the start and end date of the simulation so we can compare with the real solar system
start_date = sim_gr.ephemeris_date
tstop_d = sim_gr.param['TSTOP'] * sim_gr.param['TU2S'] / swiftest.JD2S

stop_date = (datetime.datetime.fromisoformat(start_date) + datetime.timedelta(days=tstop_d)).isoformat()

#Get the ephemerides of Mercury for the same timeframe as the simulation
obj = Horizons(id='1', location='@sun',
               epochs={'start':start_date, 'stop':stop_date,
                       'step':'10y'})
el = obj.elements()
t = (el['datetime_jd']-el['datetime_jd'][0]) / 365.25
varpi_obs = el['w'] + el['Omega']

# Compute the longitude of the periapsis
sim_gr.data['varpi'] = np.mod(sim_gr.data['omega'] + sim_gr.data['capom'],360)
sim_nogr.data['varpi'] = np.mod(sim_nogr.data['omega'] + sim_nogr.data['capom'],360)

varpisim_gr= sim_gr.data['varpi'].sel(name="Mercury")
varpisim_nogr= sim_nogr.data['varpi'].sel(name="Mercury")
tsim = sim_gr.data['time']

dvarpi_gr = np.diff(varpisim_gr)  * 3600 * 100 / tstep_out
dvarpi_nogr = np.diff(varpisim_nogr)  * 3600 * 100 / tstep_out
dvarpi_obs = np.diff(varpi_obs) / np.diff(t) * 3600 * 100


fig, ax = plt.subplots()

ax.plot(t, varpi_obs, label="JPL Horizons",linewidth=2.5)
ax.plot(tsim, varpisim_gr, label="Swiftest WHM GR",linewidth=1.5)
ax.plot(tsim, varpisim_nogr, label="Swiftest WHM No GR",linewidth=1.5)
ax.set_xlabel('Time (y)')
ax.set_ylabel('Mercury $\\varpi$ (deg)')
ax.legend()
plt.savefig("whm_gr_mercury_precession.png",dpi=300)

print('Mean precession rate for Mercury long. peri. (arcsec/100 y)')
print(f'JPL Horizons         : {np.mean(dvarpi_obs)}')
print(f'Swiftest No GR       : {np.mean(dvarpi_nogr)}')
print(f'Swiftest GR          : {np.mean(dvarpi_gr)}')
print(f'Obs - Swiftest GR    : {np.mean(dvarpi_obs - dvarpi_gr)}')
print(f'Obs - Swiftest No GR : {np.mean(dvarpi_obs - dvarpi_nogr)}')
