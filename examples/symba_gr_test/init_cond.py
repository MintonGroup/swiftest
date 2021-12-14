#!/usr/bin/env python3
import swiftest
from numpy.random import default_rng
import numpy as np

sim = swiftest.Simulation()
sim.param['PL_IN'] = "pl.swiftest.in"
sim.param['TP_IN'] = "tp.swiftest.in"
sim.param['CB_IN'] = "cb.swiftest.in"
sim.param['BIN_OUT'] = "bin.swiftest.nc"

sim.param['MU2KG'] = swiftest.MSun
sim.param['TU2S'] = swiftest.YR2S
sim.param['DU2M'] = swiftest.AU2M
sim.param['T0'] = 0.0
sim.param['DT'] = 0.25 * swiftest.JD2S / swiftest.YR2S
sim.param['TSTOP'] = 1000.0
sim.param['ISTEP_OUT']  = 1461
sim.param['ISTEP_DUMP'] = 1461
sim.param['CHK_QMIN_COORD'] = "HELIO"
sim.param['CHK_QMIN'] = swiftest.RSun / swiftest.AU2M
sim.param['CHK_QMIN_RANGE'] = f"{swiftest.RSun / swiftest.AU2M} 1000.0"
sim.param['CHK_RMIN'] = swiftest.RSun / swiftest.AU2M
sim.param['CHK_RMAX'] = 1000.0
sim.param['CHK_EJECT'] = 1000.0
sim.param['OUT_STAT'] = "UNKNOWN"
sim.param['IN_FORM'] = "EL"
sim.param['OUT_FORM'] = "XVEL"
sim.param['OUT_TYPE'] = "NETCDF_DOUBLE"
sim.param['RHILL_PRESENT'] = "YES"
sim.param['GR'] = 'YES'
sim.param['GMTINY'] = '1e-12'

bodyid = {
   "Sun": 0,
   "Mercury": 1,
   "Venus": 2,
   "Earth": 3,
   "Mars": 4,
   "Jupiter": 5,
   "Saturn": 6,
   "Uranus": 7,
   "Neptune": 8,
}

for name, id in bodyid.items():
   sim.add(name, idval=id, date="2027-04-30")

Me_a = sim.ds.sel(id=1)['a'].values
Me_e = sim.ds.sel(id=1)['e'].values
Me_i = sim.ds.sel(id=1)['inc'].values

capom_pl = sim.ds.sel(id=1)['capom'].values
omega_pl = sim.ds.sel(id=1)['omega'].values
capm_pl = sim.ds.sel(id=1)['capm'].values + 1.0

capom_tp = sim.ds.sel(id=1)['capom'].values
omega_tp = sim.ds.sel(id=1)['omega'].values
capm_tp = sim.ds.sel(id=1)['capm'].values - 1.0

GMcb = sim.ds.isel(id=0)['Gmass'].values
GU = swiftest.GC / (sim.param['DU2M']**3 / (sim.param['MU2KG'] * sim.param['TU2S']**2))
dens = 3000.0 / (sim.param['MU2KG'] / sim.param['DU2M']**3) # Assume a bulk density of 3 g/cm^3
GM_pl = 1e-11
M_pl = GM_pl / GU
R_pl = (3 * M_pl / (4 * np.pi * dens))**(1.0 / 3.0)
Rh_pl = Me_a * (GM_pl / (3 * GMcb))**(1.0/3.0)

sim.addp([9], ['Planetesimal'], Me_a, Me_e, Me_i, capom_pl, omega_pl, capm_pl, GMpl=GM_pl, Rpl=R_pl, rhill=Rh_pl)
sim.addp([10], ['TestParticle'], Me_a, Me_e, Me_i, capom_tp, omega_tp, capm_tp)

sim.save("param.swiftest.in")


