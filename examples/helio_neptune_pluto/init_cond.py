#!/usr/bin/env python3
import swiftest

sim = swiftest.Simulation()
sim.param['IN_TYPE'] = "NETCDF_DOUBLE"
sim.param['NC_IN'] = "init_cond.nc"
sim.param['IN_FORM'] = "XV"
sim.param['BIN_OUT'] = "bin.nc"
sim.param['OUT_FORM'] = "XVEL"
sim.param['OUT_TYPE'] = "NETCDF_DOUBLE"

sim.param['MU2KG'] = swiftest.MSun
sim.param['TU2S'] = swiftest.YR2S
sim.param['DU2M'] = swiftest.AU2M
sim.param['T0'] = 0.0
sim.param['DT'] = 5.0 
sim.param['TSTOP'] = 1.e5
sim.param['ISTEP_OUT']  = 100
sim.param['ISTEP_DUMP'] = 100
sim.param['CHK_QMIN_COORD'] = "HELIO"
sim.param['CHK_QMIN'] = swiftest.RSun / swiftest.AU2M
sim.param['CHK_QMIN_RANGE'] = f"{swiftest.RSun / swiftest.AU2M} 1000.0"
sim.param['CHK_RMIN'] = swiftest.RSun / swiftest.AU2M
sim.param['CHK_RMAX'] = 1000.0
sim.param['CHK_EJECT'] = 1000.0
sim.param['OUT_STAT'] = "REPLACE"
sim.param['RHILL_PRESENT'] = "NO"
sim.param['GR'] = 'NO'
sim.param['CHK_CLOSE'] = "NO"

bodyid = {
   "Sun": 0,
   "Neptune": 8,
   "Pluto" : 9
}

for name, id in bodyid.items():
   sim.add(name, idval=id, date="2027-04-30")
   
sim.save("param.in")


