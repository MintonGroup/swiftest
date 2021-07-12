#!/usr/bin/env python3
import swiftest

sim = swiftest.Simulation()

sim.param['MU2KG'] = swiftest.MSun
sim.param['TU2S'] = swiftest.YR2S
sim.param['DU2M'] = swiftest.AU2M
sim.param['T0'] = 0.0
sim.param['DT'] = 0.25 * swiftest.JD2S / swiftest.YR2S
sim.param['TSTOP'] = 100.0
sim.param['ISTEP_OUT']  = 1461
sim.param['ISTEP_DUMP'] = 1461
sim.param['CHK_QMIN_COORD'] = "HELIO"
sim.param['CHK_QMIN'] = swiftest.RSun / swiftest.AU2M
sim.param['CHK_QMIN_RANGE'] = f"{swiftest.RSun / swiftest.AU2M} 1000.0"
sim.param['CHK_RMIN'] = swiftest.RSun / swiftest.AU2M
sim.param['CHK_RMAX'] = 1000.0
sim.param['CHK_EJECT'] = 1000.0
sim.param['OUT_FORM'] = "EL"
sim.param['OUT_STAT'] = "UNKNOWN"
sim.param['GR'] = 'YES'

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
   sim.add(name, idval=id)
   
sim.param['PL_IN'] = "pl.swiftest.in"
sim.param['TP_IN'] = "tp.swiftest.in"
sim.param['CB_IN'] = "cb.swiftest.in"
sim.param['BIN_OUT'] = "bin.swiftest.dat"
sim.param['ENC_OUT'] = "enc.swiftest.dat"
sim.save("param.swiftest.in")
sim.param['PL_IN'] = "pl.swifter.in"
sim.param['TP_IN'] = "tp.swifter.in"
sim.param['BIN_OUT'] = "bin.swifter.dat"
sim.param['ENC_OUT'] = "enc.swifter.dat"
sim.save("param.swifter.in", codename="Swifter")


