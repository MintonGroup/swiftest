#!/usr/bin/env python3
import swiftest
import numpy as np
import os

# Initialize simulation object
sim = swiftest.Simulation()

# Set unit conversion factors
MU2KG = swiftest.MSun
TU2S = swiftest.YR2S
DU2M = swiftest.AU2M
sim.param['MU2KG'] = MU2KG
sim.param['TU2S'] = TU2S
sim.param['DU2M'] = DU2M

# Simulation time parameters
sim.param['T0'] = 0.0
sim.param['TSTOP'] = 300e6
sim.param['DT'] = 6 * swiftest.JD2S / sim.param['TU2S']
t_print = 1000.0
iout = int(np.ceil(t_print / sim.param['DT']))
sim.param['ISTEP_OUT']  = iout
sim.param['ISTEP_DUMP'] = iout

# Optional output file names
sim.param['PARTICLE_OUT'] = "particle.dat"
sim.param['ENERGY'] = "YES"
sim.param['ENERGY_OUT'] = "energy.dat"
sim.param['PL_IN'] = "pl_chambers_2013.in"
sim.param['CB_IN'] = "sun_MsunAUYR.in"

# Simulation parameters
sim.param['FRAGMENTATION'] = "YES"
sim.param['ROTATION'] = "YES"
sim.param['CHK_RMAX'] = 1000.0
sim.param['CHK_EJECT'] = 1000.0

# Add central body
sim.add("Sun")
sim.add("Earth")

# Add bodies described in Chambers (2013) Sec. 2.1, with the uniform spatial distribution and two bodies sizes (big and small)
Nb = 14
Ns = 140
Mb = 2.8e-7
Ms = 2.8e-8
dens = 3000.0 / (MU2KG / DU2M**3)
Rb = (3 * Mb / (4 * np.pi * dens) )**(1.0 / 3.0)
Rs = (3 * Ms / (4 * np.pi * dens) )**(1.0 / 3.0)


sim.save("param.in")



