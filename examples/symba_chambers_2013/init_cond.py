#!/usr/bin/env python3
import swiftest
import numpy as np
from numpy.random import default_rng

# Initialize simulation object
sim = swiftest.Simulation()

# Set unit conversion factors
MU2KG = swiftest.MSun
TU2S = swiftest.YR2S
DU2M = swiftest.AU2M
GU = swiftest.GC / (DU2M**3 / (MU2KG * TU2S**2))
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
sim.param['ENC_OUT'] = ""

# Simulation parameters
sim.param['FRAGMENTATION'] = "YES"
sim.param['ROTATION'] = "YES"
sim.param['CHK_RMAX'] = 1000.0
sim.param['CHK_EJECT'] = 1000.0
sim.param['IN_FORM'] = 'EL'
sim.param['OUT_FORM'] = 'EL'

# Add central body
sim.add("Sun")
GMcb = sim.ds['GMass'].values[0]
sim.add("Jupiter")
sim.add("Saturn")

# Add bodies described in Chambers (2013) Sec. 2.1, with the uniform spatial distribution and two bodies sizes (big and small)
Nb = 14
Ns = 140
Mb = 2.8e-7
Ms = 2.8e-8
dens = 3000.0 / (MU2KG / DU2M**3)
Rb = (3 * Mb / (4 * np.pi * dens) )**(1.0 / 3.0)
Rs = (3 * Ms / (4 * np.pi * dens) )**(1.0 / 3.0)
sim.param['GMTINY'] = 1e-2 * GU * Ms
sim.param['MIN_GMFRAG'] = 1e-4 * GU * Ms

# Define the initial orbital elements of the big and small bodies
avalb = default_rng().uniform(0.3, 2.0, Nb)
avals = default_rng().uniform(0.3, 2.0, Ns)
evalb = default_rng().uniform(0.0, 0.01, Nb)
evals = default_rng().uniform(0.0, 0.01, Ns)
incvalb = default_rng().uniform(0.0, 0.5, Nb)
incvals = default_rng().uniform(0.0, 0.5, Ns)
capomvalb = default_rng().uniform(0.0, 360.0, Nb)
capomvals = default_rng().uniform(0.0, 360.0, Ns)
omegavalb = default_rng().uniform(0.0, 360.0, Nb)
omegavals = default_rng().uniform(0.0, 360.0, Ns)
capmvalb = default_rng().uniform(0.0, 360.0, Nb)
capmvals = default_rng().uniform(0.0, 360.0, Ns)
GMvalb = np.full(Nb, Mb * GU)
GMvals = np.full(Ns, Ms * GU)
Rvalb = np.full(Nb, Rb)
Rvals = np.full(Ns, Rs)
Rhb = avalb * (GMvalb / (3 * GMcb))**(1.0/3.0)
Rhs = avals * (GMvals / (3 * GMcb))**(1.0/3.0)

# Give the bodies unique ids
idb = np.arange(100, 100 + Nb)
ids = np.arange(100 + Nb, 100 + Nb + Ns)

# Populate the simulation object with the two types of bodies
sim.addp(idb, avalb, evalb, incvalb, capomvalb, omegavalb, capmvalb, GMpl=GMvalb, Rpl=Rvalb, Rhill=Rhb)
sim.addp(ids, avals, evals, incvals, capomvals, omegavals, capmvals, GMpl=GMvals, Rpl=Rvals, Rhill=Rhs)

# Save everything to a set of initial conditions files
sim.save("param.in")



