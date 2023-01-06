#!/usr/bin/env python3
import swiftest
import numpy as np
from numpy.random import default_rng

# Initialize simulation object
sim = swiftest.Simulation(compute_conservation_values=True,  rotation=True, init_cond_format="EL",collision_model="fraggle",encounter_save="none")

# Add bodies described in Chambers (2013) Sec. 2.1, with the uniform spatial distribution and two bodies sizes (big and small)
Nb = 14
Ns = 140
Mb = 2.8e-7 * 14 / Nb
Ms = 2.8e-8 * 140 / Ns
dens = 3000.0 / (sim.param['MU2KG'] / sim.param['DU2M']**3)
Rb = (3 * Mb / (4 * np.pi * dens) )**(1.0 / 3.0)
Rs = (3 * Ms / (4 * np.pi * dens) )**(1.0 / 3.0)
mtiny = 1e-2 * Ms
mininum_fragment_mass = 1e-4 * Ms
rng = default_rng(seed=170834)

# Define the initial orbital elements of the big and small bodies
avalb = rng.uniform(0.3, 2.0, Nb)
avals = rng.uniform(0.3, 2.0, Ns)
evalb = rng.uniform(0.0, 0.01, Nb)
evals = rng.uniform(0.0, 0.01, Ns)
incvalb = rng.uniform(0.0, 0.005 * 180 / np.pi, Nb)
incvals = rng.uniform(0.0, 0.005 * 180 / np.pi, Ns)
capomvalb = rng.uniform(0.0, 360.0, Nb)
capomvals = rng.uniform(0.0, 360.0, Ns)
omegavalb = rng.uniform(0.0, 360.0, Nb)
omegavals = rng.uniform(0.0, 360.0, Ns)
capmvalb = rng.uniform(0.0, 360.0, Nb)
capmvals = rng.uniform(0.0, 360.0, Ns)
Ipvalb = np.full((Nb,3), 0.4)
Ipvals = np.full((Ns,3), 0.4)
rotvalb = np.zeros_like(Ipvalb)
rotvals = np.zeros_like(Ipvals)
GMvalb = np.full(Nb, Mb * sim.GU)
GMvals = np.full(Ns, Ms * sim.GU)
Rvalb = np.full(Nb, Rb)
Rvals = np.full(Ns, Rs)

# Give the bodies unique names
nameb = [f"Big{i:03}" for i in range(Nb)]
names = [f"Small{i:03}" for i in range(Ns)]

# Add the modern planets and the Sun using the JPL Horizons Database.
sim.add_solar_system_body(["Sun","Jupiter","Saturn","Uranus","Neptune"])
sim.add_body(name=nameb, a=avalb, e=evalb, inc=incvalb, capom=capomvalb, omega=omegavalb, capm=capmvalb, Gmass=GMvalb, radius=Rvalb, rot=rotvalb, Ip=Ipvalb)
sim.add_body(name=names, a=avals, e=evals, inc=incvals, capom=capomvals, omega=omegavals, capm=capmvals, Gmass=GMvals, radius=Rvals, rot=rotvals, Ip=Ipvals)
sim.set_parameter(mtiny=mtiny, minimum_fragment_mass=mininum_fragment_mass)

sim.run(tstop=3e8, dt=6.0875/365.25, istep_out=60000, dump_cadence=1)

