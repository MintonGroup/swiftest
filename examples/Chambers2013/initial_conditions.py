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
Generates a set of Swiftest input files from initial conditions with the SyMBA integrator. All simulation 
outputs are stored in the /simdata subdirectory.

Input
------
None.

Output
------
inital_conditions.png : A .png file depicting the simulation initial configuration.
init_cond.nc          : A NetCDF file containing the initial conditions for the simulation.
param.in              : An ASCII file containing the parameters for the simulation.
""" 

import swiftest
import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt

# Initialize simulation object
sim = swiftest.Simulation(compute_conservation_values=True, rotation=True, init_cond_format="EL",collision_model="fraggle",encounter_save="none")

# Add bodies described in Chambers (2013) Sec. 2.1, with the uniform spatial distribution and two bodies sizes (big and small)
Nb = 14
Ns = 140
Mb = 2.8e-7 * 14 / Nb
Ms = 2.8e-8 * 140 / Ns
dens = 3000.0 * sim.KG2MU / sim.M2DU**3
rot = 1e-6 / sim.TU2S # Use a small but non-zero value for the initial rotation state to prevent divide-by-zero errors in analysis

mtiny = 1e-2 * Ms
minimum_fragment_mass = 1e-5 * Ms
nfrag_reduction = 30.0
rng = default_rng(seed=3031179)

runname = "Chambers (2013)"
def a_profile(n_bodies):
    """
    Generates the profile described in Sec. 2.1 of Chambers:  
    
    *In all cases, the surface density R = 8 g/cm2 at 1 AU, varying as a**(-3/2), where a is the orbital semi-major axis. 
    The region with a < 0.7 AU deviates from this law, declining linearly with decreasing distance until R = 0 at 0.3 AU. 
    The outer edge of the disk is 2 AU in all cases.
    """
    def sample(r_inner, r_break, r_outer, slope1, slope2, size):
        """ 
        Define the functions to pull random semi-major axes from a distribution using a rejection sampling method
        This defines a 2-slope model with a break at r_break
        Based on (https://stackoverflow.com/questions/66874819/random-numbers-with-user-defined-continuous-probability-distribution)
        """
        y=np.ones([size])
        pdf=np.zeros([size])
        x=np.empty_like(y)
        while np.any(y>pdf): 
            x = np.where(y>pdf,rng.uniform(r_inner, r_outer, size=size),x)
          
            # The proportionality factor A ensures that the PDF approaches the same value on either side of the break point
            # Assumes the break point is the max of the PDF
            A=np.where(x < r_break,1.0,r_break**(slope1-slope2))
            slope=np.where(x < r_break,slope1+1,slope2+1)
            y = np.where(y>pdf,rng.uniform(0, A*r_break**slope, size=size),y)
            pdf = np.where(y>pdf,A*x**(slope),pdf)
        return x
    
    a_inner = 0.3
    a_break = 0.7
    a_outer = 2.0
    slope1 = 1.0
    slope2 = -1.5
   
    a_vals = sample(a_inner, a_break, a_outer, slope1, slope2, n_bodies)
    return a_vals

# Define the initial orbital elements of the big and small bodies
avalb = a_profile(Nb)
avals = a_profile(Ns)

esigma = 0.01
isigma = np.rad2deg(0.5 * esigma)
evalb = rng.rayleigh(scale=esigma, size=Nb)
evals = rng.rayleigh(scale=esigma, size=Ns)
incvalb = rng.rayleigh(scale=isigma, size=Nb)
incvals = rng.rayleigh(scale=isigma, size=Ns)

capomvalb = rng.uniform(0.0, 360.0, Nb)
capomvals = rng.uniform(0.0, 360.0, Ns)
omegavalb = rng.uniform(0.0, 360.0, Nb)
omegavals = rng.uniform(0.0, 360.0, Ns)
capmvalb = rng.uniform(0.0, 360.0, Nb)
capmvals = rng.uniform(0.0, 360.0, Ns)
Ipvalb = np.full((Nb,3), 0.4)
Ipvals = np.full((Ns,3), 0.4)
rotvalb = np.full_like(Ipvalb,rot)
rotvals = np.full_like(Ipvals,rot)

noise_digits = 4 # Approximately the number of digits of precision to vary the mass values to avoid duplicate masses
epsilon = np.finfo(float).eps
Mnoiseb = 1.0 + 10**noise_digits * rng.uniform(-epsilon,epsilon, Nb)
Mnoises = 1.0 + 10**noise_digits * rng.uniform(-epsilon,epsilon, Ns)

Mvalb = Mb * Mnoiseb
Mvals = Ms * Mnoises

Rvalb = (3 * Mvalb / (4 * np.pi * dens) )**(1.0 / 3.0)
Rvals = (3 * Mvals / (4 * np.pi * dens) )**(1.0 / 3.0)

# Give the bodies unique names
nameb = [f"Big{i:03}" for i in range(Nb)]
names = [f"Small{i:03}" for i in range(Ns)]

# Add the modern planets and the Sun using the JPL Horizons Database.
sim.add_solar_system_body(["Sun","Jupiter","Saturn","Uranus","Neptune"])
sim.add_body(name=nameb, a=avalb, e=evalb, inc=incvalb, capom=capomvalb, omega=omegavalb, capm=capmvalb, mass=Mvalb, radius=Rvalb, rot=rotvalb, Ip=Ipvalb)
sim.add_body(name=names, a=avals, e=evals, inc=incvals, capom=capomvals, omega=omegavals, capm=capmvals, mass=Mvals, radius=Rvals, rot=rotvals, Ip=Ipvals)
sim.set_parameter(mtiny=mtiny, minimum_fragment_mass=minimum_fragment_mass, nfrag_reduction=nfrag_reduction)

sim.set_parameter(tstop=3e8, dt=6.0875/365.25, istep_out=60000, dump_cadence=10)
sim.save()

# Plot the initial conditions
fig = plt.figure(figsize=(8.5, 11))
ax1 = plt.subplot(2, 1, 1)
ax2 = plt.subplot(2, 1, 2)
fig.suptitle(runname)

ic = sim.init_cond.isel(time=0)
radius = ic['radius'].values
markersize = radius * 4e5
markercolor = np.where(radius < 2.0e-5, np.full_like(markersize, "black",dtype=object), np.full_like(markersize,"red",dtype=object))
a = ic['a'].values
e = ic['e'].values
inc = ic['inc'].values

ax1.scatter(a, e, s=markersize, color=markercolor, label=None)
ax1.set_xlabel("Semimajor Axis (AU)")
ax1.set_ylabel("Eccentricity")
ax1.set_xlim([0.0, 2.5])
ax1.set_ylim([0.0,4*esigma])

ax2.scatter(a, inc, s=markersize, color=markercolor, label=None)
ax2.set_xlabel("Semimajor Axis (AU)")
ax2.set_ylabel("Inclination (deg)")
ax2.set_xlim([0.0, 2.5])
ax2.set_ylim([0.0,8*isigma])

ax1.scatter(-1, -1, label='planetesimal', color='black')
ax1.scatter(-1, -1, label='embryo', color='red') 
ax1.legend(loc='upper right', frameon=False)

plt.tight_layout()
plt.show()
fig.savefig('initial_conditions.png',dpi=300)

