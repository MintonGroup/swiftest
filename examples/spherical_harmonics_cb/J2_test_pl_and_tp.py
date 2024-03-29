#!/usr/bin/env python3

"""
 Copyright 2024 - The Minton Group at Purdue University
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest. 
 If not, see: https://www.gnu.org/licenses. 
"""

"""
Generates and runs a set of Swiftest input files from initial conditions for the Spherical Harmonics features with the WHM integrator. 
"""

import swiftest
import numpy as np

seed = 123
rng = np.random.default_rng(seed=seed)


# Central Body Parameters (just an oblate sphere to test) 
cb_mass = 6.1e18 # kg
cb_a = 160 # km
cb_b = 160 # km
cb_c = 90 # km
cb_volume = 4.0 / 3 * np.pi * cb_a*cb_b*cb_c**3 # km^3
cb_density = cb_mass / cb_volume 
cb_T_rotation = 7.004 / 24.0 # converting from hours to julian days (TU)
cb_rot = [[0, 0, 360.0 / cb_T_rotation]] # degrees/d

# Add 1 user-defined test particle.
ntp = 1

name_tp     = ["TestParticle_01"]
a_tp        = 400
e_tp        = 0.05
inc_tp      = 10
capom_tp    = 0.0
omega_tp    = 0.0
capm_tp     = 0.0


# Add 1 user-defined massive particle
npl = 1
density_pl = cb_density

name_pl     = ["MassiveBody_01"]
a_pl        = 300.0
e_pl        = 0.03
inc_pl      = 0.001
capom_pl    = 90.0
omega_pl    = 90.0
capm_pl     = 90.0
R_pl        = 1.0
M_pl        = 4.0 / 3 * np.pi * R_pl**3 * density_pl
Ip_pl       = np.full((npl,3),0.4,)
rot_pl      = np.zeros((npl,3))
mtiny       = 0.1 * np.max(M_pl)


# Extract the spherical harmonics coefficients (c_lm) from axes measurements
#
# The user can pass an optional reference radius at which the coefficients are calculated. If not provided, SHTOOLS 
# calculates the reference radius. If lref_radius = True, the function returns the reference radius used. 
# We recommend setting passing and setting a reference radius. Coefficients are geodesy (4-pi) normalised.

c_lm, cb_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, c = cb_c, lmax = 6, lref_radius = True)

# extracting only the J2 terms
tmp20 = c_lm[0, 2, 0] # c_20 = -J2
c_lm = np.zeros(np.shape(c_lm))
c_lm[0, 2, 0] = tmp20

J2 = -tmp20 * np.sqrt(5) # unnormalised J2 term
j2rp2 = J2 * cb_radius**2

# set up swiftest simulation with relevant units (here they are km, days, and kg)
sim_shgrav = swiftest.Simulation(simdir="shgrav",DU2M = 1e3, TU = 'd', MU = 'kg')

sim_shgrav.clean()
# Use the shgrav version where you input a set of spherical harmonics coefficients
sim_shgrav.add_body(name = 'OblateBody', mass = cb_mass, rot = cb_rot, radius = cb_radius, c_lm = c_lm)
sim_shgrav.add_body(name=name_tp, a=a_tp, e=e_tp, inc=inc_tp, capom=capom_tp, omega=omega_tp, capm=capm_tp)
sim_shgrav.add_body(name=name_pl, a=a_pl, e=e_pl, inc=inc_pl, capom=capom_pl, omega=omega_pl, capm=capm_pl, mass=M_pl, radius=R_pl,  Ip=Ip_pl, rot=rot_pl)
sim_shgrav.run(tstart=0.0, tstop=10.0, dt=0.01, tstep_out=10.0, dump_cadence=0, mtiny=mtiny, integrator='symba')

# Use the original "oblate" version where you pass J2 (and/or J4)
sim_obl = swiftest.Simulation(simdir="obl", DU2M = 1e3, TU='d', MU='kg')
sim_obl.clean()
sim_obl.add_body(name = 'OblateBody', mass = cb_mass, rot = cb_rot, radius = cb_radius, j2rp2 = j2rp2)
sim_obl.add_body(name=name_tp, a=a_tp, e=e_tp, inc=inc_tp, capom=capom_tp, omega=omega_tp, capm=capm_tp)
sim_obl.add_body(name=name_pl, a=a_pl, e=e_pl, inc=inc_pl, capom=capom_pl, omega=omega_pl, capm=capm_pl, mass=M_pl, radius=R_pl,  Ip=Ip_pl, rot=rot_pl)
sim_obl.run(tstart=0.0, tstop=10.0, dt=0.01, tstep_out=10.0, dump_cadence=0, mtiny=mtiny, integrator='symba')

diff_vars = ['a','e','inc','capom','omega','capm','rh','vh']
ds_diff = sim_shgrav.data[diff_vars] - sim_obl.data[diff_vars]
ds_diff /= sim_obl.data[diff_vars]

print(ds_diff.isel(time=-1,name=-2))
print(ds_diff.isel(time=-1,name=-1))

