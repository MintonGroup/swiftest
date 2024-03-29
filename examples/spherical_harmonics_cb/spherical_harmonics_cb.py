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
Generates and runs a set of Swiftest input files from initial conditions for the Spherical Harmonics features with the 
SyMBA integrator. Using Chariklo as the example body with axes measurements taken from Leiva, et al (2017) (Jacobi 
Ellipsoid model). All simulation outputs are stored in the /simdata subdirectory.

"""

import swiftest
import numpy as np

seed = 123
rng = np.random.default_rng(seed=seed)

# set up swiftest simulation with relevant units (here they are km, days, and kg)
sim = swiftest.Simulation(DU2M = 1e3, TU = 'd', MU = 'kg')

# Central Body Parameters (Chariklo parameters from Leiva, et al (2017) (Jacobi Ellipsoid model))
cb_mass = 6.1e18 # kg
cb_radius = 123 # km
cb_a = 157 # km
cb_b = 139 # km
cb_c = 86 # km
cb_volume = 4.0 / 3 * np.pi * cb_radius**3 # km^3
cb_density = cb_mass / cb_volume 
cb_T_rotation = 7.004 / 24.0 # converting from hours to julian days (TU)
cb_rot = [[0, 0, 360.0 / cb_T_rotation]] # degrees/d

# Extract the spherical harmonics coefficients (c_lm) from axes measurements
#
# The user can pass an optional reference radius at which the coefficients are calculated. If not provided, SHTOOLS 
# calculates the reference radius. If lref_radius = True, the function returns the reference radius used. 
# We recommend setting passing and setting a reference radius. Coefficients are geodesy (4-pi) normalised.

c_lm, cb_radius = swiftest.clm_from_ellipsoid(mass = cb_mass, density = cb_density, a = cb_a, b = cb_b, c = cb_c, lmax = 6, lref_radius = True, ref_radius = cb_radius)

# Add the central body
# The user can pass the c_lm coefficients directly to the add_body method if they do not wish to use the clm_from_ellipsoid method.
sim.add_body(name = 'Chariklo', mass = cb_mass, rot = cb_rot, radius = cb_radius, c_lm = c_lm)

# Add user-defined massive bodies
npl         = 5
density_pl  = cb_density

name_pl     = ["SemiBody_01", "SemiBody_02", "SemiBody_03", "SemiBody_04", "SemiBody_05"]
a_pl        = rng.uniform(250, 400, npl)
e_pl        = rng.uniform(0.0, 0.05, npl)
inc_pl      = rng.uniform(0.0, 10, npl)
capom_pl    = rng.uniform(0.0, 360.0, npl)
omega_pl    = rng.uniform(0.0, 360.0, npl)
capm_pl     = rng.uniform(0.0, 360.0, npl)
R_pl        = np.array([0.5, 1.0, 1.2, 0.75, 0.8])
M_pl        = 4.0 / 3 * np.pi * R_pl**3 * density_pl
Ip_pl       = np.full((npl,3),0.4,)
rot_pl      = np.zeros((npl,3))
mtiny       = 1.1 * np.max(M_pl)

sim.add_body(name=name_pl, a=a_pl, e=e_pl, inc=inc_pl, capom=capom_pl, omega=omega_pl, capm=capm_pl, mass=M_pl, radius=R_pl,  Ip=Ip_pl, rot=rot_pl)

# Add 10 user-defined test particles.
ntp = 10

name_tp     = ["TestParticle_01", "TestParticle_02", "TestParticle_03", "TestParticle_04", "TestParticle_05", "TestParticle_06", "TestParticle_07", "TestParticle_08", "TestParticle_09", "TestParticle_10"]
a_tp        = rng.uniform(250, 400, ntp)
e_tp        = rng.uniform(0.0, 0.05, ntp)
inc_tp      = rng.uniform(0.0, 10, ntp)
capom_tp    = rng.uniform(0.0, 360.0, ntp)
omega_tp    = rng.uniform(0.0, 360.0, ntp)
capm_tp     = rng.uniform(0.0, 360.0, ntp)

sim.add_body(name=name_tp, a=a_tp, e=e_tp, inc=inc_tp, capom=capom_tp, omega=omega_tp, capm=capm_tp)
sim.set_parameter(tstart=0.0, tstop=10.0, dt=0.01, istep_out=10, dump_cadence=0, compute_conservation_values=True, mtiny=mtiny)

# Display the run configuration parameters.
sim.get_parameter()

# Run the simulation. Arguments may be defined here or thorugh the swiftest.Simulation() method.
sim.run()
