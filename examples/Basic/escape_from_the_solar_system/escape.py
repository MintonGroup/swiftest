#!/usr/bin/env python3
"""
Copyright 2025 - David Minton.

This file is part of Swiftest.
Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Swiftest.
If not, see: https://www.gnu.org/licenses.
"""

import numpy as np

import swiftest

# Initialize the simulation object as a variable. Arguments may be defined here or through the sim.run() method.
sim = swiftest.Simulation(compute_conservation_values=True, rotation=True, integrator="symba")

# Add the modern planets and the Sun using the JPL Horizons Database.
sim.add_solar_system_body(["Sun", "Jupiter", "Saturn", "Uranus", "Neptune"])

density = 3000.0 * sim.KG2MU / sim.M2DU**3

# Make a hyperbolic body
q = 1.0
a = 0.1
e = 1.1
M = 1e-4 * swiftest.MEarth * sim.KG2MU
R = (3 * M / (4 * np.pi * density)) ** (1.0 / 3.0)
rot = 2 * sim.init_cond.sel(name="Saturn")["rot"]
sim.add_body(name="Escapee", a=a, e=e, inc=0.0, capom=0.0, omega=0.0, capm=0.0, mass=M, radius=R, Ip=[0.4, 0.4, 0.4], rot=rot)
sim.set_parameter(tstart=0.0, tstop=1000.0, dt=0.05, istep_out=200, dump_cadence=0, mtiny=2 * M)
sim.get_parameter()

# Run the simulation. Arguments may be defined here or thorugh the swiftest.Simulation() method.
sim.run()
