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

"""
Generates and runs a set of Swiftest input files from initial conditions with the SyMBA integrator. All simulation 
outputs are stored in the /simdata subdirectory.


"""

from numpy.random import default_rng

import swiftest

# Initialize the simulation object as a variable. Arguments may be defined here or through the sim.run() method.
sim = swiftest.Simulation()
sim.clean()

# Add the modern planets and the Sun using the JPL Horizons Database.
sim.add_solar_system_body(["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"])

# Add 1000 user-defined test particles.
ntp = 100

a_tp = default_rng().uniform(2.0, 4.0, ntp)
e_tp = default_rng().uniform(0.0, 0.3, ntp)
inc_tp = default_rng().uniform(0.0, 45.0, ntp)
capom_tp = default_rng().uniform(0.0, 360.0, ntp)
omega_tp = default_rng().uniform(0.0, 360.0, ntp)
capm_tp = default_rng().uniform(0.0, 360.0, ntp)

sim.add_body(a=a_tp, e=e_tp, inc=inc_tp, capom=capom_tp, omega=omega_tp, capm=capm_tp)

# Run the simulation. Arguments may be defined here or thorugh the swiftest.Simulation() method.
sim.set_parameter(tstart=0.0, tstop=1.0e2, dt=0.01, tstep_out=1.0, dump_cadence=0, integrator="rmvs", coarray=True)
# Display the run configuration parameters.
sim.get_parameter()
sim.save()
