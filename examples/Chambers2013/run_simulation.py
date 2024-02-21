#!/usr/bin/env python3

"""
 Copyright 2023 - David Minton
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest. 
 If not, see: https://www.gnu.org/licenses. 
"""

"""
This will run the simulation from a set of initial conditions. The simulation parameters in this file are set to generate
a very short simulation for testing purposes. Edit the values passed to the run() function as necessary.

Input
------
simdata/param.in    : ASCII Swiftest parameter input file.

Output
------
Outputs are stored in the /simdata subdirectory.

"""
import swiftest
sim = swiftest.Simulation(read_param=True)

# Original run parameters
# sim.run(tstop=3e8, dt=6.0875/365.25, istep_out=60000, dump_cadence=10,integreator="symba")
# 
sim.run(tstop=10000.0, dt=6.0875/365.25, istep_out=10000, dump_cadence=0, integrator="symba")
