"""
 Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest.
 If not, see: https://www.gnu.org/licenses.
"""
#!/usr/bin/env python3
"""
Generates a set of Swiftest input files from initial conditions, runs Swiftest, and generates fragmentation movies.
Returns
-------
Updates sim.data with the simulation data, produces fragmentation movies.
"""
import swiftest
import numpy as np
from numpy.random import default_rng

sim_disruption = swiftest.Simulation(simdir="disruption", tstart=0.0, tstop=1.0e-5, dt=1.0e-8, istep_out=1.0, fragmentation=True, minimum_fragment_gmass=1.0e-11, gmtiny=1.0e-11, output_file_format="XVEL", init_cond_format="XV")
sim_disruption.add_solar_system_body(["Sun"])
sim_disruption.add_body(name="Target", rh=[[1.0, -1.807993e-05, 0.0]], vh=[[-2.562596e-04, 6.280005, 0.0]], Gmass=1e-7, radius=7e-6, rhill=9e-4, Ip=[[0.4, 0.4, 0.4]], rot=[[0.0, 0.0, 0.0]])
sim_disruption.add_body(name="Projectile", rh=[[1.0, 1.807993e-05, 0.0]], vh=[[-2.562596e-04, -6.280005, 0.0]], Gmass=7e-10, radius=3.25e-6, rhill=4e-4, Ip=[[0.4, 0.4, 0.4]], rot=[[0.0, 0.0, 0.0]])
sim_disruption.get_parameter()
sim_disruption.run()

sim_hitandrun = swiftest.Simulation(simdir="hitandrun", tstart=0.0, tstop=1.0e-5, dt=1.0e-8, istep_out=1.0, fragmentation=True, minimum_fragment_gmass=1.0e-11, gmtiny=1.0e-11, output_file_format="XVEL", init_cond_format="XV")
sim_hitandrun.add_solar_system_body(["Sun"])
sim_hitandrun.add_body(name="Target", rh=[[1.0, -4.2e-05, 0.0]], vh=[[0.0, 6.28, 0.0]], Gmass=1e-7, radius=7e-6, rhill=9e-4, Ip=[[0.4, 0.4, 0.4]], rot=[[0.0, 0.0, 6.0e4]])
sim_hitandrun.add_body(name="Projectile", rh=[[1.0, 4.2e-05, 0.0]], vh=[[-1.5, -6.28, 0.0]], Gmass=7e-10, radius=3.25e-6, rhill=4e-4, Ip=[[0.4, 0.4, 0.4]], rot=[[0.0, 0.0, 1.0e5]])
sim_hitandrun.get_parameter()
sim_hitandrun.run()

sim_supercat = swiftest.Simulation(simdir="supercat", tstart=0.0, tstop=1.0e-5, dt=1.0e-8, istep_out=1.0, fragmentation=True, minimum_fragment_gmass=1.0e-11, gmtiny=1.0e-11, output_file_format="XVEL", init_cond_format="XV")
sim_supercat.add_solar_system_body(["Sun"])
sim_supercat.add_body(name="Target", rh=[[1.0, -4.2e-05, 0.0]], vh=[[0.0, 6.28, 0.0]], Gmass=1e-7, radius=7e-6, rhill=9e-4, Ip=[[0.4, 0.4, 0.4]], rot=[[0.0, 0.0, -6.0e4]])
sim_supercat.add_body(name="Projectile", rh=[[1.0, 4.2e-05, 0.0]], vh=[[1.0, -6.28, 0.0]], Gmass=1e-8, radius=3.25e-6, rhill=4e-4, Ip=[[0.4, 0.4, 0.4]], rot=[[0.0, 0.0, 1.0e5]])
sim_supercat.get_parameter()
sim_supercat.run()
