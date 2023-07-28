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
Tests that energy and momentum errors are within tolerances in a Swiftest simulation

Input
------

Output
------
None
"""

import swiftest
import unittest
import os
import numpy as np
from numpy.random import default_rng

rng = default_rng(seed=123)

major_bodies = ["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"]
param = {}


class TestSwiftest(unittest.TestCase):
   
   def test01_gen_ic(self):
      """
      Tests that Swiftest is able to successfully generate a set of initial conditions in a file without any exceptions being raised
      """
      print("\ntest_gen_ic: Test whether we can generate simulation initial conditions test")
      # Files that are expected to be generated:
      simdir = "simdata"
      file_list = [simdir, os.path.join(simdir,"param.in"), os.path.join(simdir,"init_cond.nc")]
      
      sim = swiftest.Simulation()
      sim.clean()

      # Add the modern planets and the Sun using the JPL Horizons Database.
      # Add the modern planets and the Sun using the JPL Horizons Database.
      sim.add_solar_system_body(major_bodies)
   
      # Display the run configuration parameters.
      param = sim.get_parameter(verbose=False)
      sim.save()
      
      for f in file_list:
         self.assertTrue(os.path.exists(f))
         
   def test02_read_ic(self):
      """
      Tests that Swiftest is able to read a set of pre-existing initial conditions files and that they contain the correct data
      """
      print("\ntest_read_ic: Test whether we can read back initial conditions files created by test_gen_ic")
      sim = swiftest.Simulation(read_param=True)
      # Check if all names in Dataset read in from file match the expected list of names
      self.assertTrue((major_bodies == sim.init_cond['name']).all(), msg="Name mismatch in Dataset")
      
      # Check to see if all parameter values read in from file match the expected parameters saved when generating the file
      self.assertTrue(all([v == param[k] for k,v in sim.param.items() if k in param]))
     
   def test03_integrators(self):
      """
      Tests that Swiftest is able to integrate a collection of massive bodies and test particles with all available integrators
      """ 
      print("\ntest_integrators: Tests that Swiftest is able to integrate a collection of massive bodies and test particles with all available integrators")
      sim = swiftest.Simulation(read_param=True)
      
      # Add 10 user-defined test particles.
      ntp = 10

      name_tp     = [f"TestParticle_{i:02}" for i in range(1,ntp+1)]
      a_tp        = rng.uniform(0.3, 1.5, ntp)
      e_tp        = rng.uniform(0.0, 0.2, ntp)
      inc_tp      = rng.uniform(0.0, 10, ntp)
      capom_tp    = rng.uniform(0.0, 360.0, ntp)
      omega_tp    = rng.uniform(0.0, 360.0, ntp)
      capm_tp     = rng.uniform(0.0, 360.0, ntp)

      integrators= ["whm","helio","rmvs","symba"]
      sim.add_body(name=name_tp, a=a_tp, e=e_tp, inc=inc_tp, capom=capom_tp, omega=omega_tp, capm=capm_tp)
      sim.set_parameter(tstart=0.0, tstop=0.02, dt=0.01, istep_out=1, dump_cadence=0)
      integrators= ["whm","helio","rmvs","symba"]
      for i in integrators:
         try:
            sim.run(integrator=i)
         except:
            self.fail(f"Failed with integrator {i}")
         
         
   def test04_conservation(self):
      """
      Tests that Swiftest conserves mass, energy, and momentum to within acceptable tolerances.
      """
      print("\ntest_conservation: Tests that Swiftest conserves mass, energy, and momentum to within acceptable tolerances.")
      sim = swiftest.Simulation(read_param=True) 
      
      # Add 5 user-defined semi-interacting massive bodies.
      npl         = 5
      density_pl  = 3000.0 / (sim.param['MU2KG'] / sim.param['DU2M'] ** 3)

      name_pl     = [f"SemiBody_{i:02}" for i in range(1,npl+1)]
      a_pl        = rng.uniform(0.3, 1.5, npl)
      e_pl        = rng.uniform(0.0, 0.2, npl)
      inc_pl      = rng.uniform(0.0, 10, npl)
      capom_pl    = rng.uniform(0.0, 360.0, npl)
      omega_pl    = rng.uniform(0.0, 360.0, npl)
      capm_pl     = rng.uniform(0.0, 360.0, npl)
      M_pl        = np.array([6e20, 8e20, 1e21, 3e21, 5e21]) * sim.KG2MU
      R_pl        = np.full(npl, (3 * M_pl/ (4 * np.pi * density_pl)) ** (1.0 / 3.0))
      Ip_pl       = np.full((npl,3),0.4,)
      rot_pl      = np.zeros((npl,3))
      mtiny       = 1.1 * np.max(M_pl)
      
      sim.add_body(name=name_pl, a=a_pl, e=e_pl, inc=inc_pl, capom=capom_pl, omega=omega_pl, capm=capm_pl, mass=M_pl, radius=R_pl,  Ip=Ip_pl, rot=rot_pl)
           
      sim.run(start=0.0, tstop=1.0e3, dt=0.01, istep_out=100, dump_cadence=0, compute_conservation_values=True, mtiny=mtiny, integrator="symba") 
      # Calculate the angular momentum error
      sim.data['L_tot'] = sim.data['L_orbit'] + sim.data['L_spin'] + sim.data['L_escape']
      sim.data['DL'] = sim.data['L_tot'] - sim.data['L_tot'].isel(time=0)
      sim.data['L_error'] = swiftest.tool.magnitude(sim.data,'DL') / swiftest.tool.magnitude(sim.data.isel(time=0), 'L_tot')
      L_final = sim.data['L_error'][-1].values

      # Calculate the energy error
      sim.data['E_error'] = (sim.data['TE'] - sim.data['TE'].isel(time=0)) / sim.data['TE'].isel(time=0)
      E_final = sim.data['E_error'][-1].values

      # Calculate the mass error
      sim.data['GMtot'] = sim.data['Gmass'].sum(dim='name',skipna=True) + sim.data['GMescape']
      sim.data['GM_error'] = (sim.data['GMtot'] - sim.data['GMtot'].isel(time=0)) / sim.data['GMtot'].isel(time=0)
      GM_final = sim.data['GM_error'][-1].values

      # Print the final errors
      print("Final Angular Momentum Error: ", L_final)
      print("Final Energy Error: ", E_final)
      print("Final Mass Error: ", GM_final)

      # Determine if the errors are within bounds
      L_limit = 1e-10
      E_limit = 1e-5
      GM_limit = 1e-14

      self.assertLess(L_final,L_limit, msg=f"Angular Momentum Error of {L_final} higher than threshold value of {L_limit}")
      self.assertLess(E_final,E_limit, msg=f"Energy Error of {E_final} higher than threshold value of {E_limit}")
      self.assertLess(GM_final,GM_limit, msg=f"Mass Error of {GM_final} higher than threshold value of {GM_limit}")
      
if __name__ == '__main__':
   unittest.main()

