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
import swiftest
import unittest
import os
import tempfile
import numpy as np
class TestCollisions(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir=tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name
        
    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup() 
        
    def test_solar_impact(self):
        '''
        Tests that impacts into the central body work correctly
        '''
        sim = swiftest.Simulation(simdir=self.simdir,compute_conservation_values=True, rotation=True)

        # Add the modern planets and the Sun using the JPL Horizons Database.
        sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune","Pluto"])

        density  = 3000.0 * sim.KG2MU / sim.M2DU**3

        # Make a massive body with a periapsis inside the Sun's radius
        q = 0.01 * swiftest.RSun * sim.M2DU
        a = 0.1
        e = 1.0 - q / a
        M = 2e0 * swiftest.MEarth * sim.KG2MU
        R = (3 * M  / (4 * np.pi * density)) ** (1.0 / 3.0)
        rot = 4 * sim.init_cond.sel(name="Earth")['rot']
        sim.add_body(name="Sundiver", a=a, e=e, inc=0.0, capom=0.0, omega=0.0, capm=180.0, mass=M, radius=R, Ip=[0.4,0.4,0.4], rot=rot)
        
        sim.run(tstart=0.0, tstop=5e-2, dt=0.0001, istep_out=1, dump_cadence=0, integrator="symba")
        
        # Check that the collision actually happened
        self.assertEqual(sim.collisions.collision_id.size,1) 
        
        # Check that angular momentum is conserved
        ds=sim.collisions.sel(collision_id=1)
        ds['Ltot']=ds.L_orbit+ds.L_spin
        ds['Ltot_mag']=ds.Ltot.magnitude()
        dLtot=ds.Ltot_mag.diff('stage').values[0]
        self.assertAlmostEqual(dLtot,0,places=8)
        
        # Check that energy was lost
        dEtot=ds.TE.diff('stage').values[0]
        self.assertLess(dEtot,0)
        
        
        # Now run the same test but with a massless body using both the RMVS and Symba integrators
        sim = swiftest.Simulation(simdir=self.simdir,compute_conservation_values=False, integrator="symba")
        sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune","Pluto"])
        sim.add_body(name="Sundiver", a=a, e=e, inc=0.0, capom=0.0, omega=0.0, capm=180.0)
        sim.run(tstart=0.0, tstop=5e-2, dt=0.0001, istep_out=1, dump_cadence=0)
        
        

        return 
         
if __name__ == '__main__':
    os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
    unittest.main()