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
import warnings
warnings.simplefilter('error', RuntimeWarning)

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
       sim = swiftest.Simulation(simdir=self.simdir,compute_conservation_values=True, rotation=True, integrator="symba")

       # Add the modern planets and the Sun using the JPL Horizons Database.
       sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune","Pluto"])

       density  = 3000.0 * sim.KG2MU / sim.M2DU**3

       # Make a body with a periapsis inside the Sun's radius
       q = 0.01 * swiftest.RSun * sim.M2DU
       a = 0.1
       e = 1.0 - q / a
       M = 2e0 * swiftest.MEarth * sim.KG2MU
       R = (3 * M  / (4 * np.pi * density)) ** (1.0 / 3.0)
       rot = 4 * sim.init_cond.sel(name="Earth")['rot']
       sim.add_body(name="Sundiver", a=a, e=e, inc=0.0, capom=0.0, omega=0.0, capm=180.0, mass=M, radius=R, Ip=[0.4,0.4,0.4], rot=rot)
 
       return 
         
if __name__ == '__main__':
    os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
    unittest.main()