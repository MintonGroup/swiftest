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
import astropy.constants as const

class TestUnits(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir=tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name
        
    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup() 
        
    def test_xv2el2xv_collisions(self):
        '''
        Check that the xv2el and el2xv converts correctly for sim.collisions
        '''

        print('\ntest_xv2el2xv_conversion_dims: Check that the xv2el and el2xv converts with the right dimensions for sim.collisions')

        sim = swiftest.Simulation(simdir=self.simdir, tstop = 0.5, dt = 0.1, init_cond_format = 'XV')

        sim.add_solar_system_body(['Sun', 'Mars'])

        # add massive impactors that hit the planet
        sim.add_body(name = ['pl1'], mass = 1e22 / sim.MU2KG, radius = 1e6 / sim.DU2M,
                        rh = [sim.data.sel(name = 'Mars', time = 0)['rh'].values], vh = [sim.data.sel(name = 'Mars', time = 0)['vh'].values])

        sim.run()
        GMcb = sim.data.sel(name = 'Mars', time = 0)['Gmass'].values

        try:
            sim.collisions.xv2el(GMcb)
        except:
            self.fail('xv2el failed with sim.collisions')
        
        try:
            sim.collisions.el2xv(GMcb)
        except:
            self.fail('el2xv failed with sim.collisions after xv2el')
        

        return 
         
if __name__ == '__main__':
    os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
    unittest.main()