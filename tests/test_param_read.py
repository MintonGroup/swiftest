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

class TestSwiftestUnits(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir=tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name
        
    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup() 

    def test_param_read(self):
        """
        Test that a user defined param_file name is read correctly
        """

        user_param_file = 'user_param.in'

        sim = swiftest.Simulation(simdir = self.simdir, param_file = user_param_file)

        self.assertEqual(str(sim.param_file), user_param_file)

        # new param_file name
        new_param_file = 'new_param.in'
        sim.set_parameter(param_file = new_param_file)

        self.assertEqual(str(sim.param_file), new_param_file)
    
    
    def test_param_read_restart(self):
        """
        Test that a param.restart.in file is read correctly
        """

        user_param_file = 'param.in'
        sim = swiftest.Simulation(simdir = self.simdir, integrator = 'symba', 
                                read_param = True, param_file = user_param_file,
                                dump_cadence = 1, tstep_out = 0.5, 
                                tstop = 1, dt = 0.01)
        sim.clean()

        sim.add_solar_system_body(['Sun', 'Earth', 'Mars'])

        sim.run()

        # restart run as a second check
        restarted_param_file = 'param.restart.in'
        sim_restart = swiftest.Simulation(simdir = self.simdir, 
                                        read_param = True, read_collisions=True, read_data = True,
                                        param_file = restarted_param_file)

        self.assertEqual(str(sim_restart.param_file), restarted_param_file)

if __name__ == '__main__':
    os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
    unittest.main()








