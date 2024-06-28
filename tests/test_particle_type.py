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

class TestParticleType(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir=tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name
        
    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup() 
        
    def test_correct_particle_types(self):
        '''
        Check that a variety of particles are set up with the valid particle types
        '''

        print('\ntest_correct_particle_types: Check that a variety of particles are set up with the valid particle types')

        valid_particle_types = [swiftest.constants.PL_TINY_TYPE_NAME, swiftest.constants.PL_TYPE_NAME, swiftest.constants.TP_TYPE_NAME, swiftest.constants.CB_TYPE_NAME]
        sim = swiftest.Simulation(simdir=self.simdir)

        sim.add_solar_system_body(['Sun'])

        # semi-interacting body
        sim.add_body(name=['pl_massive', 'pl_semi_int'], a=[1.5, 2.0], e=[0.1, 0.1], mass=[1e-11, 1e-6], radius=[1e-5, 1e-6])

        # test particle
        sim.add_body(name='tp', a=2.0, e=0.01)


        for particle_type in sim.init_cond.particle_type.values:
            self.assertTrue(particle_type in valid_particle_types, msg = f'Invalid particle type set up for {particle_type}')

        return 
         
if __name__ == '__main__':
    os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
    unittest.main()