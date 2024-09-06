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
        
    def test_rotation_collision(self):
        '''
        Check that the rotation rate is conserved correctly when a disruptioncollision occurs

        TEMPORARY TEST: to be incorporated in test_collision.py or test_fraggle.py
        '''

        sim = swiftest.Simulation(simdir=self.simdir, rotation=True, compute_conservation_values=True,
                                  MU2KG = 1e15, DU2M = 1e3, TU2S = 24 * 3600.0,
                                  tstop = 0.03,
                                  dt = 0.003)
        sim.add_solar_system_body(['Mars', 'Phobos'])

        # add Deimos and impactor
        # parameters taken from 10.0inc_5000y_semi_int_p simulation
        sim.add_body(name = ['Deimos', 'impactor'], 
                     Gmass = [8.96820886e+05, 4.78022953e+00], 
                     radius = [6.20305087, 0.11515388],
                     rh = [[-21619.36333433,  -8245.07366662,  -3855.71806009], [-21618.29800762,  -8240.27915225,  -3855.79044167]], 
                     vh = [[  40486.26555708, -109312.40507476,    6532.01416311], [  40418.61067766, -108069.91067059,    6485.34056135]], 
                     rot = [[-13.38759186,  -0.42146959, 284.72074881], [  0.        ,   0.        ,   0.        ]])
        
        sim.run()

        # check that the rotation rate is conserved
        rot_before = sim.collisions.sel(collision_id = 1, collision_body = 1, stage = 'before').rot.magnitude().values # rot_before = 285.03562945368185
        rot_after = sim.collisions.sel(collision_id = 1, collision_body = 1, stage = 'after').rot.magnitude().values # 285.03804937240477

        print(rot_before) 
        print(rot_after)

        self.assertAlmostEqual(rot_before, rot_after, places=5)

        return


if __name__ == '__main__':
    os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
    unittest.main()