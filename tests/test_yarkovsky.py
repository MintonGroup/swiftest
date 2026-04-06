"""
Copyright 2026 - David Minton.

This file is part of Swiftest.
Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Swiftest.
If not, see: https://www.gnu.org/licenses.
"""

import os
import tempfile
import unittest

import numpy as np

import swiftest


class TestCollisions(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir = tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name

    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup()

    def test_yarkovsky(self):
        """
        Tests that Yarkovsky forces are handles correctly.
        """

        # Asteroid 7231 Veritas clone with values from Ferich, et al (2022) and Carruba, et al (2017)
        a = 3.165802 # AU
        e = 0
        inc = 0
        density = 1300 # kg/m^3
        radius = 100.0 # m
        mass = 4.0/3.0 * np.pi * density * radius**3 # kg
        rot_period = 6.0 # h
        rot_period *= 60.0 * 60.0 # h to s
        rot_mag = 360.0 / rot_period # deg/s
        # obliquities = [0.0, 60.0, 90.0, 180.0] # degrees
        albedo = 0.07
        emissivity = 0.9
        gamma = 100.0 # SI units
        rot_k = 0.25

        # set up the simulation object
        sim = swiftest.Simulation(simdir = self.simdir,
                                  integrator = 'symba',
                                  DU = 'AU', TU = 'y', MU = 'Msun',
                                  tstop = 1e3,
                                  dt = 0.05,
                                  yarkovsky = True,
                                  general_relativity = False,
                                  verbose = False,
                                  dump_cadence = 0)
        sim.add_solar_system_body(name = ['Sun'], align_to_central_body_rotation=True)

        # check that data inputs are handled correctly 
        # if one value is missing, an error should be thrown. If multiple values are missing, this test should extend to cover that.
        
        # for this case, we take an obliquity = 0 deg
        obliquity = 0.0
        rot = np.array([np.sin(np.deg2rad(obliquity)), 0, np.cos(np.deg2rad(obliquity))])

        # missing albedo

        with self.assertRaises(Exception):
            sim.add_body(name = 'Veritas', radius = radius / sim.DU2M, mass = mass / sim.MU2KG,
                            rot = rot * rot_mag * sim.TU2S,
                            a = a, # already in AU
                            e = 0, inc = 0, capom = 0.0, omega = 0.0, capm = 0.0, 
                            emissivity = emissivity, 
                            gamma = gamma * (sim.TU2S**(5.0/2)) / sim.MU2KG, 
                            rot_k = rot_k)

        # missing emissivity
        with self.assertRaises(Exception):
            sim.add_body(name = 'Veritas', radius = radius / sim.DU2M, mass = mass / sim.MU2KG,
                            rot = rot * rot_mag * sim.TU2S,
                            a = a, # already in AU
                            e = 0, inc = 0, capom = 0.0, omega = 0.0, capm = 0.0,
                            albedo = albedo, 
                            gamma = gamma * (sim.TU2S**(5.0/2)) / sim.MU2KG, 
                            rot_k = rot_k)

        # missing rotational k constant
        with self.assertRaises(Exception):
            sim.add_body(name = 'Veritas', radius = radius / sim.DU2M, mass = mass / sim.MU2KG,
                            rot = rot * rot_mag * sim.TU2S,
                            a = a, # already in AU
                            e = 0, inc = 0, capom = 0.0, omega = 0.0, capm = 0.0,
                            albedo = albedo, 
                            emissivity = emissivity, 
                            gamma = gamma * (sim.TU2S**(5.0/2)) / sim.MU2KG) 

        # missing thermal inertia (gamma)
        with self.assertRaises(Exception):
            sim.add_body(name = 'Veritas', radius = radius / sim.DU2M, mass = mass / sim.MU2KG,
                            rot = rot * rot_mag * sim.TU2S,
                            a = a, # already in AU
                            e = 0, inc = 0, capom = 0.0, omega = 0.0, capm = 0.0,
                            albedo = albedo, 
                            emissivity = emissivity, 
                            rot_k = rot_k)

        # correct yarkovsky inputs, but missing rotation/obliquity
        obliquities = [0.0, 60.0, 90.0, 180.0] # degrees

        # check da/dt has the correct direction and value for various obliquities and sizes
    
        return


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    unittest.main()