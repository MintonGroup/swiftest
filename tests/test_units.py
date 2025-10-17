"""
Copyright 2025 - David Minton.

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

import astropy.constants as const
import numpy as np

import swiftest


class TestSwiftestUnits(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir = tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name

    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup()

    def test_unit_conversions(self):
        """
        Test setting standard and custom units and checks that the conversions are correct
        """
        sim = swiftest.Simulation(simdir=self.simdir)
        MSun = const.M_sun.value
        sim.set_parameter(MU="Msun")
        self.assertAlmostEqual(sim.MU_name, "MSun")
        self.assertAlmostEqual(sim.MU2KG / MSun, np.longdouble(1.0))
        self.assertAlmostEqual(sim.KG2MU * MSun, np.longdouble(1.0))
        self.assertAlmostEqual(sim.MU2KG, sim.param["MU2KG"])

        sim.set_parameter(MU2KG=1e-3)
        self.assertAlmostEqual(sim.MU_name, "MU")
        self.assertAlmostEqual(sim.MU2KG, np.longdouble(1e-3))
        self.assertAlmostEqual(sim.KG2MU, np.longdouble(1000.0))
        self.assertAlmostEqual(sim.MU2KG, sim.param["MU2KG"])

        sim.set_parameter(DU="cm")
        self.assertAlmostEqual(sim.DU_name, "cm")
        self.assertAlmostEqual(sim.DU2M, np.longdouble(1e-2))
        self.assertAlmostEqual(sim.M2DU, np.longdouble(100.0))

        sim.set_parameter(DU="km")
        self.assertAlmostEqual(sim.DU_name, "km")
        self.assertAlmostEqual(sim.DU2M, np.longdouble(1e3))


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    unittest.main()
