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

import numpy as np

import swiftest


class TestXV2EL2XV(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir = tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name

    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup()

    def test_xv2el2xv_dims(self):
        """
        Check that the xv2el and el2xv converts correctly with the right dimensions
        """
        print("\ntest_xv2el2xv_dims: Check that the xv2el and el2xv converts with the right dimensions")

        sim = swiftest.Simulation(simdir=self.simdir, tstop=0.5, dt=0.1)

        sim.add_solar_system_body(["Sun", "Mars"])
        xv_dims = sim.data.rh.dims
        el_dims = sim.data.a.dims

        # make dummy changes to the bodies
        sim.modify_body(name="Sun", c_lm=np.ones([2, 3, 3]))
        sim.modify_body(name="Mars", e=0.1)
        sim.data.el2xv()
        sim.data.xv2el()

        # check that rh and vh have the same dimensions
        self.assertTrue(sim.data.rh.dims == xv_dims)
        self.assertTrue(sim.data.vh.dims == xv_dims)

        # check that all orbital elements have the same dimensions
        self.assertTrue(sim.data.a.dims == el_dims)
        self.assertTrue(sim.data.e.dims == el_dims)
        self.assertTrue(sim.data.inc.dims == el_dims)
        self.assertTrue(sim.data.omega.dims == el_dims)
        self.assertTrue(sim.data.capom.dims == el_dims)
        self.assertTrue(sim.data.capm.dims == el_dims)
        self.assertTrue(sim.data.capf.dims == el_dims)
        self.assertTrue(sim.data.cape.dims == el_dims)
        self.assertTrue(sim.data.varpi.dims == el_dims)
        self.assertTrue(sim.data.lam.dims == el_dims)

        return

    def test_type_mismatch(self):
        """
        Checks that dtypes are converted properly before passing arguments to the el2xv or xv2el methods
        """
        sim = swiftest.Simulation(simdir=self.simdir)
        sim.add_solar_system_body(["Sun", "Mercury"])
        GMcb = swiftest.GMSun * sim.M2DU**3 / sim.S2TU**2
        sim.data.xv2el(GMcb)
        sim.data.el2xv(GMcb)


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    unittest.main()
