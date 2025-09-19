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


class TestCollisions(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir = tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name

    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup()

    def test_solar_impact(self):
        """
        Tests that impacts into the central body work correctly.
        """
        sim = swiftest.Simulation(simdir=self.simdir, compute_conservation_values=True, rotation=True, collision_model="merge")

        # Add the modern planets and the Sun using the JPL Horizons Database.
        sim.add_solar_system_body(["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"])
        runargs = {"tstart": 0.0, "tstop": 5e-2, "dt": 0.0001, "istep_out": 1, "dump_cadence": 0}

        density = 3000.0 * sim.KG2MU / sim.M2DU**3

        # Make a massive body with a periapsis inside the Sun's radius
        q = 0.01 * swiftest.RSun * sim.M2DU
        a = 0.1
        e = 1.0 - q / a
        M = 0.1 * sim.init_cond.sel(name="Earth")["mass"].values[0]
        R = (3 * M / (4 * np.pi * density)) ** (1.0 / 3.0)
        rot = 4 * sim.init_cond.sel(name="Earth")["rot"]
        sim.add_body(
            name="Sundiver", a=a, e=e, inc=0.0, capom=0.0, omega=0.0, capm=180.0, mass=M, radius=R, Ip=[0.4, 0.4, 0.4], rot=rot
        )

        for mtiny in [2 * M, 0.5 * M]:
            sim.clean()
            sim.run(**runargs, mtiny=mtiny, integrator="symba")

            # Check that the collision actually happened
            self.assertEqual(sim.collisions.collision_id.size, 1)

            # Check that angular momentum is conserved
            ds = sim.collisions.sel(collision_id=1)
            ds["Ltot"] = ds.L_orbit + ds.L_rot
            ds["Ltot_mag"] = ds.Ltot.magnitude()
            dLtot = ds.Ltot_mag.diff("stage").values[0]
            self.assertLess(dLtot, 1e-6, msg=f"Angular momentum not conserved: {dLtot}")

            # Check that energy was lost
            dEtot = ds.TE.diff("stage").values[0]
            self.assertLess(dEtot, 0, msg=f"Energy not lost: {dEtot}")

        # Test that massive bodies can be discarded in RMVS
        sim.run(**runargs, integrator="rmvs")
        # Check that the collision actually happened
        self.assertEqual(sim.collisions.collision_id.size, 1, msg="Collision not detected in RMVS")

        # Test that test particles are discarded in RMVS and SyMBA
        sim = swiftest.Simulation(simdir=self.simdir)
        sim.add_solar_system_body(["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"])
        sim.add_body(name="Sundiver", a=a, e=e, inc=0.0, capom=0.0, omega=0.0, capm=180.0)
        for integrator in ["rmvs", "symba", "whm", "helio"]:
            sim.clean()
            sim.run(**runargs, integrator=integrator)
            self.assertEqual(sim.collisions.collision_id.size, 1, msg=f"Collision not detected in {integrator}")
            self.assertEqual(
                sim.collisions.sel(collision_id=1).regime.values,
                "Central Body Impact",
                msg=f"{integrator}: Wrong regime: {sim.collisions.sel(collision_id=1).regime.values}",
            )

        return

    def test_escape(self):
        """
        Tests that escaping bodies are handled correctly.
        """
        sim = swiftest.Simulation(simdir=self.simdir, compute_conservation_values=True, rotation=True)

        # Add the modern planets and the Sun using the JPL Horizons Database.
        sim.add_solar_system_body(["Sun", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"])

        density = 3000.0 * sim.KG2MU / sim.M2DU**3

        # Make a hyperbolic body
        a = 0.1
        e = 1.1
        M = 1e-4 * swiftest.MEarth * sim.KG2MU
        R = (3 * M / (4 * np.pi * density)) ** (1.0 / 3.0)
        rot = 2 * sim.init_cond.sel(name="Saturn")["rot"]
        sim.add_body(
            name="Escapee", a=a, e=e, inc=0.0, capom=0.0, omega=45.0, capm=180.0, mass=M, radius=R, Ip=[0.4, 0.4, 0.4], rot=rot
        )

        runargs = {"tstart": 0.0, "tstop": 1000.0, "dt": 0.05, "istep_out": 200, "dump_cadence": 0}

        for mtiny in [2 * M, 0.5 * M]:
            # Try with semi-interacting body
            sim.clean()
            sim.run(**runargs, mtiny=mtiny, integrator="symba")

            # Check that the escape event was recorded
            self.assertEqual(sim.collisions.collision_id.size, 1)
            self.assertEqual(
                sim.collisions.sel(collision_id=1).regime.values,
                "Ejected",
                msg=f"mtiny/M: {mtiny / M}: Wrong regime: {sim.collisions.sel(collision_id=1).regime.values}",
            )

            # Check that angular momentum is conserved
            ds = sim.collisions.sel(collision_id=1)
            ds["Ltot"] = ds.L_orbit + ds.L_rot
            ds["Ltot_mag"] = ds.Ltot.magnitude()
            dLtot = ds.Ltot_mag.diff("stage").values[0]
            self.assertLess(dLtot, 1e-6, msg=f"Mtiny/M: {mtiny / M} Angular momentum not conserved: {dLtot}")

            # Check that energy was lost
            dEtot = ds.TE.diff("stage").values[0]
            self.assertLessEqual(dEtot, 0, msg=f"Energy not lost: {dEtot}")

        # Test that massive bodies can be discarded in all integrators
        for integrator in ["rmvs", "symba", "whm", "helio"]:
            sim.clean()
            sim.run(**runargs, integrator=integrator)
            # Check that the collision actually happened
            self.assertEqual(sim.collisions.collision_id.size, 1, msg=f"Collision not detected in {integrator}")

        # Test that test particles are discarded in all integrators
        sim = swiftest.Simulation(simdir=self.simdir)
        sim.add_solar_system_body(["Sun", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"])
        sim.add_body(name="Escapee", a=a, e=e, inc=0.0, capom=0.0, omega=45.0, capm=180.0)
        # Test that massive bodies can be discarded in all integrators
        for integrator in ["rmvs", "symba", "whm", "helio"]:
            sim.clean()
            sim.run(**runargs, integrator=integrator)
            # Check that the collision actually happened
            self.assertEqual(sim.collisions.collision_id.size, 1, msg=f"Collision not detected in {integrator}")
            self.assertEqual(
                sim.collisions.sel(collision_id=1).regime.values,
                "Ejected",
                msg=f"{integrator}: Wrong regime: {sim.collisions.sel(collision_id=1).regime.values}",
            )

        return

    def test_merge(self):
        """
        Tests that merging bodies are handled correctly.
        """
        sim = swiftest.Simulation(simdir=self.simdir, compute_conservation_values=True, rotation=True, collision_model="merge")
        sim.add_solar_system_body("Sun")

        name = ["Target", "Projectile"]
        rh = [np.array([1.0, -5.0e-05, 0.0]), np.array([1.0, 5.0e-05, 0.0])]
        vh = [np.array([0.00, 6.280005, 0.0]), np.array([0.00, 3.90, 0.0])]
        rot = [np.array([0.0, 0.0, 1.0e5]), np.array([0.0, 0.0, -5e5])]
        Gmass = [1e-7, 1e-9]
        density = 3000.0 * sim.KG2MU / sim.M2DU**3
        radius = [((GM / sim.GU) / (4.0 / 3.0 * np.pi * density)) ** (1.0 / 3.0) for GM in Gmass]
        runargs = {"tstart": 0.0, "tstop": 2e-3, "dt": 5e-4, "istep_out": 1, "dump_cadence": 0}

        # Test that massive body merge works in SyMBA
        sim.add_body(name=name, rh=rh, vh=vh, rot=rot, Gmass=Gmass, radius=radius)

        for gmtiny in [2 * Gmass[1], 0.5 * Gmass[1]]:
            sim.clean()
            sim.run(**runargs, gmtiny=gmtiny, integrator="symba")

            # Check that the collision actually happened
            self.assertEqual(sim.collisions.collision_id.size, 1)

            # Check that the collision was a merge
            self.assertEqual(
                sim.collisions.sel(collision_id=1).regime.values,
                "Merge",
                msg=f"Wrong regime: {sim.collisions.sel(collision_id=1).regime.values}",
            )

            # Check that angular momentum is conserved
            ds = sim.collisions.sel(collision_id=1)
            ds["Ltot"] = ds.L_orbit + ds.L_rot
            ds["Ltot_mag"] = ds.Ltot.magnitude()
            dLtot = ds.Ltot_mag.diff("stage").values[0]
            self.assertAlmostEqual(dLtot, 0, places=8, msg=f"Angular momentum not conserved: {dLtot}")

            # Check that energy was lost
            dEtot = ds.TE.diff("stage").values[0]
            self.assertLess(dEtot, 0, msg=f"Energy not lost: {dEtot}")

    def test_multi_collision(self):
        """
        Tests that multiple collisions are handled correctly.
        """
        sim = swiftest.Simulation(simdir=self.simdir, collision_model="merge")

        sim.add_solar_system_body(["Sun", "Earth"])

        # add bodies that both impact the CB at the same time
        sim.add_body(name=["close"], rh=[0, 0, sim.param["CHK_RMIN"] * 1.1], vh=[0, 0, 0.0001])
        sim.add_body(name=["close2"], rh=[0, 0, sim.param["CHK_RMIN"] * 1.11], vh=[0, 0, 0.0001])

        try:
            sim.run(tstop=10, dt=0.01, tstep_out=1)
        except Exception as e:
            self.fail(f"Multiple collisions not handled correctly: {e}")

        self.assertEqual(sim.collisions.collision_id.size, 2)

    def test_multi_pltp_collision(self):
        """
        Tests that multiple pl-tp collisions are handled correctly.
        """
        sim = swiftest.Simulation(simdir=self.simdir)  #

        sim.add_solar_system_body(["Sun", "Earth"])
        ncoll = 10
        nnotcoll = 10

        # Create bodies that will collide
        rh_E = sim.init_cond.sel(time=0, name="Earth")["rh"].values
        vh_E = sim.init_cond.sel(time=0, name="Earth")["vh"].values
        rad_E = sim.init_cond.sel(time=0, name="Earth")["radius"].values.item()

        rh_tp = []
        vh_tp = []
        for i in range(ncoll):
            rh_tp.append(rh_E + [0, -20 * (i + 1) * rad_E, 0])
            vh_tp.append(vh_E + [0.0, 1e-2, 0])

        sim.add_body(rh=rh_tp, vh=vh_tp)

        # Create bodies that will not collide
        a_tp = []
        a_tp = [10.0 + 0.1 * i for i in range(nnotcoll)]
        sim.add_body(a=a_tp)

        try:
            sim.run(tstop=100, dt=0.001, tstep_out=10, integrator="rmvs")
        except Exception as e:
            self.fail(f"Multiple pl-tp collisions not handled correctly: {e}")

        self.assertEqual(sim.collisions.collision_id.size, ncoll)


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    unittest.main()
