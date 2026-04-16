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


class TestRadiation(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir = tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name

    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup()

    def test_yarkovsky_forces(self):
        """
        Tests that Yarkovsky forces are handles correctly.
        """
        # Asteroid 7231 Veritas clone with values from Ferich, et al (2022) and Carruba, et al (2017)
        a = 3.165802  # AU
        e = 0.0
        inc = 0.0
        capom = (0.0,)
        omega = (0.0,)
        capm = (0.0,)
        density = 1300  # kg/m^3
        radius = 100.0  # m
        mass = 4.0 / 3.0 * np.pi * density * radius**3  # kg
        rot_period = 6.0  # h
        rot_period *= 60.0 * 60.0  # h to s
        rot_mag = 360.0 / rot_period  # deg/s
        # obliquities = [0.0, 60.0, 90.0, 180.0] # degrees
        albedo = 0.07
        emissivity = 0.9
        gamma = 4000.0  # SI units
        rot_k = 0.25

        # set up the simulation object
        sim = swiftest.Simulation(
            simdir=self.simdir,
            integrator="symba",
            DU="AU",
            TU="y",
            MU="Msun",
            tstop=1e3,
            dt=0.05,
            yarkovsky=True,
            general_relativity=False,
            verbose=False,
            dump_cadence=0,
        )
        sim.add_solar_system_body(name=["Sun"], align_to_central_body_rotation=True)

        # check that data inputs are handled correctly
        # if one value is missing, an error should be thrown. If multiple values are missing, this test should extend to cover that.

        # for this case, we take an obliquity = 0 deg
        obliquity = 0.0
        rot_unit_vector = np.array([np.sin(np.deg2rad(obliquity)), 0, np.cos(np.deg2rad(obliquity))])

        # missing albedo

        with self.assertRaises(ValueError):
            sim.add_body(
                name="Veritas",
                radius=radius / sim.DU2M,
                mass=mass / sim.MU2KG,
                rot=rot_unit_vector * rot_mag * sim.TU2S,
                a=a,  # already in AU
                e=e,
                inc=inc,
                capom=capom,
                omega=omega,
                capm=capm,
                emissivity=emissivity,
                gamma=gamma * (sim.TU2S ** (5.0 / 2)) / sim.MU2KG,
                rot_k=rot_k,
            )

        # missing emissivity
        with self.assertRaises(ValueError):
            sim.add_body(
                name="Veritas",
                radius=radius / sim.DU2M,
                mass=mass / sim.MU2KG,
                rot=rot_unit_vector * rot_mag * sim.TU2S,
                a=a,  # already in AU
                e=e,
                inc=inc,
                capom=capom,
                omega=omega,
                capm=capm,
                albedo=albedo,
                gamma=gamma * (sim.TU2S ** (5.0 / 2)) / sim.MU2KG,
                rot_k=rot_k,
            )

        # missing rotational k constant
        with self.assertRaises(ValueError):
            sim.add_body(
                name="Veritas",
                radius=radius / sim.DU2M,
                mass=mass / sim.MU2KG,
                rot=rot_unit_vector * rot_mag * sim.TU2S,
                a=a,  # already in AU
                e=0,
                inc=0,
                capom=0.0,
                omega=0.0,
                capm=0.0,
                albedo=albedo,
                emissivity=emissivity,
                gamma=gamma * (sim.TU2S ** (5.0 / 2)) / sim.MU2KG,
            )

        # missing thermal inertia (gamma)
        with self.assertRaises(ValueError):
            sim.add_body(
                name="Veritas",
                radius=radius / sim.DU2M,
                mass=mass / sim.MU2KG,
                rot=rot_unit_vector * rot_mag * sim.TU2S,
                a=a,  # already in AU
                e=e,
                inc=inc,
                capom=capom,
                omega=omega,
                capm=capm,
                albedo=albedo,
                emissivity=emissivity,
                rot_k=rot_k,
            )

        # check da/dt has the correct direction and value for various obliquities and sizes

        obliquities = [60.0, 90.0]  # degrees
        radii = [10.0, 100.0]  # m

        # pre-tested results for da at t = 1,000 y in AU for the given obliquities and radii

        da_pre_tested = {
            "10.0": {
                "0.0": 2.807739638921447e-05,
                "60.0": 1.562845768354748e-06,
                "90.0": -2.302262337128269e-05,
                "180.0": -6.640618802133957e-05,
            },
            "100.0": {
                "0.0": 2.8077498379630583e-06,
                "60.0": 1.5627818239494218e-07,
                "90.0": -2.302270235343684e-06,
                "180.0": -6.640592883133678e-06,
            },
        }

        delta = {"10.0": 1e-6, "100.0": 1e-7}

        sim.clean()

        for radius in radii:
            for obliquity in obliquities:
                mass = 4.0 / 3.0 * np.pi * density * radius**3  # kg
                rot_unit_vector = np.array([np.sin(np.deg2rad(obliquity)), 0, np.cos(np.deg2rad(obliquity))])

                sim = swiftest.Simulation(
                    simdir=self.simdir,
                    integrator="symba",
                    DU="AU",
                    TU="y",
                    MU="Msun",
                    tstop=1e3,
                    dt=0.05,
                    yarkovsky=True,
                    general_relativity=False,
                    verbose=False,
                    dump_cadence=0,
                )

                sim.add_solar_system_body(name=["Sun"], align_to_central_body_rotation=True)

                sim.add_body(
                    name="Veritas",
                    radius=radius / sim.DU2M,
                    mass=mass / sim.MU2KG,
                    rot=rot_unit_vector * rot_mag * sim.TU2S,
                    a=a,  # already in AU
                    e=e,
                    inc=inc,
                    capom=capom,
                    omega=omega,
                    capm=capm,
                    albedo=albedo,
                    emissivity=emissivity,
                    gamma=gamma * (sim.TU2S ** (5.0 / 2)) / sim.MU2KG,
                    rot_k=rot_k,
                )

                sim.run()

                da = sim.data.sel(name="Veritas", time=sim.data.time.values[-1]).a - sim.data.sel(name="Veritas", time=0).a
                # compare with pre-tested results
                self.assertAlmostEqual(da, da_pre_tested[f"{radius}"][f"{obliquity}"], delta=delta[f"{radius}"])
                sim.clean()
        return


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    unittest.main()
