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

import datetime
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from astroquery.jplhorizons import Horizons
from numpy.random import default_rng

import swiftest

# if python version is >=3.11 use contextlib.chdir, otherwise we need to define it manually
if sys.version_info >= (3, 11):
    from contextlib import chdir
else:
    from contextlib import contextmanager

    @contextmanager
    def chdir(path):
        oldpwd = Path.cwd()
        os.chdir(str(path))
        try:
            yield
        finally:
            os.chdir(str(oldpwd))


rng = default_rng(seed=123)

major_bodies = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
param = {}


class TestSwiftestIO(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir = tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name

    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup()

    def test_jpl_parser(self):
        """
        Tests that the JPL Horizons parser is able to read data and process it correctly for a variety of different cases.
        """
        print("\ntest_jpl_parser")
        names = major_bodies + [
            "Kleopatra",
            "Ceres",
            "Vesta",
            "Pallas",
            "Hygiea",
            "Eris",
            "Dysnomia",
            "Pluto",
            "Charon",
            "Haumea",
            "Quaoar",
            "Orcus",
            "Io",
            "Europa",
            "Titan",
            "Triton",
            "Umbriel",
            "Moon",
            "Phobos",
            "Deimos",
        ]
        test_keys = ["mass", "radius", "rot"]

        body_props = {
            "Sun": {
                "mass": 1.988409870698051e30,
                "radius": 695700000.0,
                "rot": [2.035154582041381e-05, -7.037127139847693e-05, 0.000149334181813172],
            },
            "Mercury": {
                "mass": 3.3010006367708975e23,
                "radius": 2439400.0,
                "rot": [6.4872186125818495e-06, -5.800032565782477e-06, 7.051241078428856e-05],
            },
            "Venus": {
                "mass": 4.867305814842006e24,
                "radius": 6051840.0,
                "rot": [-3.2045928611476063e-07, -1.8640854693792636e-07, -1.7141180404335365e-05],
            },
            "Earth": {
                "mass": 6.045626290427461e24,
                "radius": 6371010.0,
                "rot": [1.1179382676404839e-05, 0.001662031095337799, 0.0038332533328792147],
            },
            "Mars": {
                "mass": 6.416908921385015e23,
                "radius": 3389920.0,
                "rot": [0.0018117998589411537, -0.00022719827850226695, 0.003627605423316273],
            },
            "Jupiter": {
                "mass": 1.8981246258034552e27,
                "radius": 69911000.0,
                "rot": [-0.0001470359129430991, -0.00036060109828205247, 0.010067934171122136],
            },
            "Saturn": {
                "mass": 5.683173701212113e26,
                "radius": 58232000.0,
                "rot": [0.0008024029853603659, 0.004339610305746439, 0.008281723231138592],
            },
            "Uranus": {
                "mass": 8.680987153709004e25,
                "radius": 25362000.0,
                "rot": [0.0012296933085501616, 0.005614774719295097, -0.0007793665749130923],
            },
            "Neptune": {
                "mass": 1.024092409690904e26,
                "radius": 24624000.0,
                "rot": [0.0022373803628068065, -0.0019444371136576489, 0.005453805355586545],
            },
            "Kleopatra": {"mass": None, "radius": 61000.0, "rot": [0.0, 0.0, 0.018570102135561744]},
            "Ceres": {
                "mass": 9.383515874323901e20,
                "radius": 469700.0,
                "rot": [0.0015876541706678936, 0.00031457177647549044, 0.010900790208524158],
            },
            "Vesta": {
                "mass": 2.5902701406889116e20,
                "radius": 261385.0,
                "rot": [0.008727916415949765, -0.004872733218667549, 0.01582674491187599],
            },
            "Pallas": {
                "mass": 2.0421617248250754e20,
                "radius": 256500.0,
                "rot": [0.010719281540412906, 0.006120315206702517, -0.0033835644344280815],
            },
            "Hygiea": {"mass": 1.048799125001873e20, "radius": 203560.0, "rot": [0.0, 0.0, 0.007231703789412787]},
            "Eris": {"mass": 1.6018159207707176e22, "radius": 1163000.0, "rot": [0.0, 0.0, 0.0038610038610038607]},
            "Dysnomia": {"mass": 1.4323599478597008e20, "radius": 350000.0, "rot": [0.0, 0.0, 0.0]},
            "Pluto": {
                "mass": 1.461465621862967e22,
                "radius": 1188300.0,
                "rot": [-0.0004422702542154748, 0.00040738901587759955, -0.00025295769402782636],
            },
            "Charon": {
                "mass": 1.5896798166099818e21,
                "radius": 606000.0,
                "rot": [-0.00044226841369145176, 0.0004073873205131401, -0.0002529566413350254],
            },
            "Haumea": {"mass": 3.96165890055886e21, "radius": 1050000.0, "rot": [0.0, 0.0, 0.02554055955790313]},
            "Quaoar": {"mass": 1.4323599478597007e21, "radius": 543500.0, "rot": [0.0, 0.0, 0.011312985044233772]},
            "Orcus": {"mass": 5.777329911751045e20, "radius": 455000.0, "rot": [0.0, 0.0, 0.0075826508947528055]},
            "Io": {
                "mass": 8.929648802121572e22,
                "radius": 1821490.0,
                "rot": [-3.4523895639215165e-05, -8.533680868081025e-05, 0.0023522483405567646],
            },
            "Europa": {
                "mass": 4.7985737830184444e22,
                "radius": 1560800.0,
                "rot": [-2.2939528727418718e-05, -3.458615882284886e-05, 0.0011729749379113096],
            },
            "Titan": {
                "mass": 1.3451807680206165e23,
                "radius": 2575500.0,
                "rot": [2.308320022417276e-05, 0.00012070667598529457, 0.00023060562577272475],
            },
            "Triton": {
                "mass": 2.1402918658136437e22,
                "radius": 1352600.0,
                "rot": [0.00037532312265244825, -0.0004021490666056507, 0.0004472809262462142],
            },
            "Umbriel": {
                "mass": 1.2795349325022852e21,
                "radius": 584700.0,
                "rot": [-0.000214692403316975, -0.0009726715883856238, 0.00013706415384694278],
            },
            "Moon": {
                "mass": 7.345789170399893e22,
                "radius": 1737530.0,
                "rot": [2.842122641579648e-06, 3.0323633838439856e-06, 0.00015244753493999662],
            },
            "Phobos": {
                "mass": 1.08e16,
                "radius": 11166.666666666666,
                "rot": [0.005818491849485793, -0.0004851980251664484, 0.011684025938372589],
            },
            "Deimos": {
                "mass": 1800000000000000.2,
                "radius": 6300.0,
                "rot": [0.001331405578408647, -0.0001912658137505407, 0.0030123632841272847],
            },
        }

        for name in names:
            body = swiftest.init_cond.get_solar_system_body(name=name)
            test_prop = {k: body[k] for k in test_keys if k in body}
            for k, v in body_props[name].items():
                if k in test_keys:
                    if v is None:
                        self.assertIsNone(test_prop[k], msg=f"Error in {name} {k}: {test_prop[k]} != {v}")
                    else:
                        self.assertTrue(np.allclose(test_prop[k], v), msg=f"Error in {name} {k}: {test_prop[k]} != {v}")
        return

    def test_xv2el2xv(self):
        """
        Tests that the functions xv2el and el2xv are able to convert between position-velocity and orbital elements without any exceptions being raised.
        """
        print("\ntest_xv2el2xv")
        # Generate a set of random position-velocity vectors
        from swiftest.core import el2xv, xv2el

        # Test that we can reliably convert between orbital elements and state vectors
        n = 1000
        mu = np.ones(n)
        a = rng.uniform(0.1, 1.5, n)
        e = rng.uniform(0.0, 2.0, n)
        inc = rng.uniform(0.0, 180.0, n)
        capom = rng.uniform(0.0, 360.0, n)
        omega = rng.uniform(0.0, 360.0, n)
        capm = rng.uniform(0.0, 360.0, n)

        rh, vh = el2xv(mu, a, e, inc, capom, omega, capm)
        a2, e2, inc2, capom2, omega2, capm2, varpi, lam, f, cape, capf = xv2el(mu, rh, vh)

        # Check that the original and converted position-velocity vectors are the same
        self.assertTrue(np.allclose(a, a2), msg=f"Error converting a: {a}, {a2}")
        self.assertTrue(np.allclose(e, e2), msg=f"Error converting a: {e}, {e2}")
        self.assertTrue(np.allclose(inc, inc2), msg=f"Error converting a: {inc}, {inc2}")
        self.assertTrue(np.allclose(capom, capom2), msg=f"Error converting a: {capom}, {capom2}")
        self.assertTrue(np.allclose(omega, omega2), msg=f"Error converting a: {omega}, {omega2}")
        self.assertTrue(np.allclose(capm, capm2), msg=f"Error converting a: {capm}, {capm2}")
        return

    def test_gen_ic(self):
        """
        Tests that Swiftest is able to successfully generate a set of initial conditions in a file without any exceptions being raised.
        """
        print("\ntest_gen_ic")
        # Files that are expected to be generated:
        simdir_path = Path(self.simdir)
        file_list = [simdir_path, simdir_path / "param.in", simdir_path / "init_cond.nc"]

        sim = swiftest.Simulation(simdir=self.simdir)

        # Add the modern planets and the Sun using the JPL Horizons Database.
        sim.add_solar_system_body(major_bodies)
        sim.save()

        for f in file_list:
            self.assertTrue(Path(f).exists())

        print("\ntest_read_ic")
        sim2 = swiftest.Simulation(simdir=self.simdir, read_init_cond=True)
        # Add the modern planets and the Sun using the JPL Horizons Database.
        # Check if all names in Dataset read in from file match the expected list of names
        self.assertTrue((major_bodies == sim2.init_cond["name"]).all(), msg="Name mismatch in Dataset")

        # Check to see if all parameter values read in from file match the expected parameters saved when generating the file
        self.assertTrue(all(v == param[k] for k, v in sim2.param.items() if k in param))
        return

    def test_add_body_combinations(self):
        """
        Tests various combinations of arguments to the Simulation.add_body and Simulation.add_solar_system_body methods to ensure argument checking is working as intended.
        """
        print("\ntest_add_body_combinations")
        sim = swiftest.Simulation(simdir=self.simdir)

        # Test that we can can pass a single body:
        sim.add_solar_system_body(name="Sun")

        # Test that we can add multiple bodies:
        sim.add_solar_system_body(name=["Mercury", "Venus", "Earth"])

        # Test that we can pass by ephemeris_id alone
        sim.add_solar_system_body(ephemeris_id=["Mars", "500"])

        # Test that we can initialize a body with only the semimajor axis
        sim.add_body(a=1.0)
        self.assertEqual(sim.data.isel(name=-1).e.values[0], 0.0, msg="Failed to initialize body with only semimajor axis")
        self.assertEqual(sim.data.isel(name=-1).inc.values[0], 0.0, msg="Failed to initialize body with only semimajor axis")
        self.assertEqual(sim.data.isel(name=-1).capom.values[0], 0.0, msg="Failed to initialize body with only semimajor axis")
        self.assertEqual(sim.data.isel(name=-1).omega.values[0], 0.0, msg="Failed to initialize body with only semimajor axis")
        self.assertEqual(sim.data.isel(name=-1).capm.values[0], 0.0, msg="Failed to initialize body with only semimajor axis")

        # Test that we can input cartesian coordinates
        sim.add_body(mass=1.0, radius=1.0, rh=[1.0, 0.0, 0.0], vh=[0.0, 1.0, 0.0])

        # orbital elements without semimajor axis
        with self.assertRaises(ValueError):
            sim.add_body(e=0.0, inc=0.0, capom=0.0, omega=0.0, capm=0.0)

        # Mix of elements and cartesian inputs
        with self.assertRaises(ValueError):
            sim.add_body(a=1.0, e=0.0, rh=[1.0, 0.0, 0.0], vh=[0.0, 1.0, 0.0])

        # Position but not velocity
        with self.assertRaises(ValueError):
            sim.add_body(rh=[1.0, 0.0, 0.0])

        # Velocity but not position
        with self.assertRaises(ValueError):
            sim.add_body(vh=[0.0, 1.0, 0.0])

        # Add J2 and c_lm values
        with self.assertRaises(ValueError):
            sim.add_body(mass=1.0, radius=1.0, j2rp2=1.0e-6, c_lm=np.ones([2, 7]))

        # Wrong shape of c_lm
        with self.assertRaises(ValueError):
            sim.add_body(mass=1.0, radius=1.0, c_lm=[1.0, 0.0, 0.0])

        # Mismatched lengths of input arguments
        with self.assertRaises(ValueError):
            sim.add_solar_system_body(name=["Mercury", "Venus", "Earth"], ephemeris_id=["Mars", "500"])
        with self.assertRaises(ValueError):
            sim.add_body(a=[1.0, 2.0], e=0.0)
        with self.assertRaises(ValueError):
            sim.add_body(rh=[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], vh=[0.0, 1.0, 0.0])
        with self.assertRaises(ValueError):
            sim.add_body(rh=[1.0, 0.0], vh=[0.0, 1.0])

        # mass and Gmass at the same time
        with self.assertRaises(ValueError):
            sim.add_body(a=1.0, mass=1.0, radius=1.0, Gmass=4 * np.pi**2)

        # mass without radius
        with self.assertRaises(ValueError):
            sim.add_body(a=1.0, mass=1.0)

        return

    def test_mixed_element_cart_input(self):
        """
        Tests that we can mix orbital element and cartesian input.
        """
        print("\ntest_mixed_element_cart_input")
        sim = swiftest.Simulation(simdir=self.simdir)

        # Test that we can can pass a single body:
        sim.add_solar_system_body(name="Sun")

        sim.add_body(a=[1.0])
        sim.add_body(rh=[2.0, 0.0, 0.0], vh=[0.0, 1.0, 0.0])
        try:
            sim.run(tstop=0.02, dt=0.01)
        except Exception as e:
            self.fail(f"Failed to run simulation with mixed element and cartesian input: {e}")

        return

    def test_read_multi_dir(self):
        """
        Tests that Swiftest can generate a set of initial conditions, copy them into a new directory, and then run those from a different Simulation object (test inspired by Kaustub Anand's workflow).
        """
        print("\ntest_read_multi_dir")

        def copy_folder_contents(src, dst):
            """
            Copy the contents of the folder 'src' to the folder 'dst'.

            Parameters
            ----------
            src : str
                The path of the source directory.
            dst : str
                The path of the destination directory.

            Returns
            -------
            None
            """
            import shutil

            # Ensure the destination directory exists
            Path(dst).mkdir(parents=True)

            src_path = Path(src)
            dst_path = Path(dst)

            # Iterate over all items in the source directory
            for item in src_path.iterdir():
                target = dst_path / item.name
                if item.is_dir():
                    shutil.copytree(item, target, dirs_exist_ok=True)
                else:
                    shutil.copy(item, target)

            return

        base = Path(self.simdir)
        simdir1 = base / "sim1"
        simdir2 = base / "sim2"
        sim1 = swiftest.Simulation(simdir=str(simdir1), tstop=1.0, dt=0.01)
        sim1.add_solar_system_body(["Sun", "Mercury", "Venus", "Earth", "Mars"])
        sim1.save()
        copy_folder_contents(simdir1, simdir2)
        sim2 = swiftest.Simulation(read_init_cond=True, simdir=str(simdir2))
        try:
            sim2.run()
        except Exception as e:
            self.fail(f"Failed to run simulation from copied directory: {e}")

        return

    def test_planetocentric(self):
        """
        Tests that Swiftest is able to set up a simulation in a planetocentric frame and that the results are consistent with the expected values.
        """
        print("\ntest_planetocentric")

        cbname = "Mars"
        cbcenter = "@4"
        cb = swiftest.get_solar_system_body(cbname)
        sim = swiftest.Simulation(simdir=self.simdir, MU2KG=cb["mass"], DU2M=cb["radius"], TU="d")
        satname = ["Phobos", "Deimos"]
        sim.add_solar_system_body(cbname)
        sim.add_solar_system_body(satname)

        def get_jpl_data(id):
            tstart = datetime.date.fromisoformat(sim.ephemeris_date)
            tstep = datetime.timedelta(days=1)
            tend = tstart + tstep
            ephemerides_start_date = sim.ephemeris_date
            ephemerides_end_date = tend.isoformat()
            ephemerides_step = "1d"

            jpl = Horizons(
                id=id,
                location=cbcenter,
                epochs={"start": ephemerides_start_date, "stop": ephemerides_end_date, "step": ephemerides_step},
            )
            DCONV = swiftest.AU2M / sim.param["DU2M"]
            VCONV = (swiftest.AU2M / swiftest.JD2S) / (sim.param["DU2M"] / sim.param["TU2S"])

            vec_tab = jpl.vectors()
            rx = vec_tab["x"][0] * DCONV
            ry = vec_tab["y"][0] * DCONV
            rz = vec_tab["z"][0] * DCONV
            vx = vec_tab["vx"][0] * VCONV
            vy = vec_tab["vy"][0] * VCONV
            vz = vec_tab["vz"][0] * VCONV

            elem_tab = jpl.elements()

            a = elem_tab["a"][0] * DCONV
            e = elem_tab["e"][0]
            inc = elem_tab["incl"][0]
            capom = elem_tab["Omega"][0]
            omega = elem_tab["w"][0]
            capm = elem_tab["M"][0]

            rh = np.array([rx, ry, rz])
            vh = np.array([vx, vy, vz])
            elem = np.array([a, e, inc, capom, omega, capm])
            return rh, vh, elem

        for sat in satname:
            sim_rh = sim.data.sel(name=sat, time=0).rh.values
            sim_vh = sim.data.sel(name=sat, time=0).vh.values
            sim_elem = np.array(
                [sim.data[var].sel(name=sat, time=0).values.item() for var in ["a", "e", "inc", "capom", "omega", "capm"]]
            )

            jpl_rh, jpl_vh, jpl_elem = get_jpl_data(sat)
            self.assertTrue(np.allclose(sim_rh, jpl_rh), msg=f"Error in rh for {sat}: {sim_rh - jpl_rh}")
            self.assertTrue(np.allclose(sim_vh, jpl_vh), msg=f"Error in vh for {sat}: {sim_vh - jpl_vh}")
            self.assertTrue(np.allclose(sim_elem, jpl_elem, rtol=1e-3), msg=f"Error in elements for {sat}: {sim_elem - jpl_elem}")

        return

    def test_central_body_rotation(self):
        """
        Tests that Swiftest is able to rotate solar system bodies into the equatorial frame of a central body correctly.
        """
        print("\ntest_central_body_rotation")
        sim = swiftest.Simulation(simdir=self.simdir)
        # Precomputed values for [Mars, Phobos, and Deimos] in the ecliptic and Mars equatorial frames
        inc_vals = {
            "ecliptic": np.array([np.nan, 26.55406151, 24.06273241]),
            "mars equator": np.array([np.nan, 1.0786523, 2.69218629]),
        }
        rot_vals = {
            "ecliptic": np.array(
                [
                    [57176.05522852, -7169.83239366, 114478.52090685],
                    [183617.63838933, -15311.68519899, 368719.81695279],
                    [42015.96468119, -6035.89004401, 95062.95557518],
                ]
            ),
            "mars equator": np.array(
                [
                    [-3.08375547e-12, -9.09494702e-13, 1.28163331e05],
                    [-3.79488273e02, 7.76140989e03, 4.12121104e05],
                    [-4.88981793e03, -1.53941766e02, 1.03994254e05],
                ]
            ),
        }

        # No rotation should keep it in the ecliptic frame
        sim.add_solar_system_body(["Mars", "Phobos", "Deimos"])  # Default should be ecliptic
        sim_inc = sim.data.isel(time=0).inc.values
        sim_rot = sim.data.isel(time=0).rot.values

        inc_close = np.allclose(sim_inc, inc_vals["ecliptic"], rtol=1e-2, equal_nan=True)
        rot_close = np.allclose(sim_rot, rot_vals["ecliptic"], rtol=1e-2)

        self.assertTrue(inc_close, msg="Error in inclination 1")
        self.assertTrue(rot_close, msg="Error in rotation 1")

        # Rotating to match the Mars pole
        sim = swiftest.Simulation(simdir=self.simdir)

        sim.add_solar_system_body(["Mars", "Phobos", "Deimos"], align_to_central_body_rotation=True)
        sim_inc = sim.data.isel(time=0).inc.values
        sim_rot = sim.data.isel(time=0).rot.values

        inc_close = np.allclose(sim_inc, inc_vals["mars equator"], rtol=1e-2, equal_nan=True)
        rot_close = np.allclose(sim_rot, rot_vals["mars equator"], rtol=1e-2)

        self.assertTrue(inc_close, msg="Error in inclination 2")
        self.assertTrue(rot_close, msg="Error in rotation 2")

        sim = swiftest.Simulation(simdir=self.simdir)

        # Mix and match
        sim.add_solar_system_body("Mars")  # ecliptic
        sim.add_solar_system_body(["Phobos", "Deimos"], align_to_central_body_rotation=True)  # mars equator
        sim_inc = sim.data.isel(time=0).inc.values
        sim_rot = sim.data.isel(time=0).rot.values

        inc_close = np.allclose(sim_inc[1:], inc_vals["mars equator"][1:], rtol=1e-2)
        rot_close = np.allclose(sim_rot[0:1], rot_vals["ecliptic"][0:1], rtol=1e-2) and np.allclose(
            sim_rot[1:], rot_vals["mars equator"][1:], rtol=1e-2
        )

        self.assertTrue(inc_close, msg="Error in inclination 3")
        self.assertTrue(rot_close, msg="Error in rotation 3")

        sim = swiftest.Simulation(simdir=self.simdir)

        sim.add_solar_system_body("Mars", align_to_central_body_rotation=True)
        sim.add_solar_system_body(["Phobos", "Deimos"])  # ecliptic
        sim_inc = sim.data.isel(time=0).inc.values
        sim_rot = sim.data.isel(time=0).rot.values

        inc_close = np.allclose(sim_inc[1:], inc_vals["ecliptic"][1:], rtol=1e-2)
        rot_close = np.allclose(sim_rot[0:1], rot_vals["mars equator"][0:1], rtol=1e-2) and np.allclose(
            sim_rot[1:], rot_vals["ecliptic"][1:], rtol=1e-2
        )

        self.assertTrue(inc_close, msg="Error in inclination 4")
        self.assertTrue(rot_close, msg="Error in rotation 4")

        sim = swiftest.Simulation(simdir=self.simdir)

        sim.add_solar_system_body("Mars", align_to_central_body_rotation=True)
        sim.add_solar_system_body("Phobos")  # ecliptic
        sim.add_solar_system_body("Deimos", align_to_central_body_rotation=True)  # mars equator
        sim_inc = sim.data.isel(time=0).inc.values
        sim_rot = sim.data.isel(time=0).rot.values

        inc_close = np.allclose(sim_inc[1:2], inc_vals["ecliptic"][1:2], rtol=1e-2) and np.allclose(
            sim_inc[2:], inc_vals["mars equator"][2:], rtol=1e-2
        )
        rot_close = (
            np.allclose(sim_rot[0:1], rot_vals["mars equator"][0:1], rtol=1e-2)
            and np.allclose(sim_rot[1:2], rot_vals["ecliptic"][1:2], rtol=1e-2)
            and np.allclose(sim_rot[2:], rot_vals["mars equator"][2:], rtol=1e-2)
        )

        self.assertTrue(inc_close, msg="Error in inclination 5")
        self.assertTrue(rot_close, msg="Error in rotation 5")

        sim = swiftest.Simulation(simdir=self.simdir)

        sim.add_solar_system_body("Mars", align_to_central_body_rotation=True)
        sim.add_solar_system_body("Phobos", align_to_central_body_rotation=True)  # mars equator
        sim.add_solar_system_body("Deimos")  # ecliptic
        sim_inc = sim.data.isel(time=0).inc.values
        sim_rot = sim.data.isel(time=0).rot.values

        inc_close = np.allclose(sim_inc[1:2], inc_vals["mars equator"][1:2], rtol=1e-2) and np.allclose(
            sim_inc[2:], inc_vals["ecliptic"][2:], rtol=1e-2
        )
        rot_close = np.allclose(sim_rot[0:2], rot_vals["mars equator"][0:2], rtol=1e-2) and np.allclose(
            sim_rot[2:2], rot_vals["ecliptic"][2:], rtol=1e-2
        )

        self.assertTrue(inc_close, msg="Error in inclination 6")
        self.assertTrue(rot_close, msg="Error in rotation 6")

        return

    def test_remove_body(self):
        """
        Tests that Swiftest is able to remove a body from the simulation.
        """
        print("\ntest_remove_body")
        sim = swiftest.Simulation(simdir=self.simdir)
        sim.add_solar_system_body(["Sun", "Mercury", "Venus", "Earth", "Mars"])
        sim.remove_body("Mars")
        self.assertTrue("Mars" not in sim.data.name.values)

        with self.assertWarns(Warning) as cm:  # Using Warning as the base class, adjust if needed
            sim.remove_body("Arrakis")
        self.assertIn("Arrakis not found in the Dataset", str(cm.warning))

        with self.assertWarns(Warning) as cm:
            sim.remove_body(["Sun", "Mercury", "Venus", "Earth"])  # Remove all bodies
        self.assertIn("No bodies left in the Dataset", str(cm.warning))

        sim.add_solar_system_body(["Sun", "Mercury", "Venus", "Earth", "Mars"])
        sim.remove_body(id=4)
        self.assertTrue("Mars" not in sim.data.name.values)
        with self.assertWarns(Warning) as cm:  # Using Warning as the base class, adjust if needed
            sim.remove_body(id=10)
        self.assertIn("10 not found in the Dataset", str(cm.warning))
        return

    def test_modify_body(self):
        """
        Tests that Swiftest is able to modify the properties of a body in the simulation.
        """
        print("\ntest_modify_body")
        sim = swiftest.Simulation(simdir=self.simdir)
        sim.add_solar_system_body(["Sun", "Mercury"])
        sim.modify_body(name="Mercury", a=100.0)
        self.assertGreater(sim.data.sel(name="Mercury")["a"].values.item(), 99.0)

        self.assertTrue("j2rp2" in sim.data)

        c_lm = np.ones([2, 3, 3])
        sim.modify_body(name="Sun", c_lm=c_lm)

        is_close = np.allclose(sim.data.sel(name="Sun")["c_lm"].values, c_lm)
        self.assertTrue(is_close, msg=f"Error in c_lm: {sim.data.sel(name='Sun')['c_lm'].values - c_lm}")

        self.assertTrue("c_lm" in sim.init_cond, msg="c_lm not found in sim.init_cond")
        self.assertTrue("sign" in sim.init_cond.dims, msg="sign not found in sim.init_cond.dims")
        self.assertTrue("l" in sim.init_cond.dims, msg="l not found in sim.init_cond.dims")
        self.assertTrue("m" in sim.init_cond.dims, msg="m not found in sim.init_cond.dims")
        self.assertFalse("j2rp2" in sim.init_cond, msg="j2rp2 found in sim.init_cond")
        self.assertFalse("j4rp4" in sim.init_cond, msg="j4rp4 found in sim.init_cond")
        self.assertFalse("sign" in sim.init_cond.rh.dims, msg="sign found in sim.init_cond.rh.dims")
        self.assertFalse("l" in sim.init_cond.rh.dims, msg="l found in sim.init_cond.rh.dims")
        self.assertFalse("m" in sim.init_cond.rh.dims, msg="m found in sim.init_cond.rh.dims")

        sim.modify_body(name="Sun", j2rp2=1.0e-6, j4rp4=1.0e-8)
        self.assertEqual(sim.data.sel(name="Sun")["j2rp2"].values.item(), 1.0e-6)
        self.assertEqual(sim.data.sel(name="Sun")["j4rp4"].values.item(), 1.0e-8)

        self.assertFalse("c_lm" in sim.init_cond, msg="c_lm found in sim.init_cond")
        self.assertFalse("sign" in sim.init_cond.dims, msg="sign found in sim.init_cond.dims")
        self.assertFalse("l" in sim.init_cond.dims, msg="l found in sim.init_cond.dims")
        self.assertFalse("m" in sim.init_cond.dims, msg="m found in sim.init_cond.dims")
        self.assertTrue("j2rp2" in sim.init_cond, msg="j2rp2 not found in sim.init_cond")
        self.assertTrue("j4rp4" in sim.init_cond, msg="j4rp4 not found in sim.init_cond")

        return

    def test_remove_and_modify(self):
        sim = swiftest.Simulation(simdir=self.simdir)
        sim.add_solar_system_body(["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"])

        # Add 10 user-defined test particles.
        ntp = 10

        name_tp = [
            "TestParticle_01",
            "TestParticle_02",
            "TestParticle_03",
            "TestParticle_04",
            "TestParticle_05",
            "TestParticle_06",
            "TestParticle_07",
            "TestParticle_08",
            "TestParticle_09",
            "TestParticle_10",
        ]
        a_tp = rng.uniform(0.3, 1.5, ntp)
        e_tp = rng.uniform(0.0, 0.2, ntp)
        inc_tp = rng.uniform(0.0, 10, ntp)
        capom_tp = rng.uniform(0.0, 360.0, ntp)
        omega_tp = rng.uniform(0.0, 360.0, ntp)
        capm_tp = rng.uniform(0.0, 360.0, ntp)

        sim.add_body(name=name_tp, a=a_tp, e=e_tp, inc=inc_tp, capom=capom_tp, omega=omega_tp, capm=capm_tp)
        sim.modify_body(name="TestParticle_01", a=1.0, e=0.1, inc=0.0, capom=0.0, omega=0.0, capm=0.0)

    def test_read_with_dask(self):
        """
        Test that a swiftest Simulation data and collisions can be read in with dask.

        Adapted from test_fraggle.py.
        """
        # disruption_headon parameters
        names = ["Target", "Projectile"]
        pos_vectors = [np.array([1.0, -5.0e-05, 0.0]), np.array([1.0, 5.0e-05, 0.0])]
        vel_vectors = [np.array([0.00, 6.280005, 0.0]), np.array([0.00, 3.90, 0.0])]
        rot_vectors = [np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0])]
        body_Gmass = [1e-7, 1e-9]
        tstop = 2.0e-3
        nfrag_reduction = 10.0
        density = 3000 * swiftest.AU2M**3 / swiftest.MSun
        GU = swiftest.GMSun * swiftest.YR2S**2 / swiftest.AU2M**3
        body_radius = [((Gmass / GU) / (4.0 / 3.0 * np.pi * density)) ** (1.0 / 3.0) for Gmass in body_Gmass]

        sim = swiftest.Simulation(simdir=self.simdir, rotation=True, compute_conservation_values=True)
        sim.add_solar_system_body("Sun")
        sim.add_body(name=names, Gmass=body_Gmass, radius=body_radius, rh=pos_vectors, vh=vel_vectors, rot=rot_vectors)

        # Set fragmentation parameters
        minimum_fragment_gmass = 0.01 * body_Gmass[1]
        gmtiny = 0.50 * body_Gmass[1]
        sim.set_parameter(
            collision_model="fraggle",
            encounter_save="both",
            gmtiny=gmtiny,
            minimum_fragment_gmass=minimum_fragment_gmass,
            nfrag_reduction=nfrag_reduction,
        )
        sim.run(dt=tstop / 4, tstop=tstop, istep_out=1, dump_cadence=0)

        # read in the data with dask and then save a new set of initial conditions
        try:
            sim_read = swiftest.Simulation(simdir=self.simdir, read_param=True, read_data=True, read_collisions=True, dask=True)
            sim_read.save()
        except Exception as e:
            self.fail(f"Failed to read in data with dask: {e}")

        return

    def test_newrun_from_old(self):
        """
        Test that a new set of initial conditions can be extracted from an arbitrary output point of an old run.
        """
        run_args = {"tstop": 1.0, "dt": 0.01, "istep_out": 1, "dump_cadence": 0, "integrator": "whm"}
        # Build a fresh simulation
        sim = swiftest.Simulation(simdir=self.simdir)
        sim.add_solar_system_body(["Sun"])
        sim.add_body(a=1.0)
        sim.run(**run_args)

        # Build a new simulation from the half way point of the old one
        tmpdir2 = tempfile.TemporaryDirectory()
        simdir2 = tmpdir2.name
        sim.set_parameter(simdir=simdir2)
        sim.save(framenum=50)
        sim2 = swiftest.Simulation(simdir=simdir2, read_init_cond=True)
        tstart = sim2.init_cond.time.values[0]
        sim2.run(tstart=tstart, **run_args)

        # Now check if the final states of the two simulations are approximately the same:
        s1 = sim.data.isel(name=1, time=np.arange(50, 101))
        s2 = sim2.data.isel(name=1)
        self.assertTrue(np.allclose(s1.rh.values, s2.rh.values, rtol=1e-12), msg=f"Error in rh: {s1.rh.values - s2.rh.values}")
        self.assertTrue(np.allclose(s1.vh.values, s2.vh.values, rtol=1e-12), msg=f"Error in vh: {s1.vh.values - s2.vh.values}")

        tmpdir2.cleanup()
        return

    def test_verbose_flag(self):
        """
        Tests behavior of the verbose flag that is passed to various functions.
        """
        # Set up a default system (verbose should be set to True in this case)
        sim = swiftest.Simulation(simdir=self.simdir)

        # Capture the output
        captured_output = io.StringIO()
        sys.stdout = captured_output

        # Add the modern planets and the Sun using the JPL Horizons Database.
        sim.add_solar_system_body(["Sun", "Earth"])

        # Assert that something was printed when verbose=True
        self.assertTrue(captured_output.getvalue().strip() != "", "Verbose mode should print output, but nothing was printed.")

        # Reset output capture
        captured_output.truncate(0)
        captured_output.seek(0)

        # Test with verbose turned off
        sim = swiftest.Simulation(simdir=self.simdir, clean=True, verbose=False)
        sim.add_solar_system_body(["Sun", "Earth"])

        # Assert that nothing was printed when verbose=False
        self.assertEqual(captured_output.getvalue().strip(), "", "Verbose mode is off, but output was printed.")

        # Test that setting the verbose flag using set_parameter sets the verbose attribute of the Simulation object
        sim.set_parameter(verbose=True)
        self.assertTrue(sim.verbose, "Verbose flag should be True after setting it using set_parameter.")

        # Reset output capture
        captured_output.truncate(0)
        captured_output.seek(0)

        # Test that passing a different verbose flag to a function overrides the verbose attribute of the Simulation object
        sim.add_solar_system_body("Mercury", verbose=False)

        self.assertTrue(sim.verbose, "Verbose flag should still be set to True after calling a method with verbose=False.")
        sim.run(tstop=0.1, dt=0.01, verbose=False)
        # Assert that nothing was printed despite verbose=True in the object, because verbose=False was passed to the method
        self.assertEqual(captured_output.getvalue().strip(), "", "Method call with verbose=False should suppress output.")

        captured_output.truncate(0)
        captured_output.seek(0)

        sim.run(tstart=0, tstop=0.1)
        self.maxDiff = None

        self.assertTrue(captured_output.getvalue().strip() != "", "Verbose mode should print output, but nothing was printed.")

        # Clean up by resetting stdout
        sys.stdout = sys.__stdout__

    def test_init_cond_format_change(self):
        """
        Test that the user can switch between XV and EL input types even after bodies have been added.
        """

        class SimulationError(Exception):
            """Custom exception for simulation errors."""

            pass

        def _run_simulation(start_format, end_format):
            sim = swiftest.Simulation(simdir=self.simdir, init_cond_format=start_format)
            sim.clean(deep=True)
            sim.add_solar_system_body(["Sun", "Jupiter"])
            sim.set_parameter(init_cond_format=end_format, tstop=5, dt=5)
            sim.clean()
            sim.save()
            with chdir(sim.simdir):
                res = subprocess.run(["swiftest", "symba", "param.in", "quiet"], capture_output=True).stdout.decode()
            if "error" in res:
                raise SimulationError("Error in simulation")

        try:
            _run_simulation("XV", "EL")
        except Exception as e:
            self.fail(f"Failed XV->EL {e}")

        try:
            _run_simulation("EL", "XV")
        except Exception as e:
            self.fail(f"Failed EL->XV {e}")

    def test_symba_override_options(self):
        # Tests that the `rotation=False` and `compute_conservation_values=False` options are ignored when using the SyMBA integrator, and that
        # a warning will be issued if they are set to False.
        sim = swiftest.Simulation(simdir=self.simdir, integrator="symba")
        sim.add_solar_system_body(["Sun", "Mercury"])
        with self.assertWarns(UserWarning):
            sim.set_parameter(rotation=False)
        with self.assertWarns(UserWarning):
            sim.set_parameter(compute_conservation_values=False)
        self.assertTrue(sim.param["ROTATION"], "The `rotation` parameter should be set to True when using the SyMBA integrator.")
        self.assertTrue(
            sim.param["ENERGY"],
            "The `compute_conservation_values` parameter should be set to True when using the SyMBA integrator.",
        )

        # Test that manually overriding still throws an error on the Fortran side
        sim.param["ENERGY"] = False
        sim.save()
        with chdir(self.simdir):
            res = subprocess.run(["swiftest", "symba", "param.in", "quiet"], capture_output=True).stdout.decode()
            if "error" not in res:
                self.fail('Failed to throw error when setting ENERGY to "NO" param.in')

        sim.param["ENERGY"] = True
        sim.param["ROTATION"] = False
        sim.save()
        with chdir(self.simdir):
            res = subprocess.run(["swiftest", "symba", "param.in", "quiet"], capture_output=True).stdout.decode()
            if "error" not in res:
                self.fail('Failed to throw error when setting ROTATION to "NO" param.in')


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    unittest.main()
