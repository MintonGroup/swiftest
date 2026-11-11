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

    def test_yarkovsky_schach(self):
        """
        Tests that the Yarkovsky-Schach effect is handled correctly in ringmoons.
        """
        # ring and particle parameters
        rho_p_SI = 1500.0 # disk particle density
        Mars_k2 = 0.164
        Mars_Q = 99.5

        # mass distribution parameters
        m_dist_type = 'powerlaw'
        m_dist_sigma0 = 1.05e7 #1.05e5 * 2 / 115.7181 # 2 M_Phobos # old num = 1.05e7
        m_dist_alpha = -8.0
        m_dist_nbins = 256
        m_dist_r_outer = 8
        ring_r_max = 6.4 # cb units

        # YS parameters
        albedo = 0.07
        emissivity = 0.9
        rot_k = 0.25
        thermal_inertia_scaled = 1500.0 * (sim.TU2S**(5.0/2)) / sim.MU2KG
        a_pl = 67248.67 # Mars a in R_Mars
        # ring particle properties
        obliquities = [0, 25, 50] # degrees
        r_p_SIs = [1e0, 5e0, 10e0]
        r_p = r_p_SI * sim.M2DU  # disk particle size
        rho_p = rho_p_SI * sim.KG2MU / sim.M2DU**3  # disk particle density
        m_p = 4.0 / 3.0 * np.pi * r_p**3 * rho_p  # disk particle mass

        # set up the simulation object
        cb = swiftest.get_solar_system_body("Mars")
        r_cb = cb["radius"]
        m_cb = cb["mass"]

        sim = swiftest.Simulation(simdir = self.simdir, 
                                  integrator = 'ringmoons', 
                                  MU2KG = m_cb, DU2M = r_cb, TU = "yr", 
                                  tstop = 2e6, 
                                  dt=100,
                                  tstep_out = 2e6, 
                                  dump_cadence=0,
                                  yarkovsky_schach = True)
        sim.clean()
        sim.add_solar_system_body("Mars", align_to_central_body_rotation=True)
        sim.modify_body(name = "Mars", k2 = Mars_k2, Q = Mars_Q)
        
        sim.add_ring(
                    r_p=r_p,
                    m_p=m_p,
                    mass_distribution={
                        "type": m_dist_type,
                        "sigma0": m_dist_sigma0 * sim.KG2MU / sim.M2DU**2,
                        "alpha": m_dist_alpha,
                        "nbins": m_dist_nbins,
                        "r_outer": m_dist_r_outer
                    },
                    albedo = albedo,
                    emissivity=emissivity,
                    rot_k = rot_k,
                    gamma = thermal_inertia_scaled,
                    a_pl = a_pl,
                    obliquity = obliquity
                    # Y_21 = None
                    # delta = None
                )

        # check that data inputs are handled correctly
        # if one value is missing, an error should be thrown. If multiple values are missing, this test should extend to cover that.

        # missing albedo

        with self.assertRaises(ValueError):
            sim.add_ring(
                        r_p=r_p,
                        m_p=m_p,
                        mass_distribution={
                            "type": m_dist_type,
                            "sigma0": m_dist_sigma0 * sim.KG2MU / sim.M2DU**2,
                            "alpha": m_dist_alpha,
                            "nbins": m_dist_nbins,
                            "r_outer": m_dist_r_outer
                        },
                        emissivity=emissivity,
                        rot_k = rot_k,
                        gamma = thermal_inertia_scaled,
                        a_pl = a_pl,
                        obliquity = obliquity
                        )

        # missing emissivity
        with self.assertRaises(ValueError):
            sim.add_ring(
                        r_p=r_p,
                        m_p=m_p,
                        mass_distribution={
                            "type": m_dist_type,
                            "sigma0": m_dist_sigma0 * sim.KG2MU / sim.M2DU**2,
                            "alpha": m_dist_alpha,
                            "nbins": m_dist_nbins,
                            "r_outer": m_dist_r_outer
                        },
                        albedo = albedo,
                        rot_k = rot_k,
                        gamma = thermal_inertia_scaled,
                        a_pl = a_pl,
                        obliquity = obliquity
                        )

        # missing rotational k constant
        with self.assertRaises(ValueError):
            sim.add_ring(
                        r_p=r_p,
                        m_p=m_p,
                        mass_distribution={
                            "type": m_dist_type,
                            "sigma0": m_dist_sigma0 * sim.KG2MU / sim.M2DU**2,
                            "alpha": m_dist_alpha,
                            "nbins": m_dist_nbins,
                            "r_outer": m_dist_r_outer
                        },
                        albedo = albedo,
                        emissivity=emissivity,
                        gamma = thermal_inertia_scaled,
                        a_pl = a_pl,
                        obliquity = obliquity
                        )

        # missing thermal inertia (gamma)
        with self.assertRaises(ValueError):
            sim.add_ring(
                        r_p=r_p,
                        m_p=m_p,
                        mass_distribution={
                            "type": m_dist_type,
                            "sigma0": m_dist_sigma0 * sim.KG2MU / sim.M2DU**2,
                            "alpha": m_dist_alpha,
                            "nbins": m_dist_nbins,
                            "r_outer": m_dist_r_outer
                        },
                        albedo = albedo,
                        emissivity=emissivity,
                        rot_k = rot_k,
                        a_pl = a_pl,
                        obliquity = obliquity
                        )
            
        # missing planet's semi-major axis with respect to the Sun
        with self.assertRaises(ValueError):
            sim.add_ring(
                        r_p=r_p,
                        m_p=m_p,
                        mass_distribution={
                            "type": m_dist_type,
                            "sigma0": m_dist_sigma0 * sim.KG2MU / sim.M2DU**2,
                            "alpha": m_dist_alpha,
                            "nbins": m_dist_nbins,
                            "r_outer": m_dist_r_outer
                        },
                        albedo = albedo,
                        emissivity=emissivity,
                        rot_k = rot_k,
                        gamma = thermal_inertia_scaled,
                        obliquity = obliquity
                        )

        # missing planet obliquity
        with self.assertRaises(ValueError):
            sim.add_ring(
                        r_p=r_p,
                        m_p=m_p,
                        mass_distribution={
                            "type": m_dist_type,
                            "sigma0": m_dist_sigma0 * sim.KG2MU / sim.M2DU**2,
                            "alpha": m_dist_alpha,
                            "nbins": m_dist_nbins,
                            "r_outer": m_dist_r_outer
                        },
                        albedo = albedo,
                        emissivity=emissivity,
                        rot_k = rot_k,
                        a_pl = a_pl,
                        gamma = thermal_inertia_scaled,
                        obliquity = obliquity
                        )

        sim.clean()

        # check that the ring is pushed away from the planet

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
                self.assertAlmostEqual(da.values, da_pre_tested[f"{radius}"][f"{obliquity}"], delta=delta[f"{radius}"])
                sim.clean()
        return


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    unittest.main()
