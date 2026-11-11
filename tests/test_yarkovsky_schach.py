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
        r_p_SIs = [25e-2, 1e0]
        rho_p_SI = 1500.0 # disk particle density
        Mars_k2 = 0.164
        Mars_Q = 99.5

        # mass distribution parameters
        m_dist_type = 'powerlaw'
        m_dist_sigma0 = 1.05e7 #1.05e5 * 2 / 115.7181 # 2 M_Phobos # old num = 1.05e7
        m_dist_alpha = -8.0
        m_dist_nbins = 256
        m_dist_r_outer = 8

        # YS parameters
        albedo = 0.07
        emissivity = 0.9
        thermal_inertia_SI = 1500.0
        rot_k = 0.25
        a_pl = 67248.67 # Mars a in R_Mars
        obliquities = [0, 25] # degrees

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

        obliquity = obliquities[0]
        r_p = r_p_SIs[0] * sim.M2DU  # disk particle size
        rho_p = rho_p_SI * sim.KG2MU / sim.M2DU**3  # disk particle density
        m_p = 4.0 / 3.0 * np.pi * r_p**3 * rho_p  # disk particle mass
        thermal_inertia_scaled = thermal_inertia_SI * (sim.TU2S**(5.0/2)) / sim.MU2KG

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
                        gamma = thermal_inertia_scaled
                        )

        sim.clean()

        # check that the ring is pushed away from the planet

        for r_p_SI in r_p_SIs:
            for obliquity in obliquities:

                r_p = r_p_SI * sim.M2DU  # disk particle size
                rho_p = rho_p_SI * sim.KG2MU / sim.M2DU**3  # disk particle density
                m_p = 4.0 / 3.0 * np.pi * r_p**3 * rho_p  # disk particle mass
                
                sim = swiftest.Simulation(simdir = self.simdir, 
                                            integrator = 'ringmoons', 
                                            MU2KG = m_cb, DU2M = r_cb, TU = "yr", 
                                            tstop = 5e8, 
                                            dt=100,
                                            tstep_out = 1e8, 
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
                        )

                # # check that delta = 180 degrees at a = R_pl
                # self.assertAlmostEqual(sim.ring.sel(time = 0, ringbin = 0).delta.values, 180.0, delta = 1.0)

                sim.run()

                # compare with pre-tested results that ring surface mass density --> 0 at r = 1.5 R_Mars (ringbin = 32)
                sigma = sim.ring.sel(time = 5e8, ringbin = 32).sigma 
                self.assertLessEqual(sigma, 1e-10, delta = 1e-10)

                sim.clean()
        return


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    unittest.main()
