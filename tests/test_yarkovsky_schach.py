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
        r_p_SI_name = ['0.25m', '1m']
        rho_p_SI = 1500.0 # disk particle density
        Mars_k2 = 0.164
        Mars_Q = 99.5

        # mass distribution parameters
        m_dist_type = 'arbitrary'
        m_dist_sigma0 = 1.05e5 * 2 / 115.7181 # 2 M_Phobos
        m_dist_alpha = -8.0
        m_dist_nbins = 256
        m_dist_r_inner = 0.99
        m_dist_r_outer = 8

        # pre-computed surface mass density values right before a gap appears (in code units)
        pre_computed_sigma = {'0_obl' : np.array([8.01142413e-16, 1.20993033e-10, 1.02943128e-10, 1.88278094e-12, 2.14993115e-12, 4.17531874e-12, 2.79240120e-12, 6.32563069e-13,
        1.16105228e-13, 3.81390523e-12, 5.16478425e-12, 3.32558324e-12, 6.10291450e-13, 1.88466380e-13, 2.61573271e-12, 6.00607115e-12,
        5.83756825e-12, 3.26434078e-12, 5.80030960e-13, 1.96203569e-13, 2.41295569e-12, 7.28874920e-12, 6.03780778e-12, 3.31979966e-12,
        8.27890256e-13, 5.14129925e-14, 1.11327139e-12, 8.71095377e-12, 9.55756344e-12, 7.58087527e-12, 4.48419144e-12, 1.56363451e-12,
        2.19295660e-13, 3.63113058e-14, 7.33655597e-12, 1.54006292e-11, 1.31939610e-11, 8.98925446e-12, 4.70170941e-12, 1.23284738e-12,
        1.32276481e-13, 9.24827414e-13, 1.56018978e-11, 2.98521187e-11, 2.48423908e-11, 1.90583367e-11, 1.23879996e-11, 6.74461689e-12,
        2.25931428e-12, 4.35384396e-13, 9.30025879e-13, 1.03409199e-11, 3.71224502e-10, 1.08549917e-09, 1.20153787e-09, 1.34233872e-09,
        1.52657180e-09, 1.76088362e-09, 2.03041832e-09, 2.30268127e-09, 2.55427543e-09, 2.77689447e-09, 2.96758245e-09, 3.12444500e-09,
        3.24684324e-09, 3.33451620e-09, 3.38820917e-09, 3.40968457e-09, 3.40062760e-09, 3.36320245e-09, 3.29967730e-09, 3.21259715e-09,
        3.10436600e-09, 2.97793496e-09, 2.83577071e-09, 2.68181161e-09, 2.51848700e-09, 2.35223494e-09, 2.18349588e-09, 2.02312404e-09,
        1.86262923e-09, 1.71963273e-09, 1.56731223e-09, 1.44015767e-09, 1.29124235e-09, 1.18778894e-09, 1.05600054e-09, 9.83270523e-10,
        8.71215883e-10, 8.22735761e-10, 7.22224887e-10, 6.94135141e-10, 6.00672067e-10, 5.89165145e-10, 5.00914754e-10, 5.04000962e-10,
        4.17083868e-10, 4.34461506e-10, 3.47502598e-10, 3.77012721e-10, 2.88517287e-10, 3.30268190e-10, 2.37877651e-10, 2.91673435e-10,
        1.94572755e-10, 2.60260529e-10, 1.56551320e-10, 2.34562167e-10, 1.23757935e-10, 2.14968877e-10, 9.34711995e-11, 1.97297280e-10,
        6.72006142e-11, 1.84797188e-10, 4.30739404e-11, 1.75707627e-10, 2.08812011e-11, 1.69556142e-10, 5.05592253e-14, 1.59591608e-16,
        2.07240335e-16, 2.72354907e-16, 2.59578220e-16, 2.45101179e-16, 2.32175716e-16, 2.21058776e-16, 2.09703030e-16, 1.99064727e-16,
        1.86319648e-16, 1.77004306e-16, 1.70290293e-16, 1.62120069e-16, 1.50259847e-16, 1.39436369e-16, 1.33215681e-16, 1.29798529e-16,
        1.23585880e-16, 1.12989134e-16, 1.03115997e-16, 9.77183602e-17, 9.72814645e-17, 9.45291807e-17, 8.76420834e-17, 7.75906852e-17,
        6.99637020e-17, 6.81328901e-17, 6.97277940e-17, 6.84224573e-17, 6.37028555e-17, 5.44192216e-17, 4.70767967e-17, 4.38462699e-17,
        4.60441189e-17, 4.86220646e-17, 4.93787662e-17, 4.47918974e-17, 3.76407211e-17, 2.95764305e-17, 2.50326150e-17, 2.43568362e-17,
        2.74979966e-17, 3.22845319e-17, 3.55327127e-17, 3.57452383e-17, 3.21039194e-17, 2.55027023e-17, 1.74372406e-17, 9.51496981e-18,
        3.40942992e-18, 5.99154592e-19, 5.90461545e-20, 9.88964196e-18, 2.59549027e-16, 4.72728150e-16, 3.95029752e-16, 2.25715928e-16,
        1.31833537e-16, 8.60886566e-16, 1.16010179e-14, 2.60430364e-13, 1.39077452e-12, 1.45336960e-12, 1.40209830e-12, 1.36667671e-12,
        1.30359481e-12, 1.25919231e-12, 1.19160404e-12, 1.14593457e-12, 1.07745996e-12, 1.03600304e-12, 9.69155532e-13, 9.33192355e-13,
        8.70552259e-13, 8.39146763e-13, 7.81671209e-13, 7.55038766e-13, 7.01176370e-13, 6.80434953e-13, 6.28943674e-13, 6.13350989e-13,
        5.64984506e-13, 5.52925463e-13, 5.07623241e-13, 4.99487900e-13, 4.55595821e-13, 4.51878049e-13, 4.09314455e-13, 4.08689894e-13,
        3.68032386e-13, 3.70330620e-13, 3.30356273e-13, 3.36305035e-13, 2.96691628e-13, 3.05325000e-13, 2.66629108e-13, 2.77844183e-13,
        2.39034547e-13, 2.53414572e-13, 2.14426528e-13, 2.31127826e-13, 1.92214036e-13, 2.11544786e-13, 1.71791730e-13, 1.93903601e-13,
        1.53646960e-13, 1.77992214e-13, 1.36900948e-13, 1.64067337e-13, 1.21731167e-13, 1.51329402e-13, 1.08006575e-13, 1.40207068e-13,
        9.52460454e-14, 1.30204091e-13, 8.38086914e-14, 1.21315815e-13, 7.30772454e-14, 1.13544208e-13, 6.33635236e-14, 1.06546643e-13,
        5.42815889e-14, 1.00563028e-13, 4.59103978e-14, 9.51872140e-14, 3.81139129e-14, 9.06633275e-14, 3.08062737e-14, 8.66712527e-14,
        2.39847517e-14, 8.33730700e-14, 1.75114993e-14, 8.05628920e-14, 1.14137342e-14, 7.83177595e-14, 5.57705788e-15, 7.65298759e-14]),
        '25_obl': np.array([8.00609932e-16, 1.20952802e-10, 1.01707231e-10, 1.92584619e-12, 3.18070645e-12, 4.12625118e-12, 2.46473509e-12, 2.76365048e-13,
       1.64157126e-12, 4.44724418e-12, 5.02458661e-12, 2.86031523e-12, 5.44685868e-13, 7.15965139e-14, 3.12658589e-12, 7.13293532e-12,
       4.89929459e-12, 1.44666266e-12, 1.10874185e-14, 1.51081612e-12, 6.30800331e-12, 7.16546567e-12, 5.28855894e-12, 2.23654949e-12,
       3.09482132e-13, 1.16177733e-13, 4.92496689e-12, 1.03730831e-11, 7.97915514e-12, 4.22940276e-12, 1.02817280e-12, 1.51686291e-13,
       1.11576417e-12, 1.41910149e-11, 1.52372116e-11, 1.15211444e-11, 6.33187323e-12, 1.92066060e-12, 9.12421130e-14, 1.11019243e-12,
       1.83297545e-11, 3.40915919e-11, 2.73329299e-11, 1.82955111e-11, 9.67668613e-12, 2.96146212e-12, 3.64635299e-13, 3.63879211e-13,
       3.86354608e-11, 7.89740819e-10, 1.34578578e-09, 1.54028397e-09, 1.84639918e-09, 2.24819489e-09, 2.68527860e-09, 3.09827839e-09,
       3.46971895e-09, 3.79287779e-09, 4.06239462e-09, 4.27494376e-09, 4.42940907e-09, 4.52647291e-09, 4.56811355e-09, 4.55718379e-09,
       4.49711760e-09, 4.39173000e-09, 4.24517010e-09, 4.06189691e-09, 3.84695426e-09, 3.60588275e-09, 3.34610764e-09, 3.07469431e-09,
       2.80474329e-09, 2.53849727e-09, 2.30217174e-09, 2.06071010e-09, 1.90934078e-09, 1.65357657e-09, 1.13807426e-09, 1.04508156e-09,
       9.51756169e-10, 8.65604141e-10, 7.86368654e-10, 7.11381854e-10, 6.42366014e-10, 5.80848905e-10, 5.24046467e-10, 4.70974407e-10,
       4.23679387e-10, 3.82192230e-10, 3.43027123e-10, 3.07097813e-10, 2.74679210e-10, 2.48533324e-10, 2.21192140e-10, 1.99389997e-10,
       1.74044009e-10, 1.61098609e-10, 1.39298983e-10, 1.31250716e-10, 1.08344993e-10, 1.03755549e-10, 8.49163253e-11, 8.50612970e-11,
       6.84010528e-11, 6.96290909e-11, 5.08057740e-11, 5.50690285e-11, 3.61814498e-11, 4.52972187e-11, 2.74645820e-11, 4.03129527e-11,
       2.06927779e-11, 3.70298154e-11, 1.39611548e-11, 3.39707307e-11, 6.97993399e-12, 3.22245513e-11, 3.76372774e-15, 1.36663558e-19,
       5.23753531e-19, 4.71800248e-18, 4.60807686e-18, 4.24728042e-18, 3.17737430e-18, 2.48705186e-18, 2.87162176e-18, 3.57251630e-18,
       3.44386303e-18, 2.16799266e-18, 7.24339408e-19, 2.97569093e-20, 9.68681284e-19, 1.92017577e-17, 7.45916529e-16, 3.18817354e-16,
       3.05963476e-16, 1.92491264e-15, 1.54370419e-14, 1.11401985e-12, 2.21798389e-12, 2.05184950e-12, 1.89403637e-12, 1.74489352e-12,
       1.60941313e-12, 1.48815906e-12, 1.37914138e-12, 1.28025794e-12, 1.18986884e-12, 1.10681793e-12, 1.03027243e-12, 9.59596838e-13,
       8.94262100e-13, 8.33817276e-13, 7.77850471e-13, 7.25998949e-13, 6.77937682e-13, 6.33351862e-13, 5.91970728e-13, 5.53542308e-13,
       5.17836342e-13, 4.84642228e-13, 4.53767177e-13, 4.25034562e-13, 3.98282436e-13, 3.73362199e-13, 3.50137397e-13, 3.28482637e-13,
       3.08282614e-13, 2.89431226e-13, 2.71830771e-13, 2.55391231e-13, 2.40029610e-13, 2.25669345e-13, 2.12239761e-13, 1.99675583e-13,
       1.87916488e-13, 1.76906700e-13, 1.66594618e-13, 1.56932480e-13, 1.47876052e-13, 1.39384348e-13, 1.31419376e-13, 1.23945895e-13,
       1.16931208e-13, 1.10344959e-13, 1.04158954e-13, 9.83469982e-14, 9.28847381e-14, 8.77495268e-14, 8.29202939e-14, 7.83774275e-14,
       7.41026661e-14, 7.00789989e-14, 6.62905736e-14, 6.27226120e-14, 5.93613320e-14, 5.61938753e-14, 5.32082418e-14, 5.03932273e-14,
       4.77383679e-14, 4.52338869e-14, 4.28706471e-14, 4.06401056e-14, 3.85342729e-14, 3.65456741e-14, 3.46673138e-14, 3.28926432e-14,
       3.12155296e-14, 2.96302281e-14, 2.81313556e-14, 2.67138661e-14, 2.53730288e-14, 2.41044064e-14, 2.29038362e-14, 2.17674116e-14,
       2.06914658e-14, 1.96725555e-14, 1.87074467e-14, 1.77931012e-14, 1.69266639e-14, 1.61054509e-14, 1.53269389e-14, 1.45887546e-14,
       1.38886654e-14, 1.32245706e-14, 1.25944927e-14, 1.19965706e-14, 1.14290504e-14, 1.08902846e-14, 1.03787047e-14, 9.89288922e-15,
       9.43130911e-15, 8.99313358e-15, 8.57568245e-15, 8.18207430e-15, 7.79917678e-15, 7.45823987e-15, 7.07608348e-15, 6.84991452e-15,
       6.32870281e-15, 6.47162490e-15, 5.36599011e-15, 6.54696131e-15, 3.98019961e-15, 7.09780771e-15, 2.48828738e-15, 7.37466525e-15,
       1.83963738e-15, 6.75292633e-15, 1.74716302e-15, 6.40864112e-15, 1.00479330e-15, 6.49730317e-15, 4.96949824e-16, 6.28606859e-15]) }

        # YS parameters
        albedo = 0.07
        emissivity = 0.9
        thermal_inertia_SI = 1500.0
        rot_k = 0.25
        a_pl = 67248.67 # Mars a in R_Mars
        obliquities = [0, 25] # degrees
        obliquities_names = ["0_obl", "25_obl"]

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
            for i, obliquity in enumerate(obliquities):

                sigma = pre_computed_sigma[obliquities_names[i]]
                r_p = r_p_SI * sim.M2DU  # disk particle size
                rho_p = rho_p_SI * sim.KG2MU / sim.M2DU**3  # disk particle density
                m_p = 4.0 / 3.0 * np.pi * r_p**3 * rho_p  # disk particle mass
                
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
                                "type": "arbitrary",
                                "sigma": sigma,
                                "r_inner": m_dist_r_inner,
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

                print(f'r_p = {r_p_SI}, obliquity = {obliquities_names[i]} ')
                # compare with pre-tested results that ring surface mass density --> 0 at r = 1.94 R_Mars (ringbin = 55) 32 -> 1.5 R_Mars
                sigma = sim.ring.sel(time = 2e6, ringbin = 50).sigma.values
                self.assertLessEqual(sigma, 5e-10)

                sim.clean()
        return


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    unittest.main()
