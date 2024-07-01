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
import numpy as np
import pyshtools as pysh


class TestSwiftestRestart(unittest.TestCase):
       def setUp(self):
              # Initialize a target and surface for testing
              self.tmpdir=tempfile.TemporaryDirectory()
              self.simdir = self.tmpdir.name
        
       def tearDown(self):
              # Clean up temporary directory
              self.tmpdir.cleanup() 

       def test_restart_collision_snapshot(self):
              """
              Test that a restarted simulations runs correctly with collisions
              """

              sim = swiftest.Simulation(simdir = self.simdir, integrator = 'symba',
                                   dump_cadence = 1, tstep_out = 1 * 365.25, 
                                   tstop = 365.25 * 1, dt = 0.01,
                                   DU = 'km', TU = 'd', MU2KG = 1e15)
              sim.clean()
              sim.add_solar_system_body(['Mars', 'Phobos', 'Deimos'])

              # Add Mars c_lm from GMM3
              c_lm_data = pysh.datasets.Mars.GMM3(lmax = 6) # gravitational potential coefficients
              c_lm = c_lm_data.coeffs # 4pi normalized

              sim.modify_body(name = 'Mars', c_lm = c_lm)

              # hardcoded impactors to prevent use of my git repo

              names = np.array(['100', '101', '102', '103', '104', '105', '106', '107', '108',
                     '109', '598'])
              radius = np.array([0.35705834, 0.07421946, 0.14916194, 0.13296731, 0.12915766,
                     0.41544253, 0.46550525, 0.17445848, 0.41888955, 0.45045171,
                     0.49215104 ])
              mass = 10 * np.array([2.86020662e-04, 2.56881727e-06, 2.08522999e-05, 1.47711587e-05,
                     1.35375579e-05, 4.50518599e-04, 6.33802015e-04, 3.33623100e-05,
                     4.61826081e-04, 5.74281193e-04, 7.48988414e-04])
              rh = np.array([[-22124.94631599,  -7600.61943767,  -1953.85576394],
                     [-22099.04261316,  -7594.5053911 ,  -1955.35554083],
                     [-22122.04312455,  -7615.95767634,  -1955.17148458],
                     [-22083.30057286,  -7614.06013637,  -1953.85760781],
                     [-22107.9145528 ,  -7592.30630884,  -1950.69068123],
                     [-22124.83908881,  -7600.49807585,  -1961.07379115],
                     [-22091.70625466,  -7601.39044401,  -1943.88804141],
                     [-22091.96202388,  -7598.98000554,  -1952.81119451],
                     [-22098.41518385,  -7590.14851984,  -1954.86866938],
                     [-22117.17018901,  -7620.37817046,  -1965.44576821],
                     [-22104.19743262,  -7621.76082186,  -1942.37784931]])
              vh =  0.503 * np.array([[  36796.54983816, -110084.06529807,    2987.13829165],
                     [  38083.48622951, -110102.85015476,    2921.36450056],
                     [  36869.30640377, -110484.65220524,    2935.57333344],
                     [  38278.23679524, -110732.77923   ,    2990.69525294],
                     [  37624.67760602, -109538.53678746,    3189.04717925],
                     [  36792.03452032, -109813.24269459,    2613.69041106],
                     [  38540.20902575, -110535.5996834 ,    3651.11554834],
                     [  38768.73186052, -109936.55975385,    3071.4989985 ],
                     [  38141.5822543 , -109706.74112105,    2931.3008305 ],
                     [  36896.44746876, -111001.2491618 ,    2413.7897524 ],
                     [  37671.8566089 , -110761.29152044,    3217.24051514]])

              sim.add_body(name = names, radius = radius, mass = mass, rh = rh, vh = vh)

              try:
                 sim.run() 
              except Exception as e:
                 self.fail(f'Failed initial run with Exception: {e}')

              # restart run

              sim_restart = swiftest.Simulation(simdir = self.simdir, 
                                          read_param = True, read_collisions=False, read_data = True,
                                          param_file = 'param.restart.in')

              sim_restart.set_parameter(tstop = 365.25 * 10)

              try:
                     sim_restart.run()
              except Exception as e:
                     self.fail(f'Failed restart with Exception: {e}')
                     
              return

if __name__ == '__main__':
    os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
    unittest.main()