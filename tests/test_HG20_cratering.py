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
class TestCollisions(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir=tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name
        
    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup() 

        
    def test_HG20_cratering(self):
        """
        Tests that cratering collisions are handled correctly in the HG20 model.

        HG20: Hyodo, R., Genda, H., 2020. Escape and Accretion by Cratering Impacts: Formulation of Scaling Relations for 
                High-speed Ejecta. ApJ 898 30. https://doi.org/10.3847/1538-4357/ab9897
        """
        
        impact_angles = np.array([30, 45, 60, 75, 90]) # degrees
        impact_velocities = np.array([1.0, 3.0, 5.0, 10.0, 15.0]) # v_escape
        impact_mass_changes = {impact_angle: {impact_velocity: None for impact_velocity in impact_velocities} for impact_angle in impact_angles}
        simulation_derived_mass_changes = [1.0, 1.0, -0.3938, -1.5315, -3.0704, 
                                           1.0, 1.0, -0.3082, -3.0934, -6.4454,
                                           1.0, 1.0, 1.0, -3.8059, -8.9207,
                                           1.0, 1.0, 1.0, -3.0412, -10.6972,
                                           1.0, 1.0, 1.0, -2.71889, -11.4195]

        i = 0
        for angle in impact_angles:
            for velocity in impact_velocities:
                impact_mass_changes[angle][velocity] = simulation_derived_mass_changes[i]
                i += 1

        ###################################################################
        # Setup and run through each set of impact angles and velocities

        for angle in impact_angles[:]:
            print(f'\nimpact angle = {angle} deg')
            for velocity in impact_velocities[:]:
                expected_mass_change = impact_mass_changes[angle][velocity]

                sim = swiftest.Simulation(simdir = self.simdir, tstop = 0.01, dt = 0.005, 
                                    collision_model = "FRAGGLE", MU2KG = 1e15, DU = 'km', TU = 'day')
                
                sim.add_solar_system_body(['Mars', 'Phobos', 'Deimos'])

                sim.modify_body(name = 'Deimos', radius = 6.203) # JPL Horizons only gives the longest axis
                deimos_ds = sim.data.sel(name = 'Deimos')
                deimos_rh = deimos_ds.rh.values[0]
                deimos_vh = deimos_ds.vh.values[0]
                deimos_radius = deimos_ds.radius.values[0]
                deimos_v_esc = np.sqrt(2 * deimos_ds.Gmass.values[0] / deimos_radius)

                # set up impactor "behind" deimos => r_hat = - vh_deimos_hat
                # we want the impact in the first time step
                r_mag = 1.07 * deimos_radius # 1 v_esc * 0.005 days ~ 0.38 deimos radii
                r_hat = - deimos_vh / np.linalg.norm(deimos_vh)
                r = r_mag * r_hat

                # impactor velocity
                # rotate r_hat by angle
                # for simplicity keep Z component the same
                v_mag = velocity * deimos_v_esc

                rot_angle = np.radians(90.0 - angle) # impact angle = 90 - angle between r and v
                v_hat_0 = -r_hat # v_hat at 0 degree rotation
                x_tmp = v_hat_0[0]
                y_tmp = v_hat_0[1]
                    
                x = x_tmp * np.cos(rot_angle) + y_tmp * np.sin(rot_angle)
                y = x_tmp * -1.0 * np.sin(rot_angle) + y_tmp * np.cos(rot_angle)

                v_hat = np.array([x, y, v_hat_0[2]])
                v = v_mag * v_hat

                impactor_radius = 0.1 # km
                impactor_mass = deimos_ds.mass.values * impactor_radius**3 / deimos_radius**3
                impactor_rh = r + deimos_rh
                impactor_vh = v + deimos_vh

                sim.add_body(name = 'Impactor', radius = impactor_radius, mass = impactor_mass, rh = impactor_rh, vh = impactor_vh)
                
                sim.run()
                
                ##########################################################
                # Extract impact data and check

                col = sim.collisions.sel(collision_id = 1, stage = 'before')

                names = col['name']
                if 'Deimos' not in names:
                    self.fail(f'Body Deimos not found in collision data')
                
                # Extract body velocity
                col_body = col.sel(collision_body = 1)
                v_body = col_body['vh'].values
                r_body = col_body['rh'].values
                Gmass_before = col_body['Gmass'].values

                Gmass_after = sim.collisions.sel(collision_id = 1, stage = 'after', collision_body = 1)['Gmass'].values
                Gmass_change = Gmass_after - Gmass_before

                # impact velocity

                # Extract impactor velocity and position
                impactor_gmass = col.sel(collision_body = 2)['Gmass'].values
                v_impactor = col.sel(collision_body = 2)['vh'].values - v_body
                r_impactor = col.sel(collision_body = 2)['rh'].values - r_body

                # Calculate impactor velocity and position magnitude
                vmag = np.linalg.norm(v_impactor)
                rmag = np.linalg.norm(r_impactor)

                # impact angle 
                impact_angle = np.dot(r_impactor, v_impactor) / vmag / rmag
                impact_angle = np.abs(90.0 - np.rad2deg(np.arccos(impact_angle))) # degrees 

                # append to output arrays
                impact_velocity = np.array(vmag/ deimos_v_esc)
                impact_mass_change = np.array(Gmass_change / impactor_gmass) # m_impactor units
                impact_angle = np.array(impact_angle) 

                print(f'impact angle = {impact_angle} deg\t expected = {angle} deg')
                print(f'impact velocity = {impact_velocity} v_esc\t expected = {velocity} v_esc')
                print(f'impact mass change = {impact_mass_change} m_impactor\t expected = {expected_mass_change} m_impactor')

                self.assertGreater(impact_angle, 0, msg = 'Impact angle should be positive')
                self.assertAlmostEqual(impact_angle, angle, delta = 7, msg = 'Impact angle not as expected')
                self.assertAlmostEqual(impact_velocity, velocity, delta = 0.1 * velocity, msg = 'Impact velocity not as expected')
                self.assertAlmostEqual(impact_mass_change, expected_mass_change, places=1, msg = 'Mass change due to impact not as expected')

        return
            
if __name__ == '__main__':
    os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
    unittest.main()



