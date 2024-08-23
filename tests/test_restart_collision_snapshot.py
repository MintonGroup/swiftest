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

# Colliion initial conditions taken from Fragmentation_Movie.py in the example directory

collision_type = ["disruption_headon", "disruption_off_axis", "supercatastrophic_headon", "supercatastrophic_off_axis","hitandrun_disrupt", "hitandrun_pure", "merge", "merge_spinner"]

# ----------------------------------------------------------------------------------------------------------------------
# To increase the number of bodies generated in each collision type, decrease the value of the corresponding nfrag_reduction number
# ----------------------------------------------------------------------------------------------------------------------
nfrag_reduction = {"disruption_headon" : 10.0,
        "disruption_off_axis"         : 10.0,
        "supercatastrophic_headon"    : 10.0,
        "supercatastrophic_off_axis"  : 10.0,
        "hitandrun_disrupt"           : 10.0,
        "hitandrun_pure"              : 1.0,
        "merge"                       : 1.0,
        "merge_spinner"               : 1.0,
        }


# These initial conditions were generated by trial and error
names = ["Target","Projectile"]
pos_vectors = {"disruption_headon"         : [np.array([1.0, -5.0e-05, 0.0]),
                                            np.array([1.0,  5.0e-05 ,0.0])],
            "disruption_off_axis"        : [np.array([1.0, -5.0e-05, 0.0]),
                                            np.array([1.0,  5.0e-05 ,0.0])],
            "supercatastrophic_headon":   [np.array([1.0, -5.0e-05, 0.0]),
                                            np.array([1.0,  5.0e-05, 0.0])],
            "supercatastrophic_off_axis": [np.array([1.0, -5.0e-05, 0.0]),
                                            np.array([1.0,  5.0e-05, 0.0])],
            "hitandrun_disrupt"         : [np.array([1.0, -4.2e-05, 0.0]),
                                            np.array([1.0,  4.2e-05, 0.0])],
            "hitandrun_pure"            : [np.array([1.0, -4.2e-05, 0.0]),
                                            np.array([1.0,  4.2e-05, 0.0])],
            "merge"                      : [np.array([1.0, -5.0e-05, 0.0]),
                                            np.array([1.0,  5.0e-05 ,0.0])],
            "merge_spinner"               : [np.array([1.0, -5.0e-05, 0.0]),
                                            np.array([1.0,  5.0e-05 ,0.0])]
            }

vel_vectors = {"disruption_headon"         : [np.array([ 0.00,  6.280005, 0.0]),
                                            np.array([ 0.00,  3.90,     0.0])],
            "disruption_off_axis"       : [np.array([ 0.00,  6.280005, 0.0]),
                                            np.array([ 0.05,  3.90,     0.0])],
            "supercatastrophic_headon":   [np.array([ 0.00,  6.28,     0.0]),
                                            np.array([ 0.00,  4.28,     0.0])],
            "supercatastrophic_off_axis": [np.array([ 0.00,  6.28,     0.0]),
                                            np.array([ 0.05,  4.28,     0.0])],
            "hitandrun_disrupt"         : [np.array([ 0.00,  6.28,     0.0]),
                                            np.array([-1.45, -6.28,     0.0])],
            "hitandrun_pure"            : [np.array([ 0.00,  6.28,     0.0]),
                                            np.array([-1.52, -6.28,     0.0])],
            "merge"                     : [np.array([ 0.04,  6.28,     0.0]),
                                            np.array([ 0.05,  6.18,     0.0])],
            "merge_spinner"             : [np.array([ 0.04,  6.28,     0.0]),
                                            np.array([ 0.05,  6.18,     0.0])]
            }

rot_vectors = {"disruption_headon"         : [np.array([0.0, 0.0, 1.0e5]),
                                            np.array([0.0, 0.0, -5e5])],
            "disruption_off_axis":        [np.array([0.0, 0.0, 2.0e5]),
                                            np.array([0.0, 0.0, -1.0e5])],
            "supercatastrophic_headon":   [np.array([0.0, 0.0, 1.0e5]),
                                            np.array([0.0, 0.0, -5.0e5])],
            "supercatastrophic_off_axis": [np.array([0.0, 0.0, 1.0e5]),
                                            np.array([0.0, 0.0, -1.0e4])],
            "hitandrun_disrupt"         : [np.array([0.0, 0.0, 0.0]),
                                            np.array([0.0, 0.0, 1.0e5])],
            "hitandrun_pure"            : [np.array([0.0, 0.0, 0.0]),
                                            np.array([0.0, 0.0, 1.0e5])],
            "merge"                     : [np.array([0.0, 0.0, 0.0]),
                                            np.array([0.0, 0.0, 0.0])],
            "merge_spinner"             : [np.array([0.0, 0.0, -1.2e6]),
                                            np.array([0.0, 0.0, 0.0])],
            }

body_Gmass = {"disruption_headon"        : [1e-7, 1e-9],
            "disruption_off_axis"       : [1e-7, 1e-9],
            "supercatastrophic_headon"  : [1e-7, 1e-8],
            "supercatastrophic_off_axis": [1e-7, 1e-8],
            "hitandrun_disrupt"         : [1e-7, 7e-10],
            "hitandrun_pure"            : [1e-7, 7e-10],
            "merge"                     : [1e-7, 1e-8],
            "merge_spinner"             : [1e-7, 1e-8]
            }

tstop = {"disruption_headon"         : 2.0e-3,
        "disruption_off_axis"       : 2.0e-3,
        "supercatastrophic_headon"  : 2.0e-3,
        "supercatastrophic_off_axis": 2.0e-3,
        "hitandrun_disrupt"         : 2.0e-4,
        "hitandrun_pure"            : 2.0e-4,
        "merge"                     : 5.0e-3,
        "merge_spinner"             : 5.0e-3,
        }

nfrag_minimum_expected = {"disruption_headon"         : 20,
                        "disruption_off_axis"       : 20,
                        "supercatastrophic_headon"  : 20,
                        "supercatastrophic_off_axis"  : 20,
                        "hitandrun_disrupt"         : 20,
                        "hitandrun_pure"            : 0,
                        "merge"                     : 0,
                        "merge_spinner"             : 0,
                        }

nfrag_maximum_expected = {"disruption_headon"         : 100,
                        "disruption_off_axis"       : 100,
                        "supercatastrophic_headon"  : 100,
                        "supercatastrophic_off_axis"  : 100,
                        "hitandrun_disrupt"         : 100,
                        "hitandrun_pure"            : 0,
                        "merge"                     : 0,
                        "merge_spinner"             : 0,
                        }

expected_regime = {"disruption_headon"         : "Disruption",
                "disruption_off_axis"       : "Disruption",
                "supercatastrophic_headon"  : "Supercatastrophic disruption",
                "supercatastrophic_off_axis" : "Supercatastrophic disruption",
                "hitandrun_disrupt"         : "Hit and run",
                "hitandrun_pure"            : "Hit and run",
                "merge"                   : "Merge",
                "merge_spinner"           : "Hit and run"
                }
expected_outcome =    {"disruption_headon"         : "calculation converged",
                "disruption_off_axis"       : "calculation converged",
                "supercatastrophic_headon"  : "calculation converged",
                "supercatastrophic_off_axis" : "calculation converged",
                "hitandrun_disrupt"         : "calculation converged",
                "hitandrun_pure"            : "No new fragments generated",
                "merge"                   : "Merging",
                "merge_spinner"           : "No new fragments generated"
                }

density = 3000 * swiftest.AU2M**3 / swiftest.MSun
GU = swiftest.GMSun * swiftest.YR2S**2 / swiftest.AU2M**3
body_radius = body_Gmass.copy()
for k,v in body_Gmass.items():
    body_radius[k] = [((Gmass/GU)/(4./3.*np.pi*density))**(1./3.) for Gmass in v]

body_radius["hitandrun_disrupt"] = [7e-6, 3.25e-6]
body_radius["hitandrun_pure"] = [7e-6, 3.25e-6]


class TestSwiftestRestart(unittest.TestCase):
    def setUp(self):
            # Initialize a target and surface for testing
            self.tmpdir = tempfile.TemporaryDirectory()
            self.simdir = self.tmpdir.name

            self.tmpdir2 = tempfile.TemporaryDirectory()
            self.simdir_repeat = self.tmpdir2.name

    def tearDown(self):
            # Clean up temporary directory
            self.tmpdir.cleanup()

    def test_restart_collision_snapshot(self):
            """
            Test that a simulation with collisions restarts
            """

            sim = swiftest.Simulation(simdir = self.simdir, integrator = 'symba',
                                dump_cadence = 1, tstep_out = 1 * 365.25,
                                tstop = 365.25 * 1, dt = 0.01,
                                DU = 'km', TU = 'd', MU2KG = 1e15)
            sim.add_solar_system_body(['Mars', 'Phobos', 'Deimos'])

            # Add Mars c_lm from GMM3 (hardcoded so we don't have to rely on pyshtools being installed)
            c_lm = np.array(
                    [[[ 1.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00],
                    [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00],
                    [-8.75021132e-04,  5.90314960e-10, -8.46359039e-05, 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00],
                    [-1.18960349e-05,  3.80429811e-06, -1.59480091e-05, 3.50541857e-05,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00],
                    [ 5.12940730e-06,  4.21623326e-06, -9.52539393e-07, 6.45688952e-06,  3.09718541e-07,  0.00000000e+00, 0.00000000e+00],
                    [-1.72668240e-06,  4.83798784e-07, -4.29814806e-06, 3.31265594e-06, -4.64037477e-06, -4.44982150e-06, 0.00000000e+00],
                    [ 1.34648814e-06,  1.80236722e-06,  8.61611560e-07, 9.55649930e-07,  1.00849633e-06,  1.65767871e-06, 2.76208184e-06]],
                    [[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00],
                    [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00],
                    [ 0.00000000e+00, -4.94336174e-11,  4.89346259e-05, 0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00],
                    [ 0.00000000e+00,  2.51773217e-05,  8.36151763e-06, 2.55745706e-05,  0.00000000e+00,  0.00000000e+00, 0.00000000e+00],
                    [ 0.00000000e+00,  3.76338480e-06, -8.98081277e-06, -1.93196208e-07, -1.28730468e-05,  0.00000000e+00, 0.00000000e+00],
                    [ 0.00000000e+00,  2.12322092e-06, -1.16588966e-06, 2.71791890e-07, -3.38212383e-06,  3.77983077e-06, 0.00000000e+00],
                    [ 0.00000000e+00, -1.51845907e-06,  1.46912740e-06, 3.33031890e-07,  2.63880041e-06,  1.62290102e-06, 8.21838201e-07]]]
                    )

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

            sim_restart = swiftest.Simulation(simdir = self.simdir, read_data=True, param_file = 'param.restart.in', tstop=365.25 * 10)


            try:
                sim_restart.run()
            except Exception as e:
                self.fail(f'Failed restart with Exception: {e}')

            self.assertTrue(os.path.exists(os.path.join(self.simdir,"param.00000000000000365250.in")))
            self.assertEqual(sim_restart.data.time.size, 11)
            return

    def test_restart_accurate(self):
            '''
            Test that a restarted run gives the same outputs as a full run

            '''
            sim = swiftest.Simulation(simdir = self.simdir, integrator = 'symba',
                                        tstep_out = 1, dump_cadence = 5,
                                        tstop = 10, dt = 0.01,
                                        DU = 'AU', TU = 'y', MU2KG = 1e15,
                                        compute_conservation_values = True)
            sim.add_solar_system_body(['Sun', 'Mercury', 'Venus', 'Earth', 'Mars'])

            try:
                    sim.run()
            except Exception as e:
                    self.fail(f'Failed initial run with Exception: {e}')

            # repeat run with exact same parameters

            sim_repeat = swiftest.Simulation(simdir = self.simdir_repeat, integrator = 'symba',
                                                tstep_out = 1, dump_cadence = 5,
                                                tstop = 10, dt = 0.01,
                                                DU = 'AU', TU = 'y', MU2KG = 1e15,
                                                compute_conservation_values = True)
            sim_repeat.add_solar_system_body(['Sun', 'Mercury', 'Venus', 'Earth', 'Mars'])

            try:
                    sim_repeat.run()
            except Exception as e:
                    self.fail(f'Failed repeat run with Exception: {e}')

                    # Check that the output files are the same

            for var in sim.data.data_vars:
                    if (sim.data[var].dtype.type != np.str_ and np.isnan(sim.data[var].values).any()):
                        idx = np.where(~np.isnan(sim.data[var].values))
                        self.assertTrue((sim.data[var].values[idx] == sim_repeat.data[var].values[idx]).all(), f'{var} values are not equal\n\n{sim.data[var].values[idx]}\n\nREPEAT\n{sim_repeat.data[var].values[idx]}')
                    else:
                        self.assertTrue((sim.data[var].values == sim_repeat.data[var].values).all(), f'{var} values are not equal\n\n{sim.data[var].values}\n\nREPEAT\n{sim_repeat.data[var].values}')


            # restarted run (from halfway mark in this case)

            sim_restart = swiftest.Simulation(simdir = self.simdir,
                                                read_data = True,
                                                param_file = 'param.00000000000000000500.in',
                                                compute_conservation_values = True)

            try:
                    sim_restart.run()
            except Exception as e:
                    self.fail(f'Failed restart run with Exception: {e}')

                    # check that the output files are the same

            for var in sim.data.data_vars:
                    if (sim.data[var].dtype.type != np.str_ and np.isnan(sim.data[var].values).any()):
                        idx = np.where(~np.isnan(sim.data[var].values))
                        self.assertTrue((sim.data[var].values[idx] == sim_restart.data[var].values[idx]).all(), f'{var} values are not equal\n\n{sim.data[var].values[idx]}\n\nRESTART\n{sim_restart.data[var].values[idx]}')
                    else:
                        self.assertTrue((sim.data[var].values == sim_restart.data[var].values).all(), f'{var} values are not equal\n\n{sim.data[var].values}\n\nRESTART\n{sim_restart.data[var].values}')


    def test_restart_accurate_collision(self):
            '''
            Test that a restarted run with collisions gives the same outputs as a full run.
            collision outcomes taken from test_fraggle.py

            '''
            seed = [974899978, -2041855733, 14415898, 615945619, 1818808148, 462758063, 762516551, 1276432020]
            
            for style in collision_type:
                sim = swiftest.Simulation(simdir=self.simdir, rotation=True, compute_conservation_values=True, seed = seed)
                sim.add_solar_system_body("Sun")
                sim.add_body(name=names, Gmass=body_Gmass[style], radius=body_radius[style], rh=pos_vectors[style], vh=vel_vectors[style], rot=rot_vectors[style])

                # Set fragmentation parameters
                minimum_fragment_gmass = 0.01 * body_Gmass[style][1] 
                gmtiny = 0.50 * body_Gmass[style][1] 
                sim.set_parameter(encounter_save="both", 
                                gmtiny=gmtiny, 
                                minimum_fragment_gmass=minimum_fragment_gmass, 
                                nfrag_reduction=nfrag_reduction[style])
                try:
                    sim.run(dt=tstop[style]/4, tstop=tstop[style], istep_out=1, dump_cadence=2)
                except Exception as e:
                    self.fail(f'Failed initial run with Exception: {e}')
                
                # repeat run with exact same parameters
                sim_repeat = swiftest.Simulation(simdir=self.simdir_repeat, rotation=True, compute_conservation_values=True, seed = seed)
                sim_repeat.add_solar_system_body("Sun")
                sim_repeat.add_body(name=names, Gmass=body_Gmass[style], radius=body_radius[style], rh=pos_vectors[style], vh=vel_vectors[style], rot=rot_vectors[style])

                # Set fragmentation parameters
                sim_repeat.set_parameter(encounter_save="both",
                                        gmtiny=gmtiny,
                                        minimum_fragment_gmass=minimum_fragment_gmass,
                                        nfrag_reduction=nfrag_reduction[style])
                try:
                    sim_repeat.run(dt=tstop[style]/4, tstop=tstop[style], istep_out=1, dump_cadence=2)
                except Exception as e:
                    self.fail(f'Failed repeat run with Exception: {e}')

                    # check that the output files are the same

                for var in sim.data.data_vars:
                    if (sim.data[var].dtype.type != np.str_ and np.isnan(sim.data[var].values).any()):
                        idx = np.where(~np.isnan(sim.data[var].values))
                        self.assertTrue((sim.data[var].values[idx] == sim_repeat.data[var].values[idx]).all(), f'{var} values are not equal\n\n{sim.data[var].values[idx]}\n\nREPEAT\n{sim_repeat.data[var].values[idx]}')
                    else:
                        self.assertTrue((sim.data[var].values == sim_repeat.data[var].values).all(), f'{var} values are not equal\n\n{sim.data[var].values}\n\nREPEAT\n{sim_repeat.data[var].values}')

                # restarted run (from the halfway mark in this case)

                restart_time = f'{2 * int(tstop[style]/4)}'
                param_restart = f'param.' + restart_time.zfill(20) + '.in'

                sim_restart = swiftest.Simulation(simdir=self.simdir, read_data=True, param_file=param_restart, compute_conservation_values=True, seed = seed)
                try:
                    sim_restart.run()
                except Exception as e:
                    self.fail(f'Failed restart run with Exception: {e}')

                    # check that the output files are the same
                
                for var in sim.data.data_vars:
                    if (sim.data[var].dtype.type != np.str_ and np.isnan(sim.data[var].values).any()):
                        idx = np.where(~np.isnan(sim.data[var].values))
                        self.assertTrue((sim.data[var].values[idx] == sim_restart.data[var].values[idx]).all(), f'{var} values are not equal\n\n{sim.data[var].values[idx]}\n\nRESTART\n{sim_restart.data[var].values[idx]}')
                    else:
                        self.assertTrue((sim.data[var].values == sim_restart.data[var].values).all(), f'{var} values are not equal\n\n{sim.data[var].values}\n\nRESTART\n{sim_restart.data[var].values}')
                
                # clean 
                sim.clean()
                sim_repeat.clean()
                sim_restart.clean()


            return

if __name__ == '__main__':
    os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
    unittest.main()