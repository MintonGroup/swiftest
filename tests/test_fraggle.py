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

density = 3000 * swiftest.AU2M**3 / swiftest.MSun
GU = swiftest.GMSun * swiftest.YR2S**2 / swiftest.AU2M**3
body_radius = body_Gmass.copy()
for k,v in body_Gmass.items():
    body_radius[k] = [((Gmass/GU)/(4./3.*np.pi*density))**(1./3.) for Gmass in v]

body_radius["hitandrun_disrupt"] = [7e-6, 3.25e-6] 
body_radius["hitandrun_pure"] = [7e-6, 3.25e-6] 


class TestFraggle(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir=tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name
        
    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup() 
        
    def test_disruption_headon(self):
        '''
        Check that the head on disruption collision generates fragments and conserves quantities
        '''
        
        style = "disruption_headon""
        sim = swiftest.Simulation(simdir=self.simdir, rotation=True, compute_conservation_values=True)
        sim.add_solar_system_body("Sun")
        sim.add_body(name=names, Gmass=body_Gmass[style], radius=body_radius[style], rh=pos_vectors[style], vh=vel_vectors[style], rot=rot_vectors[style])

        # Set fragmentation parameters
        minimum_fragment_gmass = 0.01 * body_Gmass[style][1] 
        gmtiny = 0.50 * body_Gmass[style][1] 
        sim.set_parameter(collision_model="fraggle", 
                          encounter_save="both", 
                          gmtiny=gmtiny, 
                          minimum_fragment_gmass=minimum_fragment_gmass, 
                          nfrag_reduction=nfrag_reduction[style])
        sim.run(dt=5e-4, tstop=tstop[style], istep_out=1, dump_cadence=0)
        
        collision_logfile = os.path.join(self.simdir, "collisions.log")
        with open(collision_logfile, "r") as f:
            content = f.read()
           
        self.assertIn("calculation converged", content, "The output file does not contain 'calculation converged'")
     

        return 
         
if __name__ == '__main__':
    os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
    unittest.main()