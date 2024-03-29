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
import numpy as np
from numpy.random import default_rng
from astroquery.jplhorizons import Horizons
import datetime
import tempfile

rng = default_rng(seed=123)

major_bodies = ["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"]
param = {}

class TestSwiftestIntegration(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir=tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name
        
    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup() 

    def test_integrators(self):
        """
        Tests that Swiftest is able to integrate a collection of massive bodies and test particles with all available integrators
        """ 
        print("\ntest_integrators: Tests that Swiftest is able to integrate a collection of massive bodies and test particles with all available integrators")
        sim = swiftest.Simulation(simdir=self.simdir)

        # Add the modern planets and the Sun using the JPL Horizons Database.
        sim.add_solar_system_body(major_bodies)
        
        # Add 10 user-defined test particles.
        ntp = 10

        name_tp  = [f"TestParticle_{i:02}" for i in range(1,ntp+1)]
        a_tp     = rng.uniform(0.3, 1.5, ntp)
        e_tp     = rng.uniform(0.0, 0.2, ntp)
        inc_tp   = rng.uniform(0.0, 10, ntp)
        capom_tp = rng.uniform(0.0, 360.0, ntp)
        omega_tp = rng.uniform(0.0, 360.0, ntp)
        capm_tp  = rng.uniform(0.0, 360.0, ntp)

        integrators= ["whm","helio","rmvs","symba"]
        sim.add_body(name=name_tp, a=a_tp, e=e_tp, inc=inc_tp, capom=capom_tp, omega=omega_tp, capm=capm_tp)
        sim.set_parameter(tstart=0.0, tstop=0.02, dt=0.01, istep_out=1, dump_cadence=0)
        integrators= ["whm","helio","rmvs","symba"]
        for i in integrators:
            try:
                sim.run(integrator=i)
            except:
                self.fail(f"Failed with integrator {i}")
        return    
            
    def test_conservation(self):
        """
        Tests that Swiftest conserves mass, energy, and momentum to within acceptable tolerances.
        """
        print("\ntest_conservation: Tests that Swiftest conserves mass, energy, and momentum to within acceptable tolerances.")
        
        # Error limits
        L_slope_limit = 1e-9
        E_slope_limit = 1e-7
        GM_limit = 1e-14
 
        sim = swiftest.Simulation(simdir=self.simdir)
        
        sim.add_solar_system_body(major_bodies)
        
        dt = 0.01
        nout = 100
        tstop = 1e3
        tstep_out = tstop / nout
              
        sim.run(tstart=0.0, tstop=tstop, dt=dt, tstep_out=tstep_out, dump_cadence=0, compute_conservation_values=True, integrator="symba")
        
        def fit_func(x,slope,b):
            """
            Linear function for curve fitting
            """
            return slope * x + b
        
        # Calculate the angular momentum error
        sim.data['L_tot'] = sim.data['L_orbit'] + sim.data['L_spin'] + sim.data['L_escape']
        sim.data['DL'] = sim.data['L_tot'] - sim.data['L_tot'].isel(time=0)
        L_error = sim.data['DL'].magnitude() / sim.data['L_tot'].isel(time=0).magnitude()

        # Calculate the energy error
        E_error = (sim.data['TE'] - sim.data['TE'].isel(time=0)) / sim.data['TE'].isel(time=0)

        # Calculate the mass error
        sim.data['GMtot'] = sim.data['Gmass'].sum(dim='name',skipna=True) + sim.data['GMescape']
        GM_error = (sim.data['GMtot'] - sim.data['GMtot'].isel(time=0)) / sim.data['GMtot'].isel(time=0)
        GM_final = GM_error.isel(time=-1).values
        
        # Compute the slope of the error vs time fit
        E_fit_result = E_error.curvefit("time",fit_func)
        L_fit_result = L_error.curvefit("time",fit_func)
        
        E_slope = E_fit_result['curvefit_coefficients'].sel(param='slope').values
        L_slope = L_fit_result['curvefit_coefficients'].sel(param='slope').values
        
        # Print the final errors
        print("\n")
        print(f"Angular momentum error slope: {L_slope:.2e}/{sim.TU_name}")
        print(f"Energy error slope: {E_slope:.2e}/{sim.TU_name}")
        print(f"Final Mass Error: {GM_final:.2e}")

        self.assertLess(np.abs(L_slope),L_slope_limit, msg=f"Angular Momentum Error of {L_slope:.2e}/{sim.TU_name} higher than threshold value of {L_slope_limit:.2e}/{sim.TU_name}")
        self.assertLess(np.abs(E_slope),E_slope_limit, msg=f"Energy Error of {E_slope:.2e}/{sim.TU_name} higher than threshold value of {E_slope_limit:.2e}/{sim.TU_name}")
        self.assertLess(np.abs(GM_final),GM_limit, msg=f"Mass Error of {GM_final:.2e} higher than threshold value of {GM_limit:.2e}")
        return 
        
    def test_gr(self):
        """
        Tests that GR is working correctly by computing the precession of Mercury's longitude of periapsis and comparing it to
        values obtained from the JPL/Horizons ephemeris service
        """
        print("\ntest_gr: Tests that GR is working correctly.")        
       
        # Error limits 
        dvarpi_limit = 1e-3
        integrators= ["whm","helio","rmvs","symba"] 
        
        # Initialize the simulation object as a variable. Define the directory in which the output will be placed.
        tstep_out = 10.0
        sim = swiftest.Simulation(simdir=self.simdir, tstop=1000.0, dt=0.005, tstep_out=tstep_out, dump_cadence=0,general_relativity=True)
        sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])

        # Get the start and end date of the simulation so we can compare with the real solar system.
        start_date = sim.ephemeris_date
        tstop_d = sim.param['TSTOP'] * sim.param['TU2S'] / swiftest.JD2S

        stop_date = (datetime.datetime.fromisoformat(start_date) + datetime.timedelta(days=tstop_d)).isoformat()

        #Get the ephemerides of Mercury for the same timeframe as the simulation.
        obj = Horizons(id='1', location='@sun',
                    epochs={'start':start_date, 'stop':stop_date,
                            'step':'10y'})
        el = obj.elements()
        t = (el['datetime_jd']-el['datetime_jd'][0]) / 365.25
        varpi_obs = el['w'] + el['Omega']
        dvarpi_obs = np.diff(varpi_obs) / np.diff(t)
        dvarpi_obs_mean = np.mean(dvarpi_obs) 
        
        for i in integrators:
            sim.run(integrator=i)
            varpi_sim = sim.data['varpi'].sel(name="Mercury")
            dvarpi_gr = np.diff(varpi_sim) / tstep_out
            dvarpi_err = np.mean(dvarpi_obs - dvarpi_gr) / dvarpi_obs_mean
            print(f'{i}: Mercury precession rate error {dvarpi_err:.2e} "/{sim.TU_name}')
            self.assertLess(np.abs(dvarpi_err),dvarpi_limit,msg=f'{dvarpi_err:.2e} /{sim.TU_name} is higher than threshold value of {dvarpi_limit:.2e} "/{sim.TU_name}')

        return
    
         
if __name__ == '__main__':
    os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
    unittest.main()