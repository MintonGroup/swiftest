#!/usr/bin/env python3
"""
 Copyright 2023 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest. 
 If not, see: https://www.gnu.org/licenses. 
"""

"""
Tests that energy and momentum errors are within tolerances in a Swiftest simulation

Input
------

Output
------
None
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

class TestSwiftest(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir=tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name
        
    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup() 
        
    def test_xv2el2xv(self):
        """
        Tests that the functions xv2el and el2xv are able to convert between position-velocity and orbital elements without any exceptions being raised
        """
        print("\ntest_xv2el2xv: Tests that the functions xv2el and el2xv are able to convert between position-velocity and orbital elements without any exceptions being raised")
        # Generate a set of random position-velocity vectors
        from swiftest.core import xv2el, el2xv
        
        # Test that we can reliably convert between orbital elements and state vectors
        n = 1000
        mu    = np.ones(n)
        a     = rng.uniform(0.1, 1.5, n)
        e     = rng.uniform(0.0, 2.0, n)
        inc   = rng.uniform(0.0, 180.0, n)
        capom = rng.uniform(0.0, 360.0, n)
        omega = rng.uniform(0.0, 360.0, n)
        capm  = rng.uniform(0.0, 360.0, n)
       
        rh, vh = el2xv(mu, a, e, inc, capom, omega, capm)
        a2, e2, inc2, capom2, omega2, capm2, varpi, lam, f, cape, capf = xv2el(mu, rh, vh) 
        
        # Check that the original and converted position-velocity vectors are the same
        self.assertTrue(np.allclose(a,a2),msg=f"Error converting a: {a}, {a2}")
        self.assertTrue(np.allclose(e,e2),msg=f"Error converting a: {e}, {e2}")
        self.assertTrue(np.allclose(inc,inc2),msg=f"Error converting a: {inc}, {inc2}")
        self.assertTrue(np.allclose(capom,capom2),msg=f"Error converting a: {capom}, {capom2}")
        self.assertTrue(np.allclose(omega,omega2),msg=f"Error converting a: {omega}, {omega2}")
        self.assertTrue(np.allclose(capm,capm2),msg=f"Error converting a: {capm}, {capm2}")
        return
    
    def test_gen_ic(self):
        """
        Tests that Swiftest is able to successfully generate a set of initial conditions in a file without any exceptions being raised
        """
        print("\ntest_gen_ic: Test whether we can generate simulation initial conditions test")
        # Files that are expected to be generated:
        file_list = [self.simdir, os.path.join(self.simdir,"param.in"), os.path.join(self.simdir,"init_cond.nc")]
        
        sim = swiftest.Simulation(simdir=self.simdir)
        sim.clean()

        # Add the modern planets and the Sun using the JPL Horizons Database.
        sim.add_solar_system_body(major_bodies)
        sim.save()
        
        for f in file_list:
            self.assertTrue(os.path.exists(f))

        print("\ntest_read_ic: Test whether we can read back initial conditions files created by test_gen_ic")
        sim2 = swiftest.Simulation(simdir=self.simdir, read_param=True, read_data=False)
        # Add the modern planets and the Sun using the JPL Horizons Database.
        # Check if all names in Dataset read in from file match the expected list of names
        self.assertTrue((major_bodies == sim2.init_cond['name']).all(), msg="Name mismatch in Dataset")
        
        # Check to see if all parameter values read in from file match the expected parameters saved when generating the file
        self.assertTrue(all([v == param[k] for k,v in sim2.param.items() if k in param]))
        return
    
    def test_add_body_combinations(self):
        """
        Tests various combinations of arguments to the Simulation.add_body and Simulation.add_solar_system_body methods to ensure argument checking is working as intended.
        """
        
        sim = swiftest.Simulation(simdir=self.simdir)
        sim.clean() 
        
        # Test that we can can pass a single body:
        sim.add_solar_system_body(name="Sun")
        
        # Test that we can add multiple bodies:
        sim.add_solar_system_body(name=["Mercury","Venus","Earth"])
        
        # Test that we can pass by ephemeris_id alone
        sim.add_solar_system_body(ephemeris_id=["Mars","500"])
      
        # Test that we can initialize a body with only the semimajor axis 
        sim.add_body(a=1.0) 
        self.assertEqual(sim.init_cond.isel(name=-1).e.values[0], 0.0, msg="Failed to initialize body with only semimajor axis")
        self.assertEqual(sim.init_cond.isel(name=-1).inc.values[0], 0.0, msg="Failed to initialize body with only semimajor axis")
        self.assertEqual(sim.init_cond.isel(name=-1).capom.values[0], 0.0, msg="Failed to initialize body with only semimajor axis")
        self.assertEqual(sim.init_cond.isel(name=-1).omega.values[0], 0.0, msg="Failed to initialize body with only semimajor axis")
        self.assertEqual(sim.init_cond.isel(name=-1).capm.values[0], 0.0, msg="Failed to initialize body with only semimajor axis")
        
        # Test that we can input cartesian coordinates
        sim.add_body(mass=1.0, radius=1.0, rh=[1.0,0.0,0.0], vh=[0.0,1.0,0.0])
       
        # orbital elements without semimajor axis 
        with self.assertRaises(ValueError):
           sim.add_body(e=0.0, inc=0.0, capom=0.0, omega=0.0, capm=0.0) 
         
        # Mix of elements and cartesian inputs
        with self.assertRaises(ValueError):
            sim.add_body(a=1.0, e=0.0, rh=[1.0,0.0,0.0], vh=[0.0,1.0,0.0])
            
        # Position but not velocity
        with self.assertRaises(ValueError):
            sim.add_body(rh=[1.0,0.0,0.0])
            
        # Velocity but not position
        with self.assertRaises(ValueError):
            sim.add_body(vh=[0.0,1.0,0.0])
            
        # Add J2 and c_lm values
        with self.assertRaises(ValueError):
            sim.add_body(mass=1.0, radius=1.0, J2=1.0e-6, c_lm=np.ones([2,7]))
            
        # Wrong shape of c_lm
        with self.assertRaises(ValueError):
            sim.add_body(mass=1.0, radius=1.0, c_lm=[1.0,0.0,0.0])
            
        # Mismatched lengths of input arguments
        with self.assertRaises(ValueError):
            sim.add_solar_system_body(name=["Mercury","Venus","Earth"], ephemeris_id=["Mars","500"])
        with self.assertRaises(ValueError):
            sim.add_body(a=[1.0,2.0], e=0.0)
        with self.assertRaises(ValueError):
            sim.add_body(rh=[[1.0,0.0,0.0],[2.0,0.0,0.0]], vh=[0.0,1.0,0.0])     
        with self.assertRaises(ValueError):
            sim.add_body(rh=[1.0,0.0], vh=[0.0,1.0])
      
        # mass and Gmass at the same time 
        with self.assertRaises(ValueError):
            sim.add_body(a=1.0, mass=1.0, radius=1.0, Gmass=4*np.pi**2)
            
        # mass without radius
        with self.assertRaises(ValueError):
            sim.add_body(a=1.0, mass=1.0)
            
        return
    
    def test_mixed_element_cart_input(self):
        """
        Tests that we can mix orbital element and cartesian input
        """
        
        sim = swiftest.Simulation(simdir=self.simdir)
        sim.clean() 
        
        # Test that we can can pass a single body:
        sim.add_solar_system_body(name="Sun")
        
        sim.add_body(a=[1.0])
        sim.add_body(rh=[2.0,0.0,0.0], vh=[0.0,1.0,0.0])
        try:
            sim.run(tstop=0.02, dt=0.01)
        except:
            self.fail("Failed to run simulation with mixed element and cartesian input")
        
        return
    
    def test_read_multi_dir(self):
        """
        Tests that Swiftest can generate a set of initial conditions, copy them into a new directory, and then run those from a
        different Simulation object (test inspired by Kaustub Anand's workflow)
        """
        
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
            os.makedirs(dst, exist_ok=True)

            # Iterate over all items in the source directory
            for item in os.listdir(src):
                src_path = os.path.join(src, item)
                dst_path = os.path.join(dst, item)

                # Copy files and directories
                if os.path.isdir(src_path):
                    # Recursively copy a directory to a destination
                    shutil.copytree(src_path, dst_path, dirs_exist_ok=True)
                else:
                    # Copy a single file
                    shutil.copy(src_path, dst_path)
            
            return        

        simdir1 = os.path.join(self.simdir, "sim1")
        simdir2 = os.path.join(self.simdir, "sim2")
        sim1 = swiftest.Simulation(simdir=simdir1, tstop=1.0,dt=0.01)
        sim1.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars"])
        sim1.save() 
        copy_folder_contents(simdir1, simdir2)
        sim2 = swiftest.Simulation(read_param=True, simdir=simdir2)
        try:
            sim2.run()
        except:
            self.fail("Failed to run simulation from copied directory")
        
        return
        
    def test_integrators(self):
        """
        Tests that Swiftest is able to integrate a collection of massive bodies and test particles with all available integrators
        """ 
        print("\ntest_integrators: Tests that Swiftest is able to integrate a collection of massive bodies and test particles with all available integrators")
        sim = swiftest.Simulation(simdir=self.simdir)

        # Add the modern planets and the Sun using the JPL Horizons Database.
        sim.add_solar_system_body(major_bodies)
        sim.clean()
        
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
        sim.clean()
        
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
    
    def test_planetocentric(self):
        """
        Tests that Swiftest is able to set up a simulation in a planetocentric frame and that the results are consistent with the expected values
        """
        print("\ntest_planetocentric: Compare an unrotated body-centric frame in Swiftest to the JPL Horizons ephemeris ")
        sim = swiftest.Simulation(simdir=self.simdir)
        sim.clean()         
        
        cbname = "Mars"
        cbcenter = "@4"
        satname = ["Phobos","Deimos"]
        sim.add_solar_system_body(cbname, is_central_body=True)
        sim.add_solar_system_body(satname)

        def get_jpl_data(id):
            tstart = datetime.date.fromisoformat(sim.ephemeris_date)
            tstep = datetime.timedelta(days=1)
            tend = tstart + tstep
            ephemerides_start_date = sim.ephemeris_date
            ephemerides_end_date = tend.isoformat()
            ephemerides_step = '1d' 
            
            jpl = Horizons(id=id, location=cbcenter,
                epochs={'start': ephemerides_start_date, 'stop': ephemerides_end_date,
                'step': ephemerides_step})
            eph = jpl.ephemerides()
            DCONV = swiftest.AU2M / sim.param['DU2M']
            VCONV = (swiftest.AU2M / swiftest.JD2S) / (sim.param['DU2M'] / sim.param['TU2S'])    
        
            vec_tab = jpl.vectors() 
            rx = vec_tab['x'][0] * DCONV
            ry = vec_tab['y'][0] * DCONV
            rz = vec_tab['z'][0] * DCONV
            vx = vec_tab['vx'][0] * VCONV
            vy = vec_tab['vy'][0] * VCONV
            vz = vec_tab['vz'][0] * VCONV
            
            elem_tab = jpl.elements()
            
            a = elem_tab['a'][0] * DCONV
            e = elem_tab['e'][0]
            inc = elem_tab['incl'][0]
            capom = elem_tab['Omega'][0]
            omega = elem_tab['w'][0]
            capm = elem_tab['M'][0]

            rh = np.array([rx,ry,rz]) 
            vh = np.array([vx,vy,vz]) 
            elem = np.array([a,e,inc,capom,omega,capm])
            return rh,vh,elem
                
        for sat in satname:
            sim_rh = sim.init_cond.sel(name=sat,time=0).rh.values
            sim_vh = sim.init_cond.sel(name=sat,time=0).vh.values
            sim_elem = np.array([sim.init_cond[var].sel(name=sat,time=0).values.item() for var in ['a','e','inc','capom','omega','capm']])
            
            jpl_rh,jpl_vh,jpl_elem = get_jpl_data(sat)
            self.assertTrue(np.allclose(sim_rh,jpl_rh),msg=f"Error in rh for {sat}: {sim_rh - jpl_rh}")
            self.assertTrue(np.allclose(sim_vh,jpl_vh),msg=f"Error in vh for {sat}: {sim_vh - jpl_vh}") 
            self.assertTrue(np.allclose(sim_elem,jpl_elem, rtol=1e-4),msg=f"Error in elements for {sat}: {sim_elem - jpl_elem}") 
            
        return
   
    def test_central_body_rotation(self):
        """
        Tests that Swiftest is able to rotate solar system bodies into the equatorial frame of a central body correctly.
        """
        print("\ntest_central_body_rotation: Tests that Swiftest is able to rotate solar system bodies into the equatorial frame of a central body correctly")
        sim = swiftest.Simulation(simdir=self.simdir)
        sim.clean()
        # Precomputed values for [Mars, Phobos, and Deimos] in the ecliptic and Mars equatorial frames
        inc_vals = {"ecliptic": np.array([ np.nan, 26.55406151, 24.06273241]),
                    "mars equator": np.array([np.nan, 1.0753048, 2.6925768])}
        rot_vals = {"ecliptic": np.array([[ 57176.72129971,  -7162.31013346, 114478.65911179],
                                                [183617.63838933, -15311.68519899, 368719.81695279], 
                                                [ 42015.96468119,  -6035.89004401,  95062.95557518]]),
                    "mars equator": np.array([[-7.60184592e-28,  0.00000000e+00,  1.28163331e+05],
                                        [-3.81633732e+02,  7.73720290e+03,  4.12121558e+05],
                                        [-4.89032977e+03, -1.60117170e+02,  1.03994220e+05]])}
        
        # No rotation should keep it in the ecliptic frame 
        sim.add_solar_system_body(["Mars","Phobos","Deimos"]) # Default should be ecliptic
        sim_inc = sim.init_cond.isel(time=0).inc.values
        sim_rot = sim.init_cond.isel(time=0).rot.values

        inc_close = np.allclose(sim_inc,inc_vals["ecliptic"],equal_nan=True)
        rot_close = np.allclose(sim_rot,rot_vals["ecliptic"])
        
        self.assertTrue(inc_close,msg=f"Error in inclination 1")
        self.assertTrue(rot_close,msg=f"Error in rotation 1")

        # Rotating to match the Mars pole
        sim = swiftest.Simulation()
        sim.clean()

        sim.add_solar_system_body(["Mars","Phobos","Deimos"],align_to_central_body_rotation=True)
        sim_inc = sim.init_cond.isel(time=0).inc.values
        sim_rot = sim.init_cond.isel(time=0).rot.values

        inc_close = np.allclose(sim_inc,inc_vals["mars equator"],equal_nan=True)
        rot_close = np.allclose(sim_rot,rot_vals["mars equator"])        
       
        self.assertTrue(inc_close,msg=f"Error in inclination 2")
        self.assertTrue(rot_close,msg=f"Error in rotation 2")
        
        sim = swiftest.Simulation()
        sim.clean()

        # Mix and match
        sim.add_solar_system_body("Mars") # ecliptic
        sim.add_solar_system_body(["Phobos","Deimos"],align_to_central_body_rotation=True) # mars equator
        sim_inc = sim.init_cond.isel(time=0).inc.values
        sim_rot = sim.init_cond.isel(time=0).rot.values

        inc_close = np.allclose(sim_inc[1:],inc_vals["mars equator"][1:]) 
        rot_close = np.allclose(sim_rot[0:1],rot_vals["ecliptic"][0:1]) and np.allclose(sim_rot[1:],rot_vals["mars equator"][1:]) 
        
        self.assertTrue(inc_close,msg=f"Error in inclination 3")
        self.assertTrue(rot_close,msg=f"Error in rotation 3")

        sim = swiftest.Simulation()
        sim.clean()

        sim.add_solar_system_body("Mars",align_to_central_body_rotation=True)
        sim.add_solar_system_body(["Phobos","Deimos"]) # ecliptic 
        sim_inc = sim.init_cond.isel(time=0).inc.values
        sim_rot = sim.init_cond.isel(time=0).rot.values

        inc_close = np.allclose(sim_inc[1:],inc_vals["ecliptic"][1:]) 
        rot_close = np.allclose(sim_rot[0:1],rot_vals["mars equator"][0:1]) and np.allclose(sim_rot[1:],rot_vals["ecliptic"][1:]) 
        
        self.assertTrue(inc_close,msg=f"Error in inclination 4")
        self.assertTrue(rot_close,msg=f"Error in rotation 4")

        sim = swiftest.Simulation()
        sim.clean()


        sim.add_solar_system_body("Mars",align_to_central_body_rotation=True)
        sim.add_solar_system_body("Phobos") # ecliptic
        sim.add_solar_system_body("Deimos",align_to_central_body_rotation=True) # mars equator
        sim_inc = sim.init_cond.isel(time=0).inc.values
        sim_rot = sim.init_cond.isel(time=0).rot.values

        inc_close = np.allclose(sim_inc[1:2],inc_vals["ecliptic"][1:2]) and np.allclose(sim_inc[2:],inc_vals["mars equator"][2:])
        rot_close = np.allclose(sim_rot[0:1],rot_vals["mars equator"][0:1]) and np.allclose(sim_rot[1:2],rot_vals["ecliptic"][1:2]) and np.allclose(sim_rot[2:],rot_vals["mars equator"][2:])
        
        self.assertTrue(inc_close,msg=f"Error in inclination 5")
        self.assertTrue(rot_close,msg=f"Error in rotation 5")

        sim = swiftest.Simulation()
        sim.clean()

        sim.add_solar_system_body("Mars",align_to_central_body_rotation=True)
        sim.add_solar_system_body("Phobos",align_to_central_body_rotation=True) # mars equator
        sim.add_solar_system_body("Deimos") # ecliptic
        sim_inc = sim.init_cond.isel(time=0).inc.values
        sim_rot = sim.init_cond.isel(time=0).rot.values

        inc_close = np.allclose(sim_inc[1:2],inc_vals["mars equator"][1:2]) and np.allclose(sim_inc[2:],inc_vals["ecliptic"][2:])
        rot_close = np.allclose(sim_rot[0:2],rot_vals["mars equator"][0:2]) and np.allclose(sim_rot[2:2],rot_vals["ecliptic"][2:]) 
        
        self.assertTrue(inc_close,msg=f"Error in inclination 6")
        self.assertTrue(rot_close,msg=f"Error in rotation 6")
        
        return 

         
if __name__ == '__main__':
    os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
    unittest.main()