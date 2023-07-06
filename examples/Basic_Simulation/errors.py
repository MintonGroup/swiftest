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
Reads and processes a Swiftest output file.

Input
------
data.nc     : A NetCDF file containing the simulation output.

Output
------
None
"""

import swiftest
import xarray as xr
import matplotlib.pyplot as plt

# Read in the simulation output and store it as an Xarray dataset.
ds = swiftest.Simulation(read_data=True)

# Calculate the angular momentum error
ds.data['L_tot'] = ds.data['L_orbit'] + ds.data['L_spin'] + ds.data['L_escape']
ds.data['DL'] = ds.data['L_tot'] - ds.data['L_tot'].isel(time=0)
ds.data['L_error'] = swiftest.tool.magnitude(ds.data,'DL') / swiftest.tool.magnitude(ds.data.isel(time=0), 'L_tot')
L_final = ds.data['L_error'][-1].values

# Calculate the energy error
ds.data['E_error'] = (ds.data['TE'] - ds.data['TE'].isel(time=0)) / ds.data['TE'].isel(time=0)
E_final = ds.data['E_error'][-1].values

# Calculate the mass error
ds.data['GMtot'] = ds.data['Gmass'].sum(dim='name',skipna=True) + ds.data['GMescape']
ds.data['GM_error'] = (ds.data['GMtot'] - ds.data['GMtot'].isel(time=0)) / ds.data['GMtot'].isel(time=0)
GM_final = ds.data['GM_error'][-1].values

# Print the final errors
print("Final Angular Momentum Error: ", L_final)
print("Final Energy Error: ", E_final)
print("Final Mass Error: ", GM_final)

# Determine if the errors are within bounds
L_limit = 1e-10
E_limit = 1e-5
GM_limit = 1e-14

lerror = 0
if (L_final > L_limit):
   lerror = 1
   print("Angular Momentum Error of", L_final, " higher than threshold value of", L_limit, ". Test failed.")
if (E_final > E_limit):
   lerror = 1
   print("Energy Error of", E_final, " higher than threshold value of", E_limit, ". Test failed.")
if (GM_final > GM_limit):
   lerror = 1
   print("Mass Error of", GM_final, " higher than threshold value of", GM_limit, ". Test failed.")
if (lerror == 0):
   print("Test successful.")
