"""
 Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest. 
 If not, see: https://www.gnu.org/licenses. 
"""

#!/usr/bin/env python3
"""
Reads and processes a Swiftest output file.

Input
-------
out.nc     : NetCDF file
    Swiftest output file.

Returns
-------
output.eps : Encapsulated PostScript file.
    A figure containing the eccentricity vs. semi-major axis for all bodies at the start of the simulation.
"""

import swiftest
import xarray as xr
import matplotlib.pyplot as plt

# Read in the simulation output and store it as an Xarray dataset
ds = swiftest.Simulation(param_file="simdata/param.in", read_old_output_file=True).data

# Plot of the data and save the output plot
colors = ['white' if x == 'Massive Body' else 'black' for x in ds['particle_type']]
sizes = [100 if x == 'Massive Body' else 10 for x in ds['particle_type']]

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,3))
ax.set(xlabel="Semimajor Axis (AU)", ylabel="Eccentricity", title="Simulation Start")
ax.scatter(ds['a'].isel(time=0), ds['e'].isel(time=0), c=colors, s=sizes, edgecolor='black')
ax.set_xlim(0, 2.0)
ax.set_ylim(0, 0.4)
ax.text(1.5, 0.35, f"t = {ds['time'].isel(time=0).values} years", size=10, ha="left")
plt.tight_layout()
plt.show()
fig.savefig("output.eps", dpi=300, facecolor='white', transparent=False, bbox_inches="tight")