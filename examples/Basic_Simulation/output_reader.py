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
ds = swiftest.Simulation(param_file="param.in").ds

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