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
Reads and processes a Swiftest output file and a Swifter output file.

Inputs
_______
param.swifter.in  : An ASCII file containing the parameters for the Swifter simulation.
param.swiftest.in : An ASCII file containing the parameters for the Swiftest simulation.
bin.dat           : A binary file containing the Swifter simulation output.
data.nc           : A NetCDF file containing the Swiftest simulation output.

Returns
-------
swifter_swiftest_comp.eps : Encapsulated PostScript containing the number of massive bodies and test particles for a Swifter and Swiftest run over time.
"""

import swiftest
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

# Pull in the Swifter Data
swifter_sim = swiftest.Simulation(param_file="param.swifter.in", codename="Swifter").ds

# Pull in the Swiftest Data
swiftest_sim = swiftest.Simulation(param_file="param.swiftest.in").ds

# Number of Bodies in Swifter System
npl_swifter = []
ntp_swifter = []
for i in range(len(swifter_sim['time'])):
    npl_swifter.append(len(swifter_sim['Gmass'][i].values[np.logical_not(np.isnan(swifter_sim['Gmass'][i].values))]))
    tp_only = swifter_sim.where(swifter_sim['id'] >= 159)
    ntp_swifter.append(len(tp_only['a'][i].values[np.logical_not(np.isnan(tp_only['a'][i].values))]))

# Calculate Slope of Decay
npl_swifter_slope = np.polyfit(swifter_sim['time'], npl_swifter, 1)[0]
ntp_swifter_slope = np.polyfit(swifter_sim['time'], ntp_swifter, 1)[0]
npl_swiftest_slope = np.polyfit(swiftest_sim['time'], swiftest_sim['npl'], 1)[0]
ntp_swiftest_slope = np.polyfit(swiftest_sim['time'], swiftest_sim['ntp'], 1)[0]

# Plot
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,4))
ax.plot(swifter_sim['time'], npl_swifter, color='crimson', linestyle='solid', linewidth=1, label=f'Swifter Massive Bodies\nslope={npl_swifter_slope}')
ax.plot(swiftest_sim['time'], swiftest_sim['npl'], color='navy', linestyle='solid', linewidth=1, label=f'Swiftest Massive Bodies\nslope={npl_swiftest_slope}')
ax.plot(swifter_sim['time'], ntp_swifter, color='crimson', linestyle='dotted', linewidth=1, label=f'Swifter Test Particles\nslope={ntp_swifter_slope}')
ax.plot(swiftest_sim['time'], swiftest_sim['ntp'], color='navy', linestyle='dotted', linewidth=1, label=f'Swiftest Test Particles\nslope={ntp_swiftest_slope}')
ax.set_title("Total Number of Bodies in the System")
ax.set_xlabel("Time (y)", fontsize=12)
ax.set_ylabel("n", fontsize=12)
ax.set_xlim(0,1e6)
ax.set_xticks([0, 2e5, 4e5, 6e5, 8e5, 1e6])
ax.set_xticklabels(['0','$2 * 10^5$','$4 * 10^5$','$6 * 10^5$','$8 * 10^5$', '$1 * 10^6$'], size=12)
ax.set_ylim(40,110)
ax.set_yticks([40, 50, 60, 70, 80, 90, 100, 110])
ax.set_yticklabels(['40', '50', '60', '70', '80', '90', '100', '110'], size=12)
ax.tick_params(direction="in", which='both')
plt.legend(fontsize=12, frameon=False, bbox_to_anchor=(1,1), loc="upper left")
plt.tight_layout()
plt.savefig(f'swifter_swiftest_comp.eps', dpi=300,facecolor='white', transparent=False)
plt.show()