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
Generates a movie of a fragmentation event from set of Swiftest output files.

Inputs
_______
param.in : ASCII text file
    Swiftest parameter input file.
out.nc   : NetCDF file
    Swiftest output file.

Returns
-------
fragmentation.mp4 : mp4 movie file
    Movide of a fragmentation event.
"""

import swiftest
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation 

# Change this to be the parameter input file correlated with the run that you
# wish to test. Swiftest will pull in the corresponding out.nc file automatically.
param_file = "param.disruption_headon.in"

# Change this to an appropriate title and filename to appear on the movie.
movie_title = "Head-on Disruption"
movie_filename = "disruption_headon.mp4"

# Pull in the Swiftest output data from the parameter file and store it as a Xarray dataset.
ds = swiftest.Simulation(read_param=True, param_file=param_file, read_old_output_file=True).data

# Calculate the number of frames in the dataset.
nframes = int(ds['time'].size)

# Define a function to calculate the center of mass of the system.
def center(xhx, xhy, xhz, Gmass):
    x_com = np.sum(Gmass * xhx) / np.sum(Gmass)
    y_com = np.sum(Gmass * xhy) / np.sum(Gmass)
    z_com = np.sum(Gmass * xhz) / np.sum(Gmass)
    return x_com, y_com, z_com

# Calculate the distance along the y-axis between the colliding bodies at the start of the simulation. 
# This will be used to scale the axis limits on the movie.
scale_frame = abs(ds['xhy'].isel(time=0, name=1).values) + abs(ds['xhy'].isel(time=0, name=2).values)

# Set up the figure and the animation.
fig, ax = plt.subplots(figsize=(4,4))
def animate(i):
    # Calculate the position and mass of all bodies in the system at time i and store as a numpy array.
    xhx = ds['xhx'].isel(time=i).dropna(dim='name').values
    xhy = ds['xhy'].isel(time=i).dropna(dim='name').values
    xhz = ds['xhx'].isel(time=i).dropna(dim='name').values
    Gmass = ds['Gmass'].isel(time=i).dropna(dim='name').values[1:] # Drop the Sun from the numpy array.

    # Calculate the center of mass of the system at time i. While the center of mass relative to the 
    # colliding bodies does not change, the center of mass of the collision will move as the bodies
    # orbit the system center of mass.
    x_com, y_com, z_com = center(xhx, xhy, xhz, Gmass)
    
    # Create the figure and plot the bodies as points.
    fig.clear()
    ax = fig.add_subplot(111)
    ax.set_title(movie_title)
    ax.set_xlabel("xhx")
    ax.set_ylabel("xhy")
    ax.set_xlim(x_com - scale_frame, x_com + scale_frame)
    ax.set_ylim(y_com - scale_frame, y_com + scale_frame)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.scatter(xhx, xhy, s = (5000000000 * Gmass))
    
    plt.tight_layout()

# Generate the movie.
ani = animation.FuncAnimation(fig, animate, interval=1, frames=nframes, repeat=False)
ani.save(movie_filename, fps=60, dpi=300, extra_args=['-vcodec', 'libx264'])