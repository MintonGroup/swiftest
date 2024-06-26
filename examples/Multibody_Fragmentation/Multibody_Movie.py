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

"""
Generates, runs, and processes a set of initial conditions for a multi-body super-catastrophic distruption collisional event. 
All Swiftest output files are stored in the /supercatastrophic_multi subdirectory.

**NOTE: You must have ffmpeg installed on your system before running this script. For instance, on MacOS:

```brew install ffmpeg```

on Ubuntu:

```sudo apt-get install ffmpeg```


Inputs
_______
None.

Returns
-------
supercatastrophic_multi.mp4 : A .mp4 file of a fragmentation event.
collisions.log   : An ASCII file containing the information of any collisional events that occured.
collisions.nc    : A NetCDF file containing the collision output.
data.nc          : A NetCDF file containing the simulation output.
encounters.nc    : A NetCDF file containing the encounters output.
init_cond.nc     : A NetCDF file containing the initial conditions for the simulation.
param.00....in   : A series of parameter input files containing the parameters for the simulation at every output stage.
param.in         : An ASCII file containing the inital parameters for the simulation.
param.restart.in : An ASCII file containing the parameters for the simulation at the last output. 
swiftest.log     : An ASCII file containing the information on the status of the simulation as it runs.
"""

import swiftest
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.spatial.transform import Rotation as R

# ----------------------------------------------------------------------------------------------------------------------
# Define the names and initial conditions of the various fragmentation simulation types
# ----------------------------------------------------------------------------------------------------------------------
available_movie_styles = ["supercatastrophic_multi"]
movie_title_list = ["Multi-body Supercatastrophic"]
movie_titles = dict(zip(available_movie_styles, movie_title_list))
num_movie_frames = 1000


# ----------------------------------------------------------------------------------------------------------------------
# To increase the number of bodies generated in each collision type, decrease the value of the corresponding nfrag_reduction number 
# ----------------------------------------------------------------------------------------------------------------------
nfrag_reduction = {
    "supercatastrophic_multi" : 5.0,
    }



# These initial conditions were generated by trial and error
names = ["Body1","Body2","Body3","Body4"]

rdistance = 1e-4/np.sqrt(2.0)
vrel = 2.0/np.sqrt(2.0)

rcb = np.array([1.0, 0.0,  0.0])
vcb = np.array([0.0, 6.28, 0.0])

r1 = rcb + np.array([-rdistance,-rdistance,0.0])
r2 = rcb + np.array([ rdistance,-rdistance,0.0])
r3 = rcb + np.array([ rdistance, rdistance,0.0])
r4 = rcb + np.array([-rdistance, rdistance,0.0])

v1 = vcb + np.array([ vrel, vrel, 0.0])
v2 = vcb + np.array([-vrel, vrel, 0.0])
v3 = vcb + np.array([-vrel,-vrel, 0.0])
v4 = vcb + np.array([ vrel,-vrel, 0.0])

rot1 = np.array([0.0,0.0,1e5])
rot2 = np.array([0.0,0.0,-2e5])
rot3 = np.array([0.0,0.0,3e5])
rot4 = np.array([0.0,0.0,0.0])

m1 = 2e-7
m2 = 2e-7
m3 = 2e-7
m4 = 2e-7

pos_vectors = {"supercatastrophic_multi"   : [r1, r2, r3, r4],
               }

vel_vectors = {"supercatastrophic_multi":   [v1, v2, v3, v4]
               }

rot_vectors = { "supercatastrophic_multi":   [rot1, rot2, rot3, rot4]
               }

body_Gmass = {"supercatastrophic_multi":   [m1, m2, m3, m4]
               }

tstop = {"supercatastrophic_multi"   : 2.0e-3,
         }

density = 3000 * swiftest.AU2M**3 / swiftest.MSun
GU = swiftest.GMSun * swiftest.YR2S**2 / swiftest.AU2M**3
body_radius = body_Gmass.copy()
for k,v in body_Gmass.items():
    body_radius[k] = [((Gmass/GU)/(4./3.*np.pi*density))**(1./3.) for Gmass in v]

# ----------------------------------------------------------------------------------------------------------------------
# Define the animation class that will generate the movies of the fragmentation outcomes
# ----------------------------------------------------------------------------------------------------------------------
def encounter_combiner(sim):
    """
    Combines simulation data with encounter data to produce a dataset that contains the position,
    mass, radius, etc. of both. It will interpolate over empty time values to fill in gaps.
    """

    # Only keep a minimal subset of necessary data from the simulation and encounter datasets
    keep_vars = ['name','rh','vh','Gmass','radius','rot']
    data = sim.data[keep_vars]
    enc = sim.encounters[keep_vars].load()

    # Remove any encounter data at the same time steps that appear in the data to prevent duplicates
    t_not_duplicate = ~enc['time'].isin(data['time'])
    enc = enc.where(t_not_duplicate,drop=True)
    tgood=enc.time.where(~np.isnan(enc.time),drop=True)
    enc = enc.sel(time=tgood)

    # The following will combine the two datasets along the time dimension, sort the time dimension, and then fill in any time gaps with interpolation
    ds = xr.combine_nested([data,enc],concat_dim='time').sortby("time").interpolate_na(dim="time")
    
    # Interpolate in time to make a smooth, constant time step dataset 
    # Add a bit of padding to the time, otherwise there are some issues with the interpolation in the last few frames.
    smooth_time = np.linspace(start=tgood.isel(time=0), stop=ds.time[-1], num=int(1.2*num_movie_frames))
    ds = ds.interp(time=smooth_time)
    ds['rotangle'] = xr.zeros_like(ds['rot'])
    ds['rot'] = ds['rot'].fillna(0.0)

    return ds

def center(Gmass, x, y):
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    Gmass = Gmass[~np.isnan(Gmass)]
    x_com = np.sum(Gmass * x) / np.sum(Gmass)
    y_com = np.sum(Gmass * y) / np.sum(Gmass)
    return x_com, y_com 
class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""

    def __init__(self, sim, animfile, title, style, nskip=1):

        self.ds = encounter_combiner(sim)
        self.npl = len(self.ds['name']) - 1
        self.title = title
        self.body_color_list = {'Initial conditions': 'xkcd:windows blue',
                      'Disruption': 'xkcd:baby poop',
                      'Supercatastrophic': 'xkcd:shocking pink',
                      'Hit and run fragmention': 'xkcd:blue with a hint of purple',
                      'Central body': 'xkcd:almost black'}

        # Set up the figure and axes...
        self.figsize = (4,4)
        self.fig, self.ax = self.setup_plot()
        self.ani = animation.FuncAnimation(self.fig, self.update_plot, init_func=self.init_func, interval=1, frames=range(0,num_movie_frames,nskip), blit=True)
        self.ani.save(animfile, fps=60, dpi=300, extra_args=['-vcodec', 'libx264'])
        print(f"Finished writing {animfile}")

    def setup_plot(self):
        fig = plt.figure(figsize=self.figsize, dpi=300)
        plt.tight_layout(pad=0)

        # Calculate the distance along the y-axis between the colliding bodies at the start of the simulation.
        # This will be used to scale the axis limits on the movie.
        rhy1 = self.ds['rh'].sel(name="Body1",space='y').isel(time=0).values[()]
        rhy2 = self.ds['rh'].sel(name="Body3",space='y').isel(time=0).values[()]

        scale_frame =   abs(rhy1) + abs(rhy2)
           
        ax = plt.Axes(fig, [0.1, 0.1, 0.8, 0.8])
        self.ax_pt_size = self.figsize[0] *  72 / scale_frame  * 0.7
        ax.set_xlim(-scale_frame, scale_frame)
        ax.set_ylim(-scale_frame, scale_frame)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(self.title)
        fig.add_axes(ax)

        return fig, ax
    
    def init_func(self):
        self.artists = [] 
        aarg = self.vec_props('xkcd:beige')
        for i in range(self.npl):
            self.artists.append(self.ax.annotate("",xy=(0,0),**aarg)) 
        
        self.artists.append(self.ax.scatter([],[],c='k', animated=True, zorder=10))
        return self.artists

    def update_plot(self, frame):
        # Define a function to calculate a reference frame for the animation
        t, Gmass, rh, radius, rotangle = next(self.data_stream(frame))
        x_ref, y_ref = center(Gmass, rh[:,0], rh[:,1]) 
        rh = np.c_[rh[:,0] - x_ref, rh[:,1] - y_ref]
        self.artists[-1].set_offsets(rh)
        point_rad = radius * self.ax_pt_size
        self.artists[-1].set_sizes(point_rad**2)
        
        sarrowend, sarrowtip = self.spin_arrows(rh, rotangle, 1.1*radius) 
        for i, s in enumerate(self.artists[:-1]):
            self.artists[i].set_position(sarrowtip[i])
            self.artists[i].xy = sarrowend[i]
        
        return self.artists

    def data_stream(self, frame=0):
        while True:
            t = self.ds.isel(time=frame)['time'].values[()]
            
            if frame > 0:
                dsold = self.ds.isel(time=frame-1)
                told = dsold['time'].values[()]
                dt = t - told
                self.ds['rotangle'][dict(time=frame)] = dsold['rotangle'] + dt * dsold['rot']
            
            ds = self.ds.isel(time=frame)
            ds = ds.where(ds['name'] != "Sun", drop=True)
            radius = ds['radius'].values
            Gmass = ds['Gmass'].values
            rh = ds['rh'].values
            rotangle = ds['rotangle'].values

            yield t, Gmass, rh, radius, rotangle
            
    def spin_arrows(self, rh, rotangle, rotlen):
        px = rh[:, 0]
        py = rh[:, 1]
        sarrowend = []
        sarrowtip = []
        for i in range(rh.shape[0]):
            endrel = np.array([0.0, -rotlen[i],  0.0])
            tiprel = np.array([0.0, rotlen[i], 0.0])
            r = R.from_rotvec(rotangle[i,:], degrees=True)
            endrel = r.apply(endrel)
            tiprel = r.apply(tiprel)
            send = (px[i] + endrel[0], py[i] + endrel[1])
            stip = (px[i] + tiprel[0], py[i] + tiprel[1])
            sarrowend.append(send)
            sarrowtip.append(stip)
        return sarrowend, sarrowtip
    
    def vec_props(self, c):
        arrowprops = {
            'arrowstyle': '-',
            'linewidth' : 1,
        }

        arrow_args = {
            'xycoords': 'data',
            'textcoords': 'data',
            'arrowprops': arrowprops,
            'annotation_clip': True,
            'zorder': 100,
            'animated' : True
        }
        aarg = arrow_args.copy()
        aprop = arrowprops.copy()
        aprop['color'] = c
        aarg['arrowprops'] = aprop
        aarg['color'] = c
        return aarg    
    
if __name__ == "__main__":

    print("Select a fragmentation movie to generate.")
    print("1. Multi-body supercatastrophic")
    print("2. All of the above")
    user_selection = int(input("? "))

    if user_selection > 0 and user_selection < 2:
        movie_styles = [available_movie_styles[user_selection-1]]
    else:
        print("Generating all movie styles")
        movie_styles = available_movie_styles.copy()

    for style in movie_styles:
        print(f"Generating {movie_titles[style]}")
        movie_filename = f"{style}.mp4"
        # Pull in the Swiftest output data from the parameter file and store it as a Xarray dataset.
        sim = swiftest.Simulation(simdir=style, rotation=True, init_cond_format = "XV", compute_conservation_values=True)
        sim.add_solar_system_body("Sun")
        sim.add_body(name=names, Gmass=body_Gmass[style], radius=body_radius[style], rh=pos_vectors[style], vh=vel_vectors[style], rot=rot_vectors[style])

        # Set fragmentation parameters
        minimum_fragment_gmass = 0.05 * body_Gmass[style][1] 
        gmtiny = 0.10 * body_Gmass[style][1] 
        sim.set_parameter(collision_model="fraggle", 
                          encounter_save="both", 
                          gmtiny=gmtiny, 
                          minimum_fragment_gmass=minimum_fragment_gmass,
                          nfrag_reduction=nfrag_reduction[style],
                          verbose=False)
        sim.run(dt=5e-4, tstop=tstop[style], istep_out=1, dump_cadence=0)

        print("Generating animation")
        anim = AnimatedScatter(sim,movie_filename,movie_titles[style],style,nskip=1)
