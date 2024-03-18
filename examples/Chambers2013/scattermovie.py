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
Creates a movie from a set of Swiftest output files. All simulation 
outputs are stored in the /simdata subdirectory.

**NOTE: You must have ffmpeg installed on your system before running this script. For instance, on MacOS:

```brew install ffmpeg```

on Ubuntu:

```sudo apt-get install ffmpeg```


Input
------
param.in    : ASCII Swiftest parameter input file.
data.nc     : A NetCDF file containing the simulation output.

Output
------
Chambers2013-aescatter.mp4  : A .mp4 file plotting eccentricity vs semimajor axis. 
OR 
Chambers2013-aiscatter.mp4  : A .mp4 file plotting inclination vs semimajor axis. 
"""

import swiftest 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.colors as mcolors
from collections import namedtuple


plt.switch_backend('agg')

titletext = "Chambers (2013)"
valid_plot_styles = ["aescatter", "aiscatter", "arotscatter"]
xlim={"aescatter" : (0.0, 3.0),
      "aiscatter" : (0.0, 3.0), 
      "arotscatter" : (0.0, 2.5)}
ylim={"aescatter" : (0.0, 1.0),
      "aiscatter" : (0.0, 40.0),
      "arotscatter" : (1.0, 10000.0)}
xlabel={"aescatter": "Semimajor axis (AU)",
        "aiscatter": "Semimajor axis (AU)",
        "arotscatter": "Semimajor axis (AU)"}
ylabel={"aescatter": "Eccentricity",
        "aiscatter": "Inclination (deg)",
        "arotscatter": "Rotation period (h)"}

YR2HR = 365.25 * 24
ROT2PERIOD = YR2HR * 360.0

framejump = 1
origin_types = ["Initial conditions", "Merger", "Disruption", "Supercatastrophic", "Hit and run fragmentation"]


class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, ds, param):

        self.radscale = 2000
        nframes = int(ds['time'].size / framejump)
        ds['rot_mag'] = ds['rot'].magnitude()
        ds['rot_mag'] = ROT2PERIOD / ds['rot_mag']
        self.Rcb = ds['radius'].sel(name="Sun").isel(time=0).values[()]
        self.ds = ds
        self.param = param
        colors = ["k", "xkcd:faded blue", "xkcd:marigold", "xkcd:shocking pink", "xkcd:baby poop green"]
        self.clist = dict(zip(origin_types,colors))

        # Setup the figure and axes...
        fig = plt.figure(figsize=(8,4.5), dpi=300)
        plt.tight_layout(pad=0)
        # set up the figure
        self.ax = plt.Axes(fig, [0.1, 0.15, 0.8, 0.75])
        fig.add_axes(self.ax)

        self.make_artists()
        
        self.ani = animation.FuncAnimation(fig, func=self.update, interval=1, frames=nframes, init_func=self.init_plot, blit=True)
        self.ani.save(animation_file, fps=60, dpi=300, extra_args=['-vcodec', 'libx264'])
        print(f'Finished writing {animation_file}')
        
    def make_artists(self):
        scatter_names = [f"s{i}" for i,k in enumerate(origin_types)]
        self.scatter_artist_names = dict(zip(origin_types,scatter_names))
        
        animated_elements = [self.ax.text(0.50, 1.05, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5}, transform=self.ax.transAxes,ha="center", animated=True)]
        element_names = ["title"]
        for key, value in self.clist.items():
            animated_elements.append(self.ax.scatter([], [], marker='o', s=[], c=value, alpha=0.75, label=key, animated=True))
            element_names.append(self.scatter_artist_names[key])
        
        Artists = namedtuple("Artists",tuple(element_names))
        self.artists = Artists(*animated_elements)
        return 

    def init_plot(self):
        self.ax.set_xlim(xlim[plot_style])
        self.ax.set_ylim(ylim[plot_style])
        
        # set up the figure
        self.ax.margins(x=10, y=1)
        self.ax.set_xlabel(xlabel[plot_style], fontsize='16', labelpad=1)
        self.ax.set_ylabel(ylabel[plot_style], fontsize='16', labelpad=1) 
        
        leg = plt.legend(loc="upper left", scatterpoints=1, fontsize=10)
        for i,l in enumerate(leg.legend_handles):
            leg.legend_handles[i]._sizes = [20]
        
        if plot_style == "arotscatter":
            self.ax.set_yscale('log')
        
        return self.artists

    def get_data(self, frame=0):
        d = self.ds.isel(time = frame)
        n=len(d['name'])
        d = d.isel(name=range(1,n))
        d['radmarker'] = (d['radius'] / self.Rcb) * self.radscale

        t = d['time'].values
        npl = d['npl'].values
        radmarker = d['radmarker'].values
        origin = d['origin_type'].values
        
        if plot_style == "aescatter":
            pl = np.c_[d['a'].values,d['e'].values]
        elif plot_style == "aiscatter":
            pl = np.c_[d['a'].values,d['inc'].values]
        elif plot_style == "arotscatter":
            pl = np.c_[d['a'].values,d['rot_mag'].values]

        return t, npl, pl, radmarker, origin

    def update(self,frame):
        """Update the scatter plot."""
        t,  npl, pl, radmarker, origin = self.get_data(framejump * frame)

        self.artists.title.set_text(f"{titletext} - Time = ${t*1e-6:6.3f}$ My with ${npl:4.0f}$ particles")

        for key,name in self.scatter_artist_names.items():
            idx = origin == key
            if any(idx) and any(~np.isnan(radmarker[idx])):
                scatter = self.artists._asdict()[name]
                scatter.set_sizes(radmarker[idx])
                scatter.set_offsets(pl[idx,:])
                scatter.set_facecolor(self.clist[key])

        return self.artists

sim = swiftest.Simulation(read_data=True, read_collisions=False, dask=True)
for plot_style in valid_plot_styles:
    animation_file = f"Chambers2013-{plot_style}.mp4"
    print('Making animation')
    anim = AnimatedScatter(sim.data,sim.param)
    print('Animation finished')
