#!/usr/bin/env python3
import swiftest 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.colors as mcolors

titletext = "Chambers (2013)"
radscale = 2000
AU = 1.0
xmin = 0.0
xmax = 2.25
ymin = 0.0
ymax = 1.0
framejump = 10
ncutoff = 1e20

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, ds, param):

        frame = 0
        nframes = int(ds['time'].size / framejump)
        self.ds = ds
        self.param = param
        self.Rcb = self.ds['radius'].sel(name="Sun").isel(time=0).values[()]

        self.clist = {'Initial conditions' : 'k',
                      'Merger' : 'xkcd:faded blue',
                      'Disruption' : 'xkcd:marigold',
                      'Supercatastrophic' : 'xkcd:shocking pink',
                      'Hit and run fragmentation' : 'xkcd:baby poop green'}

        # Setup the figure and axes...
        fig = plt.figure(figsize=(8,4.5), dpi=300)
        plt.tight_layout(pad=0)
        # set up the figure
        self.ax = plt.Axes(fig, [0.1, 0.15, 0.8, 0.75])
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)
        fig.add_axes(self.ax)
        self.ani = animation.FuncAnimation(fig, self.update, interval=1, frames=nframes, init_func=self.setup_plot, blit=True)
        self.ani.save('aescatter.mp4', fps=30, dpi=300, extra_args=['-vcodec', 'libx264'])
        print('Finished writing aescattter.mp4')

    def scatters(self, pl, radmarker, origin):
        scat = []
        for key, value in self.clist.items():
            idx = origin == key
            s = self.ax.scatter(pl[idx, 0], pl[idx, 1], marker='o', s=radmarker[idx], c=value, alpha=0.75, label=key)
            scat.append(s)
        return scat

    def setup_plot(self):
        # First frame
        """Initial drawing of the scatter plot."""
        t, name, Gmass, radius, npl, pl, radmarker, origin = next(self.data_stream(0))

        # set up the figure
        self.ax.margins(x=10, y=1)
        self.ax.set_xlabel("Semimajor Axis (AU)", fontsize='16', labelpad=1)
        self.ax.set_ylabel("Eccentricity", fontsize='16', labelpad=1)

        title = self.ax.text(0.50, 1.05, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5}, transform=self.ax.transAxes,
                        ha="center")

        title.set_text(f"{titletext} -  Time = ${t*1e-6:6.2f}$ My with ${npl:4.0f}$ particles")
        self.slist = self.scatters(pl, radmarker, origin)
        self.slist.append(title)
        leg = plt.legend(loc="upper right", scatterpoints=1, fontsize=10)
        for i,l in enumerate(leg.legendHandles):
           leg.legendHandles[i]._sizes = [20]
        return self.slist

    def data_stream(self, frame=0):
        while True:
            d = self.ds.isel(time = frame)
            name_good = d['name'].where(d['name'] != "Sun", drop=True)
            d = d.sel(name=name_good)
            d['radmarker'] = (d['radius'] / self.Rcb) * radscale
            radius = d['radmarker'].values

            radius = d['radmarker'].values
            Gmass = d['Gmass'].values
            a = d['a'].values 
            e = d['e'].values
            name = d['name'].values
            npl = d['npl'].values
            radmarker = d['radmarker']
            origin = d['origin_type']

            t = self.ds.coords['time'].values[frame]

            yield t, name, Gmass, radius, npl, np.c_[a, e], radmarker, origin

    def update(self,frame):
        """Update the scatter plot."""
        t, name, Gmass, radius, npl, pl, radmarker, origin = next(self.data_stream(framejump * frame))

        self.slist[-1].set_text(f"{titletext} - Time = ${t*1e-6:6.3f}$ My with ${npl:4.0f}$ particles")

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        for i, (key, value) in enumerate(self.clist.items()):
            idx = origin == key
            self.slist[i].set_sizes(radmarker[idx])
            self.slist[i].set_offsets(pl[idx,:])
            self.slist[i].set_facecolor(value)

        return self.slist

sim = swiftest.Simulation(read_old_output=True)
print('Making animation')
anim = AnimatedScatter(sim.data,sim.param)
print('Animation finished')
