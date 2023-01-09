#!/usr/bin/env python3
import swiftest 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.colors as mcolors

titletext = "Chambers (2013)"
valid_plot_styles = ["aescatter", "aiscatter"]
xlim={"aescatter" : (0.0, 2.25),
      "aiscatter" : (0.0, 2.25)}
ylim={"aescatter" : (0.0, 1.0),
      "aiscatter" : (0.0, 40.0)}
xlabel={"aescatter": "Semimajor axis (AU)",
        "aiscatter": "Semimajor axis (AU)"}
ylabel={"aescatter": "Eccentricity",
        "aiscatter": "Inclination (deg)"}


plot_style = valid_plot_styles[1]
framejump = 1
animation_file = f"Chambers2013-{plot_style}.mp4"

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, ds, param):

        self.radscale = 2000
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
        self.ax.set_xlim(xlim[plot_style])
        self.ax.set_ylim(ylim[plot_style])
        fig.add_axes(self.ax)
        self.ani = animation.FuncAnimation(fig, self.update, interval=1, frames=nframes, init_func=self.setup_plot, blit=True)
        self.ani.save(animation_file, fps=60, dpi=300, extra_args=['-vcodec', 'libx264'])
        print(f'Finished writing {animation_file}')

    def scatters(self, pl, radmarker, origin):
        scat = []
        for key, value in self.clist.items():
            idx = origin == key
            s = self.ax.scatter(pl[idx, 0], pl[idx, 1], marker='o', s=radmarker[idx], c=value, alpha=0.75, label=key, animated=True)
            scat.append(s)
        return scat

    def setup_plot(self):
        # First frame
        """Initial drawing of the scatter plot."""
        t, npl, pl, radmarker, origin = next(self.data_stream(0))

        # set up the figure
        self.ax.margins(x=10, y=1)
        self.ax.set_xlabel(xlabel[plot_style], fontsize='16', labelpad=1)
        self.ax.set_ylabel(ylabel[plot_style], fontsize='16', labelpad=1)

        title = self.ax.text(0.50, 1.05, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5}, transform=self.ax.transAxes,
                        ha="center", animated=True)

        title.set_text(f"{titletext} -  Time = ${t*1e-6:6.2f}$ My with ${npl:4.0f}$ particles")
        self.slist = self.scatters(pl, radmarker, origin)
        self.slist.append(title)
        leg = plt.legend(loc="upper left", scatterpoints=1, fontsize=10)
        for i,l in enumerate(leg.legendHandles):
           leg.legendHandles[i]._sizes = [20]
        return self.slist

    def data_stream(self, frame=0):
        while True:
            d = self.ds.isel(time = frame)
            name_good = d['name'].where(d['name'] != "Sun", drop=True)
            d = d.sel(name=name_good)
            d['radmarker'] = (d['radius'] / self.Rcb) * self.radscale

            t = d['time'].values
            npl = d['npl'].values
            radmarker = d['radmarker'].values
            origin = d['origin_type'].values
            
            if plot_style == "aescatter":
               pl = np.c_[d['a'].values,d['e'].values]
            elif plot_style == "aiscatter":
               pl = np.c_[d['a'].values,d['inc'].values]

            yield t, npl, pl, radmarker, origin

    def update(self,frame):
        """Update the scatter plot."""
        t,  npl, pl, radmarker, origin = next(self.data_stream(framejump * frame))

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
