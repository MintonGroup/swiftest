import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import Init_Cond as ic
from scipy.io import FortranFile
from uranian_satellites import *

# First set up the figure, the axis, and the plot element we want to animate


class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self):
        self.stream = self.data_stream()
        self.ringfilename = 'ring.dat'

        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots()
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1, frames=41000,
                                          init_func=self.setup_plot, blit=True)

        #self.ani.save('frames/uranian_ringsat.png', writer = "imagemagick")
        #self.ani.save('uranian_ringsat-S0.6e4g_cm2-480-dt10.mp4', fps=60, dpi=600, extra_args=['-vcodec', 'libx264'])

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        t, seeds, ring = next(self.stream)
        x = seeds[:,0]
        y = seeds[:,1]
        r = ring[:,0]
        s = ring[:,1]
        xmin = 1.0
        xmax = 6
        ymin = 0.01
        ymax = 5e4

        y2min = 1e11
        y2max = 1e24
        self.ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))

        #self.ax.set_xlim(xmin, xmax)
        #self.ax.set_ylim(ymin, ymax)
        self.ax.set_xlabel('Distance to Uranus (RU)')
        self.ax.set_ylabel('$\Sigma$ (g$\cdot$cm$^{-2}$)')
        self.ax.set_yscale('log')

        self.secax = self.ax.twinx()
        self.secax.set_yscale('log')
        self.secax.set_ylabel('Mass of satellite (g)')
        self.secax.set_ylim(y2min, y2max)

        self.line, = self.ax.plot(r, s, '-', color="black", linewidth=1.0, zorder=50)
        self.RRL = self.ax.plot([ic.RRL / ic.RP, ic.RRL / ic.RP], [ymin, ymax], '--', color="black", linewidth=0.5, zorder=50)
        self.RRLlab = self.ax.text(ic.RRL / ic.RP - 0.20, 0.4 * ymax, "RRL", rotation=90, fontsize="10")
        self.FRL = self.ax.plot([ic.FRL / ic.RP, ic.FRL / ic.RP], [ymin, ymax], ':', color="black", linewidth=0.5, zorder=50)
        self.FRLlab = self.ax.text(ic.FRL / ic.RP - 0.20, 0.4 * ymax, "FRL", rotation=90, fontsize="10")
        self.Rsync = self.ax.plot([ic.Rsync / ic.RP, ic.Rsync / ic.RP], [ymin, ymax], '-.', color="black", linewidth=0.5, zorder=50)
        self.Rsynclab = self.ax.text(ic.Rsync / ic.RP - 0.20, 0.4 * ymax, "$a_{sync}$", rotation=90, fontsize="10")
        #plt.axvline(x=xc, color='k', linestyle='--')

        self.title = self.ax.text(0.80, 0.1, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5},
                        transform=self.ax.transAxes, ha="center")
        #self.line.set_label(f'Time = ${t[0]*ic.TU2S/ic.year * 1e-6:5.1f}$ My')
        #self.legend = plt.legend()
        #self.legend.remove()
        self.title.set_text(f'Time = ${t[0] * ic.TU2S / ic.year * 1e-6:7.2f}$ My')
        self.usats = self.secax.scatter(Sat_r_RM, Sat_M_Mass, marker='o', color="silver", s=5, zorder=50)
        self.scat = self.secax.scatter(seeds[:,0], seeds[:,1], marker='o', color="black", s=5, zorder=50)

        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat, self.line, self.title,

    def data_stream(self):
        with FortranFile(self.ringfilename, 'r') as f:
            while True:
                try:
                    t = f.read_reals(np.float64)
                except:
                    f.close()
                    break
                Nbin = f.read_ints(np.int32)
                r = f.read_reals(np.float64)
                Gsigma = f.read_reals(np.float64)
                nu = f.read_reals(np.float64)
                kval = int(t / ic.t_print)
                Nseeds = f.read_ints(np.int32)
                a = f.read_reals(np.float64)
                Gm = f.read_reals(np.float64)

                yield t,np.c_[a / ic.RP, Gm * ic.MU2GM / ic.GU], np.c_[r / ic.RP, Gsigma * ic.MU2GM / ic.DU2CM**2 / ic.GU]


    def update(self, i):
        """Update the scatter plot."""
        t, seeds, ring = next(self.stream)
        #seeds = next(self.stream)

        # Set x and y data...
        r = ring[:,0]
        s = ring[:,1]
        self.scat.set_offsets(seeds[:, :2]) #data[:, :2])
        self.line.set_data(r, s) # :, :3])

        # Set sizes...
        #self.scat.set_sizes(300 * abs(data[:, 2])**1.5 + 100)
        # Set colors..
        #self.scat.set_array(x,y)

        self.title.set_text(f'Time = ${t[0]*ic.TU2S/ic.year * 1e-6:7.2f}$ My')
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat, self.line, self.title,


if __name__ == '__main__':
    anim = AnimatedScatter()

    plt.show()


