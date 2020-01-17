import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from Init_Cond import *
from scipy.io import FortranFile
import matplotlib.patches as patches



alind = 4**(1./3.) * FRL

SIGU2CGS = MU2GM / DU2CM**2
class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self):
        self.ringfilename = 'ring.dat'
        self.frame_skip = 10
        nframes = int(74140 / self.frame_skip)

        self.stream = self.data_stream()
        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots(figsize=(6,2.33))
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1, frames=nframes,
                                          init_func=self.setup_plot, blit=True)


        #self.ani.save('frames/charnoz2010-saturn-ringmoons.png', writer = "imagemagick")
        self.ani.save('HesselbrockMinton2017PhobosCycle-thumbnail.mp4', fps=60, dpi=300, extra_args=['-vcodec', 'libx264'])

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        t = np.array([0])
        rs = [0.0]
        ms = [0.0]
        r = [0.0]
        s = [0.0]

        xmin = 1.0
        xmax = 8.00 #1.5 * FRL  / RP
        ymin = 0.1
        ymax = 1e6

        y2min = 1e17
        y2max = 1e22

        tsize = 14
        RRLcolor = 'sandybrown'
        FRLcolor = 'slategrey'
        asynccolor = 'maroon'
        alindcolor = 'steelblue'
        satcolor = 'cadetblue'
        seedcolor = 'black'
        sigmacolor = 'black'
        satalpha = 1.0

        self.ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))

        self.ax.set_xlabel('Distance to Mars ($R_p$)', fontsize=tsize)
        self.ax.set_ylabel('$\Sigma$ (g$\cdot$cm$^{-2}$)', fontsize=tsize)
        self.ax.set_yscale('log')

        self.secax = self.ax.twinx()
        self.secax.set_yscale('log')
        self.secax.set_ylabel('$M_{satellite}$ (g)', fontsize=tsize)
        self.secax.set_ylim(y2min, y2max)

        self.ax.tick_params(axis='both', which='major', labelsize=tsize)
        self.secax.tick_params(axis='both', which='major', labelsize=tsize)

        # Time label
        #self.title = self.ax.text(0.82, 0.1, "", bbox={'facecolor': 'w', 'alpha': 1.0, 'pad': 5},
        #                transform=self.ax.transAxes, ha='center', fontsize=tsize)
        #self.title.set_text(f'Time = ${t[0] * TU2S / year * 1e-6:4.0f}$ My')

        # surface mass density and seeds
        self.line, = self.ax.plot(r, s, '-', color=sigmacolor, linewidth=1.5, zorder=50)
        self.scat = self.secax.scatter(rs, ms, marker='o', color="black", s=25, zorder=50)

        # reference lines
        self.RRL = self.ax.plot([RRL / RP, RRL / RP], [ymin, ymax], '--', color=RRLcolor, linewidth=1.5, zorder=40)
        self.RRLlab = self.ax.text(RRL / RP - 0.00, 1.5 * ymax, "RRL", color=RRLcolor, rotation=0, fontsize=tsize, ha='center')
        self.FRL = self.ax.plot([FRL / RP , FRL / RP], [ymin, ymax], ':', color=FRLcolor, linewidth=1.5, zorder=40)
        self.FRLlab = self.ax.text(FRL / RP - 0.00, 1.5 * ymax, "FRL", color=FRLcolor, rotation=0, fontsize=tsize, ha='center')
        self.Rsync = self.ax.plot([Rsync / RP , Rsync / RP], [ymin, ymax], '-.', color=asynccolor, linewidth=1.5, zorder=40)
        self.Rsynclab = self.ax.text(Rsync / RP + 0.10, 1.8 * ymax, "$a_{sync}$", color=asynccolor, rotation=0, fontsize=tsize, ha='center')
        # self.Rlind = self.ax.plot([alind / RP , alind / RP], [ymin, ymax], '-.', color=alindcolor, linewidth=1.5, zorder=40)
        # self.Rlindlab = self.ax.text(alind / RP - 0.02, 1.5 * ymax, "$a_{lind}$", color=alindcolor,rotation=0, fontsize=tsize, ha='center')
        #plt.axvline(x=xc, color='k', linestyle='--')

        xsc = Sat_semimajor / DU2CM / RP
        satmask = (xsc >= xmin) & (xsc <= xmax) & (Sat_mass >= y2min) & (Sat_mass <= y2max)
        satx = xsc[satmask]
        saty = Sat_mass[satmask]
        satn = Sat_name[satmask]

        self.fig.tight_layout()

        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat, self.line #, self.title,

    def data_stream(self):
        with FortranFile(self.ringfilename, 'r') as f:
            while True:
                for _ in range(self.frame_skip):
                    try:
                        t = f.read_reals(np.float64)
                    except:
                        f.close()
                        break
                    Nbin = f.read_ints(np.int32)
                    r = f.read_reals(np.float64)
                    Gsigma = f.read_reals(np.float64)
                    nu = f.read_reals(np.float64)
                    Q = f.read_reals(np.float64)
                    r_pdisk = f.read_reals(np.float64)
                    vrel_pdisk = f.read_reals(np.float64)
                    Nseeds = f.read_ints(np.int32)
                    a = f.read_reals(np.float64)
                    Gm = f.read_reals(np.float64)

                yield t,np.c_[a , Gm / GU], np.c_[r , Gsigma / GU]


    def update(self, i):
        """Update the scatter plot."""
        t, seeds, ring = next(self.stream)


        seeds[:,0] = seeds[:,0] / RP
        seeds[:,1] = seeds[:,1] * MU2GM
        # Set x and y data...
        r = ring[:,0] / RP
        s = ring[:,1] * SIGU2CGS
        self.scat.set_offsets(seeds)
        self.line.set_data(r, s)

        # Set sizes...
        #self.scat.set_sizes(300 * abs(data[:, 2])**1.5 + 100)
        # Set colors..
        #self.scat.set_array(x,y)

        #self.title.set_text(f'Time = ${t[0]*TU2S/year * 1e-6:4.0f}$ My')
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat, self.line #, self.title,


if __name__ == '__main__':
    anim = AnimatedScatter()

    plt.show()


