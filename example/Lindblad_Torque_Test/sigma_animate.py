import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from Init_Cond import *
from scipy.io import FortranFile



alind = 4**(1./3.) * FRL

SIGU2CGS = MU2GM / DU2CM**2
class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self):
        self.ringfilename = 'ring.dat'
        self.stream = self.data_stream()

        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots()
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1, frames=1,
                                          init_func=self.setup_plot, blit=True)

        #self.ani.save('frames/charnoz2010-saturn-ringmoons.png', writer = "imagemagick")
        #self.ani.save('LindbladTorqueTest.mp4', fps=60, dpi=600, extra_args=['-vcodec', 'libx264'])

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        t, seeds, ring = next(self.stream)
        rs = seeds[:,0] / RP
        ms = seeds[:,1] * MU2GM
        r = ring[:,0] / RP
        s = ring[:,1] * SIGU2CGS
        xmin = 1.0
        xmax = 8.00 #1.5 * FRL  / RP
        ymin = 0.1
        ymax = 1e5

        y2min = 1e16
        y2max = 1e24


        self.ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))

        #self.ax.set_xlim(xmin, xmax)
        #self.ax.set_ylim(ymin, ymax)
        self.ax.set_xlabel('Distance to Mars ($R_p$)', fontsize='12')
        self.ax.set_ylabel('$\Sigma$ (g$\cdot$cm$^{-2}$)', fontsize='12')
        self.ax.set_yscale('log')

        self.secax = self.ax.twinx()
        self.secax.set_yscale('log')
        self.secax.set_ylabel('Mass of satellite (g)', fontsize='12')
        self.secax.set_ylim(y2min, y2max)

        self.line, = self.ax.plot(r, s, '-', color="black", linewidth=1.5, zorder=50)
        self.RRL = self.ax.plot([RRL / RP, RRL / RP], [ymin, ymax], '--', color="black", linewidth=0.5, zorder=50)
        self.RRLlab = self.ax.text(RRL / RP - 0.00, 0.3 * ymax, "RRL", rotation=90, fontsize="12")
        self.FRL = self.ax.plot([FRL / RP , FRL / RP], [ymin, ymax], ':', color="black", linewidth=0.5, zorder=50)
        self.FRLlab = self.ax.text(FRL / RP - 0.00, 0.3 * ymax, "FRL", rotation=90, fontsize="12")
        self.Rsync = self.ax.plot([Rsync / RP , Rsync / RP], [ymin, ymax], '-.', color="black", linewidth=0.5, zorder=50)
        self.Rsynclab = self.ax.text(Rsync / RP + 0.02, 0.25 * ymax, "$a_{sync}$", rotation=90, fontsize="12")
        self.Rlind = self.ax.plot([alind / RP , alind / RP], [ymin, ymax], '-.', color="black", linewidth=0.5, zorder=50)
        self.Rlindlab = self.ax.text(alind / RP + 0.02, 0.25 * ymax, "$a_{lind}$", rotation=90, fontsize="12")
        #plt.axvline(x=xc, color='k', linestyle='--')

        self.title = self.ax.text(0.80, 0.1, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5},
                        transform=self.ax.transAxes, ha="center")
        #self.line.set_label(f'Time = ${t[0]*TU2S/year * 1e-6:5.1f}$ My')
        #self.legend = plt.legend()
        #self.legend.remove()
        self.title.set_text(f'Time = ${t[0] * TU2S / year * 1e-3:4.0f}$ ky')
        self.usats = self.secax.scatter(Sat_r_RM / RP / DU2CM , Sat_M_Mass, marker='o', color="lightsteelblue", s=15, zorder=50)
        self.scat = self.secax.scatter(rs, ms, marker='o', color="black", s=15, zorder=100)

        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat, self.line, self.title,

    def data_stream(self):
        with FortranFile(self.ringfilename, 'r') as f:
            while True:
                #for _ in range(1):
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
        self.scat.set_offsets(seeds[: :2])
        self.line.set_data(r, s)



        # Set sizes...
        #self.scat.set_sizes(300 * abs(data[:, 2])**1.5 + 100)
        # Set colors..
        #self.scat.set_array(x,y)

        self.title.set_text(f'Time = ${t[0]*TU2S/year * 1e-3:4.0f}$ ky')
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat, self.line, self.title,


if __name__ == '__main__':
    anim = AnimatedScatter()

    plt.show()


