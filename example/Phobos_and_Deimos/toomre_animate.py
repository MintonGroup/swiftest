import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import Init_Cond as ic
from scipy.io import FortranFile
from Init_Cond import *

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

        #self.ani.save('frames/uranian_ringsat.png', writer = "imagemag")
        #self.ani.save('uranian_ringsat-S0.6e4g_cm2-480-dt10.mp4', fps=60, dpi=600, extra_args=['-vcodec', 'libx264'])

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        t, seeds, ring = next(self.stream)
        x = seeds[:,0]
        y = seeds[:,1]
        r = ring[:,0]
        s = ring[:,1]
        nu = ring[:,2]
        Q = ring[:,3]
        r_pdisk = ring[:,4]
        xmin = 1.0
        xmax = 4.0
        ymin = 1.0e-2
        ymax = 1e4

        y2min = 1e0
        y2max = 1e6
        self.ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))

        #self.ax.set_xlim(xmin, xmax)
        #self.ax.set_ylim(ymin, ymax)
        self.ax.set_xlabel('Distance to Mars (Rp)')
        self.ax.set_ylabel('$r_{pdisk}$ (cm)')
        self.ax.set_yscale('log')

        self.secax = self.ax.twinx()
        #self.secax.set_ylabel('Toomre parameter Q', color="blue")
        self.secax.set_ylabel('Kinematic viscosity (cm$^2$ s$^{-1}$)', color="blue")
        self.secax.set_yscale('log')
        self.secax.set_ylim(y2min, y2max)

        self.rpline, = self.ax.plot(r, r_pdisk, '-', color="black", linewidth=1.0, zorder=50)
        self.RRL = self.ax.plot([RRL / RP, RRL / RP], [ymin, ymax], '--', color="black", linewidth=0.5, zorder=50)
        self.RRLlab = self.ax.text(RRL / RP - 0.20, 0.4 * ymax, "RRL", rotation=90, fontsize="10")
        self.FRL = self.ax.plot([FRL / RP, FRL / RP], [ymin, ymax], ':', color="black", linewidth=0.5, zorder=50)
        self.FRLlab = self.ax.text(FRL / RP - 0.20, 0.4 * ymax, "FRL", rotation=90, fontsize="10")
        self.Rsync = self.ax.plot([Rsync / RP, Rsync / RP], [ymin, ymax], '-.', color="black", linewidth=0.5, zorder=50)
        self.Rsynclab = self.ax.text(Rsync / RP - 0.20, 0.4 * ymax, "$a_{sync}$", rotation=90, fontsize="10")
        #plt.axvline(x=xc, color='k', linestyle='--')
        #self.Qstab = self.secax.plot(r,np.full_like(r,1.0), ':', color="blue", linewidth=0.5, zorder = 50)
        #self.Qstab = self.secax.plot(r,np.full_like(r,2.0), '-.', color="blue", linewidth=0.5, zorder = 50)

        self.Qline, = self.secax.plot(r, nu, '-', color="blue", linewidth=1.0, zorder=50)
        self.title = self.ax.text(0.80, 0.1, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5},
                        transform=self.ax.transAxes, ha="center")
        #self.line.set_label(f'Time = ${t[0]*TU2S/year * 1e-6:5.1f}$ My')
        #self.legend = plt.legend()
        #self.legend.remove()
        self.title.set_text(f'Time = ${t[0] * TU2S / year * 1e-6:7.2f}$ My')

        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.Qline, self.rpline, self.title,

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
                Q = f.read_reals(np.float64)
                r_pdisk = f.read_reals(np.float64)
                vrel_pdisk = f.read_reals(np.float64)
                kval = int(t / t_print)
                Nseeds = f.read_ints(np.int32)
                a = f.read_reals(np.float64)
                Gm = f.read_reals(np.float64)

                yield t,np.c_[a / RP,
                              Gm * MU2GM / GU],\
                        np.c_[r / RP,
                              Gsigma * MU2GM / DU2CM**2 / GU,
                              nu * DU2CM**2 / TU2S,
                              Q,
                              r_pdisk * DU2CM,
                              vrel_pdisk * DU2CM / TU2S]


    def update(self, i):
        """Update the scatter plot."""
        t, seeds, ring = next(self.stream)
        #seeds = next(self.stream)

        # Set x and y data...
        r = ring[:,0]
        s = ring[:,1]
        nu = ring[:,2]
        Q = ring[:,3]
        r_pdisk = ring[:,4]
        r_pdisk[s < 1e-6] = 0.0
        Q[s < 1e-6] = 0.0
        nu[s < 1e-8] = 0.0

        self.rpline.set_data(r, r_pdisk)

        self.Qline.set_data(r,nu)

        # Set sizes...
        #self.scat.set_sizes(300 * abs(data[:, 2])**1.5 + 100)
        # Set colors..
        #self.scat.set_array(x,y)

        self.title.set_text(f'Time = ${t[0]*TU2S/year * 1e-6:7.2f}$ My')
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.Qline, self.rpline, self.title,


if __name__ == '__main__':
    anim = AnimatedScatter()

    plt.show()


