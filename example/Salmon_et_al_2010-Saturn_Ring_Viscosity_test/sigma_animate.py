import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from Init_Cond import *
from scipy.io import FortranFile


MU2KG = MU2GM * 1e-3
DU2KMe3 = DU2CM * 1e-8


SIGU2SI = 10 * MU2GM / DU2CM**2
SIGU2CGS = MU2GM / DU2CM**2
class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self):
        self.stream = self.data_stream()
        self.ringfilename = 'ring.dat'

        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots()
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1, frames=500,
                                          init_func=self.setup_plot, blit=True)

        #self.ani.save('frames/charnoz2010-saturn-ringmoons.png', writer = "imagemagick")
        self.ani.save('salmon2010-saturn-viscosity.mp4', fps=60, dpi=600, extra_args=['-vcodec', 'libx264'])

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        #t, seeds, ring = next(self.stream)
        t = np.array([0.0])
        ring = np.array([[0,0]])

        r = ring[:,0] / RP
        s = ring[:,1] * SIGU2CGS
        xmin = 1.7
        xmax = 2.1
        ymin = 0.0
        ymax = 7000.0



        self.ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))

        #self.ax.set_xlim(xmin, xmax)
        #self.ax.set_ylim(ymin, ymax)
        self.ax.set_xlabel('Distance to Saturn ($R_p$)', fontsize='12')
        self.ax.set_ylabel('$\Sigma$ (g$\cdot$cm$^{-2}$)', fontsize='12')
        #self.ax.set_yscale('log')

        self.line, = self.ax.plot(r, s, '-', color="black", linewidth=1.5, zorder=50)
        #plt.axvline(x=xc, color='k', linestyle='--')

        self.title = self.ax.text(0.80, 0.1, "", bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5},
                        transform=self.ax.transAxes, ha="center")
        #self.line.set_label(f'Time = ${t[0]*TU2S/year * 1e-6:5.1f}$ My')
        #self.legend = plt.legend()
        #self.legend.remove()
        self.title.set_text(f'Time = ${t[0] * TU2S / year * 1e-3:4.0f}$ ky')

        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.line, self.title,

    def data_stream(self):
        with FortranFile(self.ringfilename, 'r') as f:
            while True:
                for _ in range(1):
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

        # Set x and y data...
        r = ring[:,0] / RP
        s = ring[:,1] * SIGU2CGS
        self.line.set_data(r, s)

        # Set sizes...

        self.title.set_text(f'Time = ${t[0]*TU2S/year * 1e-3:4.0f}$ ky')
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.line, self.title,


if __name__ == '__main__':
    anim = AnimatedScatter()

    plt.show()


