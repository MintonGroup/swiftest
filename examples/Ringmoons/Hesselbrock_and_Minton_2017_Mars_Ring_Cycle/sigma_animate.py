import sys

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
from tqdm import tqdm

import swiftest


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["cmbright"],
    "text.latex.preamble":r"\usepackage{cmbright}",
})

class AnimatedScatter:
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""

    def __init__(self):
        self.sim = swiftest.Simulation(read_data=True)
        nframes = self.sim.ring.time.size
        frames = tqdm(
            range(nframes),
            file=sys.stdout,
            desc="Creating animation",
            unit="frame",
            unit_scale=True,
            leave=False,
        )

        self.fig, self.ax = plt.subplots()

        xmin = 1.0
        xmax = 8
        ymin = 1.0
        ymax = 1.0e7

        y2min = 1e14
        y2max = 1e19

        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)
        self.ax.set_xlabel(r"Distance to Mars ($R_p$)", fontsize="12")
        self.ax.set_ylabel(r"$\Sigma$ (kg$\cdot$m$^{-2}$)", fontsize="12")
        self.ax.set_yscale("log")

        self.secax = self.ax.twinx()
        self.secax.set_yscale("log")
        self.secax.set_ylabel("Mass of satellite (kg)")
        self.secax.set_ylim(y2min, y2max)

        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1, frames=frames, init_func=self.init_func, blit=True)
        self.ani.save("Hesselbrock_and_Minton_2017_Mars_Ring_Cycle.mp4", fps=60, dpi=600, extra_args=["-vcodec", "libx264"])

    def init_func(self):
        """Initial drawing of the scatter plot."""
        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        (self.ring,) = self.ax.plot([], [], "-", color="black", linewidth=1.5, zorder=50)
        self.seeds = self.secax.scatter([], [], marker="o", color="black", s=5, zorder=50)

        self.title = self.ax.set_title("", fontsize="14", pad=20)

        return (
            self.ring,
            self.seeds,
            self.title,
        )

    def data_stream(self, frame=0):
        while True:
            ring = self.sim.ring.isel(time=frame)
            seeds = self.sim.data.isel(time=frame)
            seeds = seeds.where(seeds.particle_type=="Ringmoons Seed",drop=True)
            seeds = seeds.where(~np.isnan(seeds.a),drop=True)
            t = ring.time.data[()]
            r = ring.r.data
            sigma = ring.sigma.data * self.sim.MU2KG / self.sim.DU2M**2
            aseed = seeds.a.data
            mseed = seeds.mass.data * self.sim.MU2KG
            yield t, r, sigma, aseed, mseed

    def update(self, frame):
        """Update the scatter plot."""
        t, r, sigma, aseed, mseed = next(self.data_stream(frame))

        # Set x and y data with unit conversion from m to 1000 km and from kg/m^2 to 10^4 kg/m^2, respectively.
        self.ring.set_data(r, sigma)
        if len(aseed) > 0:
            self.seeds.set_offsets(np.c_[aseed,mseed])

        self.title.set_text(f"Time = ${t * 1e-6:3.3f}$ My")
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return (
            self.ring,
            self.seeds,
            self.title,
        )


if __name__ == "__main__":
    anim = AnimatedScatter()
