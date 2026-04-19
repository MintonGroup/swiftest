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


cb = swiftest.get_solar_system_body("Saturn")
r_cb = cb["radius"]
m_cb = cb["mass"]
gm_cb = cb["Gmass"]
rho_cb = m_cb / (4.0 / 3.0 * np.pi * r_cb**3)
rot_cb = np.sqrt(cb["rot"][0] ** 2 + cb["rot"][1] ** 2 + cb["rot"][2] ** 2)
sim = swiftest.Simulation(integrator="ringmoons", MU="kg", DU="m", TU="yr")
sim.add_solar_system_body(["Saturn"], align_to_central_body_rotation=True)

r_p = 1.0e-2 * sim.M2DU  # disk particle size
rho_p = 900.0 * sim.KG2MU / sim.M2DU**3  # disk particle density
m_p = 4.0 / 3.0 * np.pi * r_p**3 * rho_p  # disk particle mass

# Convert to simualtion units
r_cb *= sim.M2DU
m_cb *= sim.KG2MU
gm_cb *= sim.M2DU**3 / sim.S2TU**2
rho_cb *= sim.KG2MU / sim.M2DU**3
rot_cb = 360.0 / rot_cb * sim.S2TU

mass_distribution = {
    "type": "gaussian",
    "mu" : 110.0e6,   # radius of the center of the gaussian (m)
    "dev" : 360.0e3,  # width of the gaussian (m)
    "sigma0" : 6.15e4,  # peak of the surface mass density at the center of the gaussian (kg/m^2)
    "nbins" : 2048
}

sim.add_ring(
    r_p=r_p,
    m_p=m_p,
    mass_distribution=mass_distribution
)

dt = 100.0
istep_out = 1
tstop = 1e5

sim.set_parameter(tstop=tstop, dt=dt, istep_out=istep_out, dump_cadence=0)
sim.run()

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

        xmin = 100.0
        xmax = 120.0
        ymin = 0.0
        ymax = 8.0

        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)
        self.ax.set_xticks(np.arange(xmin, xmax+5, 5))
        self.ax.set_xticks(np.arange(xmin, xmax+1, 1), minor=True)
        self.ax.set_yticks(np.arange(ymin, ymax+2, 2))
        self.ax.set_yticks(np.arange(ymin, ymax+0.5, 0.5), minor=True)
        self.ax.set_xlabel("Distance to Saturn (1000 km)", fontsize="12")
        self.ax.set_ylabel(r"$\Sigma$ ($10^4$kg$\cdot$m$^{-2}$)", fontsize="12")


        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1, frames=frames, init_func=self.init_func, blit=True)
        self.ani.save("salmon2010-saturn-viscosity.mp4", fps=60, dpi=300, extra_args=["-vcodec", "libx264"])

    def init_func(self):
        """Initial drawing of the scatter plot."""
        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        (self.line,) = self.ax.plot([], [], "-", color="black", linewidth=1.5, zorder=50)

        self.title = self.ax.set_title("", fontsize="14", pad=20)

        return (
            self.line,
            self.title,
        )

    def data_stream(self, frame=0):
        while True:
            ds = self.sim.ring.isel(time=frame)
            t = ds.time.data[()]
            r = ds.r.data
            sigma = ds.sigma.data
            yield t, r, sigma

    def update(self, frame):
        """Update the scatter plot."""
        t, r, sigma = next(self.data_stream(frame))

        # Set x and y data with unit conversion from m to 1000 km and from kg/m^2 to 10^4 kg/m^2, respectively.
        self.line.set_data(r*1e-6, sigma*1e-4)

        self.title.set_text(f"Time = ${t * 1e-3:4.0f}$ ky")
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return (
            self.line,
            self.title,
        )


anim = AnimatedScatter()

