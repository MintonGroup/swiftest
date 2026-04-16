import numpy as np
from matplotlib import pyplot as plt

import swiftest

cb = swiftest.get_solar_system_body("Mars")
r_cb = cb["radius"]
m_cb = cb["mass"]
gm_cb = cb["Gmass"]
rho_cb = m_cb / (4.0 / 3.0 * np.pi * r_cb**3)
rot_cb = np.sqrt(cb["rot"][0] ** 2 + cb["rot"][1] ** 2 + cb["rot"][2] ** 2)
sim = swiftest.Simulation(integrator="ringmoons", MU2KG=m_cb, DU2M=r_cb, TU="yr")
sim.add_solar_system_body(["Mars", "Phobos", "Deimos"], align_to_central_body_rotation=True)

r_p = 25.0e-2 * sim.M2DU  # disk particle size
rho_p = 1500.0 * sim.KG2MU / sim.M2DU**3  # disk particle density
m_p = 4.0 / 3.0 * np.pi * r_p**3 * rho_p  # disk particle mass

# Convert to simualtion units
r_cb *= sim.M2DU
m_cb *= sim.KG2MU
gm_cb *= sim.M2DU**3 / sim.S2TU**2
rho_cb *= sim.KG2MU / sim.M2DU**3
rot_cb = 360.0 / rot_cb * sim.S2TU

sim.add_ring(
    r_p=r_p,
    m_p=m_p,
    mass_distribution={
        "type": "powerlaw",
        "sigma0": 1.05e7 * sim.KG2MU / sim.M2DU**2,
        "alpha": -8.0,
        "nbins": 256,
    },
)

tstep_out = 1.0e6
dt = 1.0
tstop = 1.0e10


sim.set_parameter(tstop=tstop, dt=dt, tstep_out=tstep_out, dump_cadence=1)
sim.save()

xmin = 1.0
xmax = 8.00
ymin = 0.1
ymax = 1e6

y2min = 1e17
y2max = 1e22

tsize = 16
rrlcolor = "sandybrown"
frlcolor = "slategrey"
asynccolor = "maroon"
alindcolor = "steelblue"
satcolor = "cadetblue"
seedcolor = "black"
sigmacolor = "black"
satalpha = 1.0

fig, ax = plt.subplots(layout="constrained")
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

ax.set_xlabel("Distance to Mars ($R_p$)", fontsize=tsize)
ax.set_ylabel(r"Ring surface mass density (g$\cdot$cm$^{-2}$)", fontsize=tsize)
ax.set_yscale("log")

secax = ax.twinx()
secax.set_yscale("log")
secax.set_ylabel("Mass of satellite (g)", fontsize=tsize)
secax.set_ylim(y2min, y2max)

ax.tick_params(axis="both", which="major", labelsize=tsize)
secax.tick_params(axis="both", which="major", labelsize=tsize)

ids = sim.data.isel(time=0)
r = ids.ring_r.values / r_cb
s = ids.ring_sigma.values * sim.MU2KG / sim.DU2M**2 * 1000.0 / 100.0

rs = ids.a.values[1:] / r_cb
ms = ids.mass.values[1:] * sim.MU2KG * 1000.0

frl = 2.456 * r_cb * (rho_cb / rho_p) ** (1.0 / 3.0)
rrl = 1.44 * r_cb * (rho_cb / rho_p) ** (1.0 / 3.0)
rsync = (gm_cb * rot_cb**2 / (4 * np.pi**2)) ** (1.0 / 3.0)
alind = 4 ** (1.0 / 3.0) * frl

# surface mass density and seeds
ax.plot(r, s, "-", color=sigmacolor, linewidth=1.5, zorder=50)
scat = secax.scatter(rs, ms, marker="o", color="black", s=25, zorder=20)

# reference lines
ax.plot([rrl / r_cb, rrl / r_cb], [ymin, ymax], "--", color=rrlcolor, linewidth=1.5, zorder=50)
ax.text(rrl / r_cb - 0.00, 1.3 * ymax, "RRL", color=rrlcolor, rotation=0, fontsize=tsize, ha="center")
ax.plot([frl / r_cb, frl / r_cb], [ymin, ymax], ":", color=frlcolor, linewidth=1.5, zorder=50)
ax.text(frl / r_cb - 0.00, 1.3 * ymax, "FRL", color=frlcolor, rotation=0, fontsize=tsize, ha="center")
ax.plot([rsync / r_cb, rsync / r_cb], [ymin, ymax], "-.", color=asynccolor, linewidth=1.5, zorder=50)
ax.text(rsync / r_cb + 0.10, 1.3 * ymax, "$a_{sync}$", color=asynccolor, rotation=0, fontsize=tsize, ha="center")
ax.plot([alind / r_cb, alind / r_cb], [ymin, ymax], "-.", color=alindcolor, linewidth=1.5, zorder=50)
ax.text(alind / r_cb - 0.02, 1.3 * ymax, "$a_{lind}$", color=alindcolor, rotation=0, fontsize=tsize, ha="center")

plt.show()
