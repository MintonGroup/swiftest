import numpy as np
from matplotlib import pyplot as plt
import xarray as xr

import swiftest

cb = swiftest.get_solar_system_body("Uranus")
r_cb = cb["radius"]
m_cb = cb["mass"]
gm_cb = cb["Gmass"]
rho_cb = m_cb / (4.0 / 3.0 * np.pi * r_cb**3)
rot_cb = np.sqrt(cb["rot"][0] ** 2 + cb["rot"][1] ** 2 + cb["rot"][2] ** 2)

sim = swiftest.Simulation(integrator="ringmoons", MU2KG=m_cb, DU2M=r_cb, TU="yr")
sim.clean()
sim.add_solar_system_body("Uranus", align_to_central_body_rotation=True)
sim.modify_body(name="Uranus",k2=0.104,Q=3000.0)

# Convert to simualtion units
r_p = 1.0 * sim.M2DU  # disk particle size
rho_p = 1200.0 * sim.KG2MU / sim.M2DU**3  # disk particle density
m_p = 4.0 / 3.0 * np.pi * r_p**3 * rho_p  # disk particle mass
r_cb *= sim.M2DU
m_cb *= sim.KG2MU
gm_cb *= sim.M2DU**3 / sim.S2TU**2
rho_cb *= sim.KG2MU / sim.M2DU**3
rot_cb = 360.0 / rot_cb * sim.S2TU

# Compute useful radii
frl = 2.456 * r_cb * (rho_cb / rho_p) ** (1.0 / 3.0)
rrl = 1.44 * r_cb * (rho_cb / rho_p) ** (1.0 / 3.0)
rsync = (gm_cb * rot_cb**2 / (4 * np.pi**2)) ** (1.0 / 3.0)
alind = 4 ** (1.0 / 3.0) * frl
sigma_frl = 8000.0 * sim.KG2MU / sim.M2DU**2
alpha = -3.0
sigma0 = sigma_frl * (frl)**(-alpha)

sim.add_ring(
    r_p=r_p,
    m_p=m_p,
    mass_distribution={
        "type": "powerlaw",
        "sigma0": sigma0,
        "alpha": alpha,
        "nbins": 1024,
        "r_outer": 1.2*frl,
    },
)
dt = 1e3
tstop = 100*dt
sim.ring["sigma"] = xr.where(sim.ring.r < frl, sim.ring.sigma, xr.zeros_like(sim.ring.sigma))


sim.set_parameter(tstop=tstop, dt=dt, tstep_out=dt, dump_cadence=1)
sim.save()
