import numpy as np
from matplotlib import pyplot as plt
import xarray as xr

import swiftest

cb = swiftest.get_solar_system_body("Mars")
r_cb = cb["radius"]
m_cb = cb["mass"]

sim = swiftest.Simulation(integrator="ringmoons", MU2KG=m_cb, DU2M=r_cb, TU="yr")
sim.clean()
sim.add_solar_system_body("Mars", align_to_central_body_rotation=True)
sim.modify_body(name="Mars",k2=0.164,Q=99.5)

# Convert to simualtion units
r_p = 25.0e-2 * sim.M2DU  # disk particle size
rho_p = 1500.0 * sim.KG2MU / sim.M2DU**3  # disk particle density
m_p = 4.0 / 3.0 * np.pi * r_p**3 * rho_p  # disk particle mass

sim.add_ring(
    r_p=r_p,
    m_p=m_p,
    mass_distribution={
        "type": "powerlaw",
        "sigma0": 1.05e7 * sim.KG2MU / sim.M2DU**2,
        "alpha": -8.0,
        "nbins": 256,
        "r_outer": 8
    },
)

# Remove ring material in the region beyond Deimos. We still have bins there in case material spreads out that far.
sim.ring["sigma"] = xr.where(sim.ring.r < 6.4*r_cb, sim.ring.sigma, xr.zeros_like(sim.ring.sigma))

tstep_out = 10000.0
dt = 100.0
tstop = 3e8

sim.set_parameter(tstop=tstop, dt=dt, tstep_out=tstep_out, dump_cadence=1)
sim.run()
