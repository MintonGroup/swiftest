import numpy as np
from matplotlib import pyplot as plt

import swiftest

cb = swiftest.get_solar_system_body("Saturn")
r_cb = cb["radius"]
m_cb = cb["mass"]
gm_cb = cb["Gmass"]
rho_cb = m_cb / (4.0 / 3.0 * np.pi * r_cb**3)
rot_cb = np.sqrt(cb["rot"][0] ** 2 + cb["rot"][1] ** 2 + cb["rot"][2] ** 2)
sim = swiftest.Simulation(integrator="ringmoons", MU="kg", DU="km", TU="yr")
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
    "mu" : 110.0e3,   # radius of the center of the gaussian
    "dev" : 3600.0,  # width of the gaussian
    "sigma0" : 6.15e4,  # peak of the surface mass density at the center of the gaussian
    "nbins" : 1024
}

sim.add_ring(
    r_p=r_p,
    m_p=m_p,
    mass_distribution=mass_distribution
)

tstep_out = 100.0
dt = 100.0
tstop = 1e6

sim.set_parameter(tstop=tstop, dt=dt, tstep_out=tstep_out, dump_cadence=0)
sim.save()

