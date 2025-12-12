import numpy as np

import swiftest

mars = swiftest.get_solar_system_body("Mars")
sim = swiftest.Simulation(integrator="ringmoons", MU2KG=mars["mass"], DU2M=mars["radius"], TU="day")
sim.add_solar_system_body(["Mars", "Phobos", "Deimos"], align_to_central_body_rotation=True)

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
    },
)
sim.save()

print(sim.init_cond)
