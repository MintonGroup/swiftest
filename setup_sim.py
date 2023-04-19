import numpy as np
import swiftest

seed = None
rng = np.random.default_rng(seed=seed)
tstop = 2e6 # years
dt = 8 / 365.25 # 8 DAYS
dt_unit = 'Years'
dt_max = 0
G = 6.6743e-11 # SI Units

# unit conversion factors
S2DAY = 1 / (60.0 * 60 * 24)
S2YR = S2DAY / 365.25
KG2MSUN = 1 / 1.98847e30
M2AU = 1 / 1.496e11

# initial central body parameters
density = 2000 # kg/m^3
radius = 1.5e6 # m
volume = 4.0 / 3.0 * np.pi * radius**3
rot = [0, 0, 0.00058] # rad/s

#initial ring particle orbital parameters

n_bodies = 500 # number of bodies
id_start = 100 # starting ID for bodies
random_2pi = 360 * np.random.rand(n_bodies) # 360 degrees
random_pi = 180 * np.random.rand(n_bodies) # 180 degrees
random_a = (rng.random(n_bodies) * (3.0 - 0.5) + 0.5) # randomize the semi-major axis between 0.5 and 3.0 (AU))
random_e = rng.random(n_bodies) # randomize eccentricity between 0 and 1
random_i = random_pi / 2 + (np.pi/2 - np.pi/3) # between (-pi/3, pi/3) + pi/2 (equatorial)
random_capom = random_2pi
random_omega = random_2pi
random_capm = random_2pi
random_radius = (rng.random(n_bodies) * (10.0 - 1.0) + 1.0) * M2AU # randomize the radii of the ring particles between 1 and 10 m
random_mass = (np.power(random_radius, 3) * 4.0 * np.pi / 3.0 * density) * KG2MSUN # randomize the mass 

# # check that dt < dt_max

# dt_max = 2.0 * np.pi / np.sqrt(G * density * volume) 
# dt_max = dt_max * np.power(np.min(random_a), 3/2) / 10.0
# dt_max = dt_max * S2DYR # convert s to days

# print(f'\ndt_max = {dt_max} {dt_unit}')
# print(f'Current dt = {dt} {dt_unit}\n')

# if(dt > dt_max):
#     print("dt is TOO BIG!")
#     raise Exception

# set up swiftest simulation

sim = swiftest.Simulation(simdir = "kaustub_test", integrator = "symba", tstop = tstop, dt = dt, istep_out = 500, dump_cadence = 1000, rotation = False, collision_model = "FRAGGLE")#, MU = 'kg', DU = 'm', TU = 'd', rmin = radius + 1e4)
sim.add_solar_system_body(["Sun", "Mercury", "Venus", "Earth", "Mars"])
# sim.add_body(name = 'Centaur', id = 1, mass = density * volume, radius = radius, a = 0, e = 0, inc = 0, capom = 0, omega = 0, capm = 0)
sim.add_body(name = np.arange(id_start, n_bodies + id_start, step = 1), id = np.arange(id_start, n_bodies + id_start, step = 1), a = random_a, e = random_e, inc = random_i, capom = random_capom, omega = random_omega, capm = random_capm, mass = random_mass, radius = random_radius)
# sim.get_parameter()
# sim.set_parameter()
sim.write_param()
sim.save()
# sim.run()
