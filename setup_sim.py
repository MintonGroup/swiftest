import numpy as np
import swiftest
from astropy import constants as const
from astroquery.jplhorizons import Horizons
import datetime

seed = None
rng = np.random.default_rng(seed=seed)
tstop = 1e3 * 365.25 # 1e3 * 365.25 # years
dt = 90.0 / 60.0 / 24 / 10 # 1/10 * 90 min to days
dt_unit = 'Years'
dt_max = 0

scratch_dr = '/scratch/bell/anand43/swiftest_runs/'
body_dr = 'tests/'
simdir = 'test' # '2e3_part_radial'
simdir = scratch_dr + body_dr + simdir

# unit conversion factors
S2DAY = 1 / (60.0 * 60 * 24)
S2YR = S2DAY / 365.25
KG2MSUN = 1 / 1.98847e30
M2AU = 1 / 1.496e11

# initial central body parameters
density = 2000 # kg/m^3
radius = 1.5e6 # m
volume = 4.0 / 3.0 * np.pi * radius**3
rot = [[0, 0, 0.00058]] # rad/s

#initial ring particle orbital parameters

n_bodies = 1 # number of bodies
id_start = 100 # starting ID for bodies
random_2pi = 360 * np.random.rand(n_bodies) # 360 degrees
random_pi = 180 * np.random.rand(n_bodies) # 180 degrees
# ISS parameters on 4/28/23
# random_a = 6.7999524e3                     # (rng.random(n_bodies) * (3.0 - 0.5) + 0.5) # randomize the semi-major axis between 0.5 and 3.0 (AU))
# random_e = 4.3357182e-4                    # rng.random(n_bodies) # randomize eccentricity between 0 and 1
# random_i = 6.7813805e1                     # random_pi / 2 + (np.pi/2 - np.pi/3) # between (-pi/3, pi/3) + pi/2 (equatorial)
# random_capom = 2.2163223e2                 # random_2pi
# random_omega = 4.4320487e1                 # random_2pi
# random_capm = 3.2698275e2                  # random_2pi
# random_radius = 100e-3 # km                # (rng.random(n_bodies) * (10.0 - 1.0) + 1.0) * M2AU # randomize the radii of the ring particles between 1 and 10 m
# random_mass = 419725 # kg                  # (np.power(random_radius, 3) * 4.0 * np.pi / 3.0 * density) * KG2MSUN # randomize the mass 

# use JPL Horizons for ISS details
# date = '2023-04-27'
ephemerides_start_date = '2023-05-18'
tstart = datetime.date.fromisoformat(ephemerides_start_date)
tstep = datetime.timedelta(days=1)
tend = tstart + tstep
ephemerides_end_date = tend.isoformat()
ephemerides_step = '1d'
iss_obj = Horizons(id = '-125544', location = '399', epochs={'start': ephemerides_start_date, 'stop': ephemerides_end_date,
                                       'step': ephemerides_step}) # 399 = Earth geocenter; -125544 = ISS

iss_el = iss_obj.elements()
iss_el['a'].convert_unit_to('km')
random_a = iss_el['a'][0]
random_e = iss_el['e'][0]
random_i = iss_el['incl'][0]
random_capom = iss_el['Omega'][0]
random_omega = iss_el['w'][0]
random_capm = iss_el['M'][0]
random_radius = 100e-3 # km
random_mass = 419725 # km

earth_radius = 6371.01 # km
earth_mass = 5.97219e24 # kg
earth_J2 = 1083e-6 # from Murray and Dermott txt. bk.
earth_J4 = -2e-6 # from Murray and Dermott txt. bk.
earth_tilt = np.deg2rad(23.4392811)
earth_rot = np.array([0.0, np.sin(earth_tilt), np.cos(earth_tilt)]) * 2 * np.pi / (24 * 60 * 60.0)
earth_rot = [earth_rot]

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

sim = swiftest.Simulation(simdir = simdir, integrator = "symba", tstop = tstop, dt = dt, istep_out = 1e7, dump_cadence = 10, rotation = True, collision_model = "FRAGGLE", MU = 'kg', DU2M = 1e3, TU = 'd', rmin = earth_radius)
# sim.add_solar_system_body(["Earth"], date = ephemerides_start_date)
# Manually add Earth as a CB because no J2 or J4 term is added otherwise
sim.add_body(name = "Earth", id = 1, a = 0, e = 0, inc = 0, capom = 0, omega = 0, capm = 0, mass = earth_mass, radius = earth_radius, J2 = earth_J2 * earth_radius**2, J4 = earth_J4 * earth_radius**4, rot = earth_rot)

# sim.add_body(name = 'Centaur', id = 1, mass = density * volume, radius = radius, a = 0, e = 0, inc = 0, capom = 0, omega = 0, capm = 0)
sim.add_body(name = "ISS", id = np.arange(id_start, n_bodies + id_start, step = 1), a = random_a, e = random_e, inc = random_i, capom = random_capom, omega = random_omega, capm = random_capm, mass = random_mass, radius = random_radius)
# sim.get_parameter()
# sim.set_parameter()
sim.write_param()
sim.save()
# sim.run()

print(f'total number of steps ={int(tstop / dt)}')
J2 = sim.init_cond.isel(name = 0)['j2rp2'].values
print(f'J2 in init_cond = {J2}')
d = sim.data.isel(name = 0)
J2 = d['j2rp2'].values
print(f'J2 for central body = {J2}')
