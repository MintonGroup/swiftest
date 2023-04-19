import numpy as np
import swiftest
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

seed = None
rng = np.random.default_rng(seed=seed)
n_bodies = int(1e4) # number of planet bodies
array_shift = np.random.randint(0, n_bodies)
tstop = 500 # rotation periods
dt = 1.0e-3
dt_unit = 'Haumea_rot_period'
# dt_unit = 'd'
dt_max = 0
G = 6.6743e-11 # SI Units
scratch_dr = '/scratch/bell/anand43/swiftest_runs/'
simdir = 'haumea_1e4_part'
simdir = scratch_dr + simdir

# unit conversion factors
S2DAY = 1 / (60.0 * 60 * 24)
S2YR = S2DAY / 365.25
KG2MSUN = 1 / 1.98847e30
M2AU = 1 / 1.496e11

# initial central body parameters
# Haumea parameters taken from P. S. Jean Carvalho (2019)

cb_mass = 4.006e21 # kg
cb_radius = 711.9097e3 # m
mass = cb_mass / cb_mass
radius = cb_radius / cb_radius
longest_axis = 1161e3 # m
rmin = longest_axis / cb_radius
volume = 4.0 / 3.0 * np.pi * radius**3
density = mass / volume # kg/m^3
T_rotation = 3.915341 * 60 * 60 # hours to sec
rot = 2 * np.pi / T_rotation # rad/s
rot = [[0, 0, rot]] 
J2 = 0.305
J4 = -0.216
time_unit = T_rotation
print(f'R_min = {rmin}')

# set up swiftest simulation

sim = swiftest.Simulation(simdir = simdir, integrator = "symba", tstop = tstop, dt = dt, istep_out = 1e5, dump_cadence = 10, rotation = True, collision_model = "FRAGGLE", rmin = rmin, MU2KG = mass, DU2M = radius, TU2S = T_rotation, init_cond_format = 'XV') # TU = dt_unit, MU_name = 'M_Haumea', DU_name = 'R_Haumea'

#initial ring particle orbital parameters

id_start = 100 # starting ID for bodies

random_radius = rng.random(n_bodies) * (50 - 5) + 5 # randomize the radii of the ring particles between 5 and 50 (m)
random_mass = np.power(random_radius, 3) * 4.0 * np.pi / 3.0 * density # randomize the mass 

z_sign = rng.random(n_bodies) - 0.5 # randomize the sign of the z-component

# set up the position vectors
random_pos_mag = (rng.random(n_bodies) * (3 - 1.65) + 1.65) # randomize the distance between 1.65 and 3 Haumea Radii
random_phi = np.deg2rad((rng.random(n_bodies)) * 90.0)
random_theta = np.deg2rad((rng.random(n_bodies)) * (30.0 - (-30)) - 30) # equatorial zone
x = random_pos_mag * np.sin(random_theta) * np.cos(random_phi)
y = random_pos_mag * np.sin(random_theta) * np.sin(random_phi)
z = random_pos_mag * np.cos(random_theta) * np.sign(z_sign)
random_pos_vec = np.array([x, y, z]).T

# set up the velocity vectors
random_vel_mag = np.sqrt(mass * sim.GU / random_pos_mag) * (rng.random(n_bodies) * (2.0 - 0.5) + 0.5) # randomize the velocity by 0.5 - 2.0 about the keplerian velocity
random_phi = np.deg2rad((rng.random(n_bodies)) * 90.0)
random_theta = np.deg2rad((rng.random(n_bodies)) * (30.0 - (-30)) - 30) # equatorial zone
x = random_vel_mag * np.sin(random_theta) * np.cos(random_phi)
y = random_vel_mag * np.sin(random_theta) * np.sin(random_phi)
z = random_vel_mag * np.cos(random_theta) * np.sign(z_sign)
random_vel_vec = np.array([x, y, z]).T

print(f'Shape of pos vec = {random_pos_vec.shape}')
print(f'Shape of vel vec = {random_vel_vec.shape}')

# random_rot = 0 * [0, 0, 1.0] * rng.random(n_bodies) # rad/s
# CHECK and FIX random_rot (sequence multiplication issue)

random_radius = random_radius / radius
random_mass = random_mass / mass

ring_mass_tot = np.sum(random_mass)

print(f'Total ring mass = {ring_mass_tot} {sim.MU_name}')
print(f'sim.GU = {sim.GU}')

sim.add_body(name = 'Centaur', id = 0, mass = mass, rot = rot, radius = radius, rh=[np.array([0,0,0])], vh = [np.array([0,0,0])], J2 = J2, J4 = J4)
sim.add_body(name = np.arange(id_start, n_bodies + id_start, step = 1), id = np.arange(id_start, n_bodies + id_start, step = 1), mass = random_mass, radius = random_radius, rh = random_pos_vec, vh = random_vel_vec)#, rot = random_rot)

# check that dt < dt_max

dt_max = 2.0 * np.pi / np.sqrt(G * mass) 
dt_max = dt_max * np.power(np.min(random_pos_mag) * sim.DU2M, 3/2)
dt_max = dt_max / 10.0
dt_max = dt_max / sim.TU2S # convert s to time units

print(f'\ndt_max = {dt_max} {dt_unit}')
print(f'Current dt = {dt} {dt_unit}\n')

if(dt > dt_max):
    print("dt is TOO BIG!")
    raise Exception

sim.write_param()
sim.save()
# sim.run()

# Make initial plots

fig, ax = plt.subplots(figsize=(8,4.5), dpi=300)

plt.scatter(np.linalg.norm(sim.data.isel(time = 0)['rh'], axis = 1), np.linalg.norm(sim.data.isel(time = 0)['vh'], axis = 1))
plt.xlabel('rh')
plt.ylabel('vh')
# plot longest axis
cb_ax_long = rmin
ymin = np.min(random_vel_mag)
ymax = np.max(random_vel_mag)
plt.ylim([0, ymax])
plt.plot([cb_ax_long, cb_ax_long], [0, ymax], '-', color = 'orange')
plt.text(cb_ax_long - 0.1, ymax / 5, s = r'Longest CB Axis', color = 'orange', rotation = 'vertical')

plt.savefig(simdir + '/r_vs_v_initial.png')

plt.clf()
fig = plt.figure(dpi = 300)
ax = fig.add_subplot(projection = '3d')
ax.scatter(random_pos_vec[:, 0], random_pos_vec[:, 1], random_pos_vec[:, 2], marker = 'o')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_zlim([-5, 5])
plt.savefig(simdir + '/pos_vector.png')

ax.view_init(elev = 0, azim = 0)
plt.savefig(simdir + '/pos_vector_eq.png')


# Extra code

# random_2pi = 360 * np.random.rand(n_bodies) # 360 degrees
# random_pi = 180 * np.random.rand(n_bodies) # 180 degrees
# random_a = (rng.random(n_bodies) * (4 - 1.65) + 1.65)  # randomize the semi-major axis between 1.65 and 4 (Haumea Radii)
#                                                        # Haumea's longest axis = 1.63 * R_Haumea (1161 km = 1.63 * 711.9097 km)
# random_e = rng.random(n_bodies) * 0.8 # randomize eccentricity between 0 and 0.8
# random_i = np.random.rand(n_bodies) * (30 - (-30)) + (-30)  # between (-30, 30) (degrees) (equatorial zone)
# random_capom = random_2pi
# random_omega = np.roll(random_2pi, array_shift)
# array_shift = np.random.randint(0, n_bodies) # compute another random shift
# random_capm = np.roll(random_2pi, array_shift)



# fig, ax = plt.subplots(figsize=(8,4.5), dpi=300)

# plt.scatter(sim.data['a'], sim.data['e'])
# plt.xlabel('a')
# plt.ylabel('e')
# plt.savefig(simdir + '/a_vs_e_initial.png')

# fig, ax = plt.subplots(figsize=(8,4.5), dpi=300)

# plt.scatter(sim.data['a'], sim.data['inc'])
# plt.xlabel('a')
# plt.ylabel('inc')
# plt.savefig(simdir + '/a_vs_i_initial.png')

# print(f'sim.MU2KG = {sim.MU2KG}')
# print(f'sim.param["MU2KG"] = {sim.param["MU2KG"]}')
# print(f'sim.MU_name = {sim.MU_name}')
# sim.add_solar_system_body(["Sun", "Mercury", "Venus", "Earth", "Mars"])
# random_radius = random_radius / sim.DU2M
# random_mass = random_mass / sim.MU2KG
# sim.add_body(name = 'Centaur', id = 0, mass = mass / mass, rot = rot, radius = radius / radius, a = 0, e = 0, inc = 0, capom = 0, omega = 0, capm = 0, J2 = J2, J4 = J4)
# sim.add_body(name = np.arange(id_start, n_bodies + id_start, step = 1), id = np.arange(id_start, n_bodies + id_start, step = 1), a = random_a, e = random_e, inc = random_i, capom = random_capom, omega = random_omega, capm = random_capm, mass = random_mass, radius = random_radius)#, rot = random_rot)

# sim = swiftest.Simulation(simdir = simdir, integrator = "symba", tstop = tstop, dt = dt, istep_out = 1e5, dump_cadence = 10, rotation = True, collision_model = "FRAGGLE", TU = dt_unit, rmin = rmin, MU2KG = mass, DU2M = radius, init_cond_format = 'XV',)# MU_name = 'M_Haumea', DU_name = 'R_Haumea')
