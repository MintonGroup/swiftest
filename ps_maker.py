import rebound
import numpy as np

def setupSimulation():
    sim = rebound.Simulation()
    sim.units = ('AU', 'yr', 'Msun')
    sim.integrator = "mercurius"
    sim.dt = 0.008
    sim.testparticle_type = 1
    sim.move_to_com()
    return sim

sim = setupSimulation()
ps = sim.particles
G_auy = 4 * np.pi * np.pi #G in units of AU^3 year^-2 M_sun^-1
M_Sun = 1
M_Sun_to_g = 1.989e33
M_Sun_to_kg = 1.989e30
AU_cubed_to_cm_cubed = 3.348e39
AU_cubed_to_km_cubed = 3.348e24
year_to_seconds = 3.154e7
Mtot_disk = 6.006e-6 #~3*M_earth

OUTFILE = open('pl.in', 'w')

sim.add( m=1.0, hash="sun")  # SUN - Adds a particle of mass 1

N_fully = 2001
N_semi = 0

d_bodies = 2.0 * AU_cubed_to_cm_cubed * (1/M_Sun_to_g) #Changes 2 g/cm^3 to 3366515.837 M_sun/AU^3

m_semi =  Mtot_disk / (N_semi + 2*N_fully)
m_fully = 2 * m_semi

r_semi = ((3*m_semi)/(4*np.pi*d_bodies))**(1/3)
r_fully = ((3*m_fully)/(4*np.pi*d_bodies))**(1/3)

np.random.seed(1)

def uniform(minimum, maximum):
    return np.random.uniform()*(maximum-minimum)+minimum

while sim.N < N_fully:
    a_fully = uniform(0.5,1)
    e_fully = uniform(0.0, 0.3)
    inc_fully = uniform(0.0, 0.3)
    O_fully = uniform(0,2*np.pi)
    o_fully = uniform(0,2*np.pi)
    M_fully = uniform(-np.pi, np.pi)
    fully = rebound.Particle(simulation=sim,primary=sim.particles[0],m=m_fully, r=r_fully, a=a_fully, e=e_fully, inc=inc_fully, Omega=O_fully, omega=o_fully, M=M_fully)
    sim.add(fully)

while sim.N < (N_fully+N_semi):
    a_semi = uniform(0.5,1)
    e_semi = uniform(0.0, 0.3)
    inc_semi = uniform(0.0, 0.3)
    O_semi = uniform(0,2*np.pi)
    o_semi = uniform(0,2*np.pi)
    M_semi = uniform(-np.pi, np.pi)
    semi = rebound.Particle(simulation=sim,primary=sim.particles[0],m=m_semi, r=r_semi, a=a_semi, e=e_semi, inc=inc_semi, Omega=O_semi, omega=o_semi, M=M_semi)
    sim.add(semi)

x = [ps[i].x for i in range(1, sim.N)]
y = [ps[i].y for i in range(1, sim.N)]
z = [ps[i].z for i in range(1, sim.N)]
vx = [ps[i].vx for i in range(1, sim.N)]
vy = [ps[i].vy for i in range(1, sim.N)]
vz = [ps[i].vz for i in range(1, sim.N)]
m = [ps[i].m for i in range(1, sim.N)]
r = [ps[i].r for i in range(1, sim.N)]
Rhill = [ps[i].a*((ps[i].m/(3*M_Sun))**(0.333333)) for i in range(1, sim.N)]

with OUTFILE as output:
    output.write("%s      ! Solar System in unit system AU, M_sun, and years\n" %(sim.N))
    output.write("1 %s\n"%"{:10.8e}".format(M_Sun*G_auy))
    output.write(".0 .0 .0        ! x y z\n")
    output.write(".0 .0 .0        !vx vy vz\n")
    for i in range (0, (sim.N-1)):
        output.write("%s %s %s       ! ID / G*Mass / Rhill\n"%((i+2),"{:10.8e}".format(m[i]*G_auy),"{:10.8e}".format(Rhill[i])))
        output.write("%s             ! Radius\n"%("{:10.8e}".format(r[i])))
        output.write("%s %s %s       ! x y z\n"%("{:10.8e}".format(x[i]),"{:10.8e}".format(y[i]),"{:10.8e}".format(z[i])))
        output.write("%s %s %s       ! vx vy vz\n"%("{:10.8e}".format(vx[i]),"{:10.8e}".format(vy[i]),"{:10.8e}".format(vz[i])))






