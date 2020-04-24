import rebound
import numpy as np

sim = rebound.Simulation()
sim.units = ('AU', 'yr', 'Msun')
sim.add(m=1.)
sim.move_to_com()
ps = sim.particles

#!!!!!!! CHANGE THESE THINGS !!!!!!!!!!!!!

sim.convert_particle_units('AU', 'd', 'Msun')

N_tp = 3000
N_ps = 2001

#TP_OUTFILE_SWIFT = open('Feb25_tp_2_swift.txt', 'w')
TP_OUTFILE_SWIFTER = open('tp.in', 'w')

np.random.seed(2)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def uniform(minimum, maximum):
    return np.random.uniform()*(maximum-minimum)+minimum

while sim.N < (1+N_tp):
    a_tp = uniform(1.7,2.0)
    e_tp = uniform(0.0, 0.3)
    inc_tp = uniform(0.0, 0.3)
    O_tp = uniform(0,2*np.pi)
    o_tp = uniform(0,2*np.pi)
    M_tp = uniform(-np.pi, np.pi)
    tp = rebound.Particle(simulation=sim,primary=sim.particles[0],m=0.0, r=0.0, a=a_tp, e=e_tp, inc=inc_tp, Omega=O_tp, omega=o_tp, M=M_tp)
    sim.add(tp)

x_tp = [ps[i].x for i in range(0, sim.N)]
y_tp = [ps[i].y for i in range(0, sim.N)]
z_tp = [ps[i].z for i in range(0, sim.N)]
vx_tp = [ps[i].vx for i in range(0, sim.N)]
vy_tp = [ps[i].vy for i in range(0, sim.N)]
vz_tp = [ps[i].vz for i in range(0, sim.N)]

#with TP_OUTFILE_SWIFT as output:
#    output.write("%s      \n" %(N_tp)) #number of particles in the system
#    for i in range (1, sim.N):
#        output.write("%s %s %s        \n" % ("{:10.8e}".format(x_tp[i]), "{:10.8e}".format(y_tp[i]), "{:10.8e}".format(z_tp[i]))) #x y z
#        output.write("%s %s %s        \n" % ("{:10.8e}".format(vx_tp[i]), "{:10.8e}".format(vy_tp[i]), "{:10.8e}".format(vz_tp[i]))) #vx vy vz
#        output.write("0 0 0 0 0 0 0 0 0 0 0 0 0\n") #flags
#        output.write("0.0d0 0.0d0 0.0d0 0.0d0 0.0d0\n") #flags
#        output.write("0.0d0 0.0d0 0.0d0 0.0d0 0.0d0\n") #flags
#        output.write("0.0d0 0.0d0 0.0d0\n") #flags

with TP_OUTFILE_SWIFTER as output:
    output.write("%s      \n" %(N_tp)) #number of particles in the system
    for i in range (1, sim.N):
        output.write("%s              \n" % ((i+N_ps) )) #ID
        output.write("%s %s %s        \n" % ("{:10.8e}".format(x_tp[i]), "{:10.8e}".format(y_tp[i]), "{:10.8e}".format(z_tp[i]))) #x y z
        output.write("%s %s %s        \n" % ("{:10.8e}".format(vx_tp[i]), "{:10.8e}".format(vy_tp[i]), "{:10.8e}".format(vz_tp[i]))) #vx vy vz