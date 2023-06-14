#!/usr/bin/env python3
import numpy as np
import swiftest 

# Simple initial conditions of a circular planet with one test particle in a close encounter state
# Simulation start, stop, and output cadence times
t_0	  = 0 # simulation start time
deltaT	= 0.25 * swiftest.JD2S / swiftest.YR2S   # simulation step size
end_sim = 0.15
t_print = deltaT  #output interval to print results

iout = int(np.ceil(t_print / deltaT))

radius = np.double(4.25875607065041e-05)
Gmass = np.double(0.00012002693582795244940133) 
apl = np.longdouble(1.0)
atp = np.longdouble(1.01)
vpl = np.longdouble(2 * np.pi)
vtp = np.longdouble(2 * np.pi / np.sqrt(atp))

p_pl = np.array([apl, 0.0, 0.0], dtype=np.double)
v_pl = np.array([0.0, vpl, 0.0], dtype=np.double)

p_tp = np.array([atp, 0.0, 0.0], dtype=np.double)
v_tp = np.array([0.0, vtp, 0.0], dtype=np.double)

rhill = np.double(apl * 0.0100447248332378922085)

sim = swiftest.Simulation(simdir="swiftest_sim",init_cond_format="XV", general_relativity=False, integrator="RMVS")
sim.clean()
sim.add_solar_system_body(["Sun"])
sim.add_body(name=["Planet"], rh=[p_pl], vh=[v_pl], Gmass=[Gmass], radius=[radius], rhill=[rhill])
sim.add_body(name=["TestParticle"], rh=[p_tp], vh=[v_tp])
sim.set_parameter(tstart=t_0, tstop=end_sim, dt=deltaT, istep_out=iout, dump_cadence=0)
sim.save()

sim.set_parameter(simdir="swifter_sim",codename="Swifter",init_cond_file_type="ASCII",output_format="XV",output_file_type="REAL8")
sim.save()