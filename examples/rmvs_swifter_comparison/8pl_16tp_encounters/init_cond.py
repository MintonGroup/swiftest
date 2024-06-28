#!/usr/bin/env python3
import numpy as np
import swiftest

tstart = 0.0
dt = 1.0
tstop = 365.25e2
tstep_out = 100*dt

sim = swiftest.Simulation(simdir="swiftest_sim",init_cond_format="XV",output_format="XV",general_relativity=False, integrator="RMVS", rhill_present=True, MU="Msun", DU="AU", TU="d")
sim.clean()
sim.add_solar_system_body(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"])

plname = sim.init_cond['name'].where(sim.init_cond['name'] != "Sun", drop=True)
pl = sim.init_cond.sel(name=plname)

for i,n in enumerate(pl.name):
    pli = pl.sel(name=n)
    rstart = 2 * pli['radius'].data[0]  # Start the test particles at a multiple of the planet radius away
    vstart = 1.5 * np.sqrt(2 * pli['Gmass'].data[0]  / rstart)  # Start the test particle velocities at a multiple of the escape speed
    rstart_vec = np.array([rstart / np.sqrt(2.0), rstart / np.sqrt(2.0), 0.0])
    vstart_vec = np.array([vstart, 0.0, 0.0])
    rp = pli['rh'].data[0]
    vp = pli['vh'].data[0]
    sim.add_body(name=[f"TestParticle{100+i}",f"TestParticle{200+i}"],rh=[rp+rstart_vec, rp-rstart_vec],vh=[vp+vstart_vec, vp-vstart_vec])


sim.set_parameter(tstart=tstart, tstop=tstop, dt=dt, tstep_out=tstep_out, dump_cadence=0)
sim.save()

sim.set_parameter(simdir="swifter_sim",codename="Swifter",init_cond_file_type="ASCII",output_file_type="REAL8")
sim.save()
