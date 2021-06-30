import swiftest
import swiftest.io as swio
param_file = "param.swiftest.in"
sim = swiftest.Simulation(source="param.swiftest.in")
ds = swio.swiftest2xr(sim.param)