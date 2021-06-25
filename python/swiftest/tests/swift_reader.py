import swiftest
inparam = "param.swift.in"
outparam = "param.swift.new"
print(f"Reading Swift parameter {inparam} and saving it to {outparam}")
sim = swiftest.Simulation(inparam, codename="Swift")
sim.write_param(outparam)