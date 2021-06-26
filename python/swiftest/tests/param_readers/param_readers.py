import swiftest
"""
 Reads in an example parameter file for each of the codes (Swift, Swifter, and Swiftest) and outputs a copy using the
 internal parsing system. Input files are param.swift[er|est].in and output files are param.swift[er|est].new.
 The contents of the parameter files should be the same, though the new version may be formatted differently and also
 contain explicit values for default parameters.
 """
inparam = "param.swift.in"
outparam = "param.swift.new"
print(f"Reading Swift parameter {inparam} and saving it to {outparam}")
sim = swiftest.Simulation(source=inparam, codename="Swift")
sim.write_param(outparam)

inparam = "param.swifter.in"
outparam = "param.swifter.new"
print(f"Reading Swifter parameter {inparam} and saving it to {outparam}")
sim = swiftest.Simulation(source=inparam, codename="Swifter")
sim.write_param(outparam)

inparam = "param.swiftest.in"
outparam = "param.swiftest.new"
print(f"Reading Swifter parameter {inparam} and saving it to {outparam}")
sim = swiftest.Simulation(source=inparam) # The default value of codename is "Swiftest"
sim.write_param(outparam)