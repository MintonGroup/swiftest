import swiftest
"""
 Reads in an example parameter file for each of the codes (Swift, Swifter, and Swiftest) and outputs an equivalent
 parameter file in another code type.
 
 Input files are param.swift[er|est].in and output files are param.[src]2[dest].new, where [src] is the original file type
 and [dest] is the new one. The user may be prompted to answer questions regarding how to convert between the different
 simulation types.
 """
inparam = "param.swift.in"
outparam = "param.swift2swifter.new"
print(f"Reading Swift parameter {inparam} and saving it to {outparam}")
sim = swiftest.Simulation(source=inparam, codename="Swift")
oldparam = sim.convert(outparam, newcodename="Swifter", plname="pl.swift2swifter.in", tpname="tp.swift2swifter.in")