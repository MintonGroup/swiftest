import swiftest
"""
 Reads in an example parameter file for each of the codes (Swift, Swifter, and Swiftest) and outputs an equivalent
 parameter file in another code type.
 
 Input files are param.swift[er|est].in and output files are param.[src]2[dest].new, where [src] is the original file type
 and [dest] is the new one. The user may be prompted to answer questions regarding how to convert between the different
 simulation types.
 """
inparam = "param.swifter.in"
outparam = "param.swifter2swiftest.new"
print(f"Reading Swift parameter {inparam} and saving it to {outparam}")
sim = swiftest.Simulation(source=inparam, codename="Swifter")
oldparam = sim.convert(outparam, newcodename="Swiftest", plname="pl.swifter2swiftest.in", tpname="tp.swifter2swiftest.in", cbname="cb.swifter2swiftest.in")