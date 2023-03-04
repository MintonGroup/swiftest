#!/usr/bin/env python3
import swiftest
sim = swiftest.Simulation(param_file="param.in", codename="Swifter")
sim.param = swiftest.io.swifter2swift_param(sim.param)
sim.write_param(param_file="param.swift.in")
outfile = open("pl.in", 'w')
with open("Jul27_32k_fully.in") as f:
   line = f.readline()
   npl = line.split()[0]
   print(npl, file=outfile) 
   line = f.readline()
   Mcb = line.split()[1]
   print(Mcb, file=outfile) 
   line = f.readline().strip()
   print(line, file=outfile)
   line = f.readline().strip()
   print(line, file=outfile)
   for i in range(int(npl)-1):
      line = f.readline()
      GMpl = line.split()[1]
      Rhill = line.split()[2]
      line = f.readline()
      Rpl = line.split()[0]
      print(GMpl, Rhill, Rpl, file=outfile)
      line = f.readline().strip()
      print(line, file=outfile)
      line = f.readline().strip()
      print(line, file=outfile)
