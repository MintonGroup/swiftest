#!/usr/bin/env python3
"""
For testing RMVS, the code generates clones of test particles based on one that is fated to impact Mercury.
To use the script, modify the variables just after the  "if __name__ == '__main__':" line
"""
import numpy as np
import swiftest 
from scipy.io import FortranFile
import sys

swifter_input  = "param.swifter.in"
swifter_pl     = "pl.swifter.in"
swifter_tp     = "tp.swifter.in"
swifter_bin    = "bin.swifter.dat"
swifter_enc    = "enc.swifter.dat"

swiftest_input = "param.swiftest.in"
swiftest_pl    = "pl.swiftest.in"
swiftest_tp    = "tp.swiftest.in"
swiftest_cb    = "cb.swiftest.in"
swiftest_bin   = "bin.swiftest.nc"
swiftest_enc   = "enc.swiftest.dat"
swiftest_dis  = "discard.swiftest.dat"

MU2KG = swiftest.MSun
TU2S = swiftest.YR2S
DU2M = swiftest.AU2M

J2 = 0.0 #swiftest.J2Sun
J4 = 0.0 #swiftest.J4Sun

GMSun = swiftest.GMSunSI * TU2S**2 / DU2M**3

# Simple initial conditions of a circular planet with one test particle in a close encounter state
# Simulation start, stop, and output cadence times
t_0	  = 0 # simulation start time
deltaT	= 0.25 * swiftest.JD2S / TU2S   # simulation step size
end_sim = 0.15
t_print = deltaT  #output interval to print results

iout = int(np.ceil(t_print / deltaT))
rmin = swiftest.RSun / swiftest.AU2M
rmax = 1000.0

npl = 1
ntp = 1
plid = 1
tpid = 2

cbname = "Sun"
plname = "Planet"
tpname = "TestParticle"

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

#Make Swifter files
with open(swifter_pl, 'w') as plfile:
   print(npl+1, f'! Planet input file generated using init_cond.py',file=plfile)
   print(0,GMSun,file=plfile)
   print('0.0 0.0 0.0',file=plfile)
   print('0.0 0.0 0.0',file=plfile)
   print(plid,"{:.23g}".format(Gmass),rhill, file=plfile)
   print(radius, file=plfile)
   print(*p_pl, file=plfile)
   print(*v_pl, file=plfile)

with open(swifter_tp, 'w') as tpfile:
   print(1,file=tpfile)
   print(tpid, file=tpfile)
   print(*p_tp, file=tpfile)
   print(*v_tp, file=tpfile)

sys.stdout = open(swifter_input, "w")
print(f'! Swifter input file generated using init_cond.py')
print(f'T0            {t_0} ')
print(f'TSTOP         {end_sim}')
print(f'DT            {deltaT}')
print(f'PL_IN         {swifter_pl}')
print(f'TP_IN         {swifter_tp}')
print(f'IN_TYPE       ASCII')
print(f'ISTEP_OUT     {iout:d}')
print(f'ISTEP_DUMP    {iout:d}')
print(f'BIN_OUT       {swifter_bin}')
print(f'OUT_TYPE      REAL8')
print(f'OUT_FORM      XV')
print(f'OUT_STAT      UNKNOWN')
print(f'J2            {J2}')
print(f'J4            {J4}')
print(f'CHK_CLOSE     yes')
print(f'CHK_RMIN      {rmin}')
print(f'CHK_RMAX      {rmax}')
print(f'CHK_EJECT     {rmax}')
print(f'CHK_QMIN      {rmin}')
print(f'CHK_QMIN_COORD HELIO')
print(f'CHK_QMIN_RANGE {rmin} {rmax}')
print(f'ENC_OUT        {swifter_enc}')
print(f'EXTRA_FORCE    no')
print(f'BIG_DISCARD    no')
print(f'RHILL_PRESENT  yes')
sys.stdout = sys.__stdout__

#Now make Swiftest files
with open(swiftest_cb, 'w') as cbfile:
   print(cbname,file=cbfile)
   print(GMSun, file=cbfile)
   print(rmin, file=cbfile)
   print(J2, file=cbfile)
   print(J4, file=cbfile)

with open(swiftest_pl, 'w') as plfile:
   print(npl, f'! Planet input file generated using init_cond.py',file=plfile)
   print(plname,"{:.23g}".format(Gmass),rhill, file=plfile)
   print(radius, file=plfile)
   print(*p_pl, file=plfile)
   print(*v_pl, file=plfile)

with open(swiftest_tp, 'w') as tpfile:
   print(ntp,file=tpfile)
   print(tpname, file=tpfile)
   print(*p_tp, file=tpfile)
   print(*v_tp, file=tpfile)

sys.stdout = open(swiftest_input, "w")
print(f'! Swiftest input file generated using init_cond.py')
print(f'T0             {t_0} ')
print(f'TSTOP          {end_sim}')
print(f'DT             {deltaT}')
print(f'CB_IN          {swiftest_cb}')
print(f'PL_IN          {swiftest_pl}')
print(f'TP_IN          {swiftest_tp}')
print(f'IN_TYPE        ASCII')
print(f'IN_FORM        XV')
print(f'ISTEP_OUT      {iout:d}')
print(f'ISTEP_DUMP     {10*iout:d}')
print(f'BIN_OUT        {swiftest_bin}')
print(f'OUT_TYPE       NETCDF_DOUBLE')
print(f'OUT_FORM       XVEL')
print(f'OUT_STAT       REPLACE')
print(f'CHK_CLOSE      yes')
print(f'CHK_RMIN       {rmin}')
print(f'CHK_RMAX       {rmax}')
print(f'CHK_EJECT      {rmax}')
print(f'CHK_QMIN       {rmin}')
print(f'CHK_QMIN_COORD HELIO')
print(f'CHK_QMIN_RANGE {rmin} {rmax}')
print(f'ENC_OUT        {swiftest_enc}')
print(f'DISCARD_OUT    {swiftest_dis}')
print(f'EXTRA_FORCE    no')
print(f'BIG_DISCARD    no')
print(f'ROTATION       no')
print(f'GR             no')
print(f'MU2KG          {MU2KG}')
print(f'DU2M           {DU2M}')
print(f'TU2S           {TU2S}')
print(f'RHILL_PRESENT  yes')
print(f'GMTINY          1e-12')




