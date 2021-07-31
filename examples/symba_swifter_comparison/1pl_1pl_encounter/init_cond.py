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
swiftest_bin   = "bin.swiftest.dat"
swiftest_enc   = "enc.swiftest.dat"

MU2KG = swiftest.MSun
TU2S = swiftest.YR2S
DU2M = swiftest.AU2M

GMSun = swiftest.GMSunSI * TU2S**2 / DU2M**3

# Simple initial conditions of a circular planet with one smaller massive body in a close encounter state
# Simulation start, stop, and output cadence times
t_0	  = 0 # simulation start time
deltaT	= 0.25 * swiftest.JD2S / TU2S   # simulation step size
end_sim = 0.15
t_print = deltaT  #output interval to print results

iout = int(np.ceil(t_print / deltaT))
rmin = swiftest.RSun / swiftest.AU2M
rmax = 1000.0

npl = 2
ntp = 0
plid1 = 2
plid2 = 100

radius1 = np.double(4.25875607065041e-05)
mass1 = np.double(0.00012002693582795244940133) 
mass2 = mass1 / 100.0
radius2 = radius1 * (mass2 / mass1)**(1.0/3.0)

apl1 = np.longdouble(1.0)
apl2 = np.longdouble(1.01)
vpl1 = np.longdouble(2 * np.pi)
vpl2 = np.longdouble(2 * np.pi / np.sqrt(apl2))

p_pl1 = np.array([apl1, 0.0, 0.0], dtype=np.double)
v_pl1 = np.array([0.0, vpl1, 0.0], dtype=np.double)

p_pl2 = np.array([apl2, 0.0, 0.0], dtype=np.double)
v_pl2 = np.array([0.0, vpl2, 0.0], dtype=np.double)

Rhill1 = np.double(apl1 * 0.0100447248332378922085)
Rhill2 = Rhill1 * (mass2 / mass1)**(1.0 / 3.0)

#Make Swifter files
plfile = open(swifter_pl, 'w')
print(npl+1, f'! Planet input file generated using init_cond.py',file=plfile)

print(1,GMSun,file=plfile)
print('0.0 0.0 0.0',file=plfile)
print('0.0 0.0 0.0',file=plfile)

print(plid1,"{:.23g}".format(mass1),Rhill1, file=plfile)
print(radius1, file=plfile)
print(*p_pl1, file=plfile)
print(*v_pl1, file=plfile)

print(plid2,"{:.23g}".format(mass2),Rhill2, file=plfile)
print(radius2, file=plfile)
print(*p_pl2, file=plfile)
print(*v_pl2, file=plfile)

plfile.close()

tpfile = open(swifter_tp, 'w')
print(0,file=tpfile)
tpfile.close()

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
#print(f'J2            {swiftest.J2Sun}')
#print(f'J4            {swiftest.J4Sun}')
print(f'J2            0.0')
print(f'J4            0.0')
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
cbfile = FortranFile(swiftest_cb, 'w')
Msun = np.double(1.0)
cbfile.write_record(0)
cbfile.write_record(np.double(GMSun))
cbfile.write_record(np.double(rmin))
#cbfile.write_record(np.double(swiftest.J2Sun))
#cbfile.write_record(np.double(swiftest.J4Sun))
cbfile.write_record(np.double(0.0))
cbfile.write_record(np.double(0.0))
cbfile.close()

plfile = FortranFile(swiftest_pl, 'w')
plfile.write_record(npl)

plfile.write_record(np.array([plid1, plid2]))
plfile.write_record(np.vstack([p_pl1[0],p_pl2[0]]))
plfile.write_record(np.vstack([p_pl1[1],p_pl2[1]]))
plfile.write_record(np.vstack([p_pl1[2],p_pl2[2]]))
plfile.write_record(np.vstack([v_pl1[0],v_pl2[0]]))
plfile.write_record(np.vstack([v_pl1[1],v_pl2[1]]))
plfile.write_record(np.vstack([v_pl1[2],v_pl2[2]]))
plfile.write_record(np.array([mass1,mass2]))
plfile.write_record(np.array([Rhill1,Rhill2]))
plfile.write_record(np.array([radius1,radius2]))
plfile.close()
tpfile = FortranFile(swiftest_tp, 'w')
tpfile.write_record(ntp)

tpfile.close()

sys.stdout = open(swiftest_input, "w")
print(f'! Swiftest input file generated using init_cond.py')
print(f'T0             {t_0} ')
print(f'TSTOP          {end_sim}')
print(f'DT             {deltaT}')
print(f'CB_IN          {swiftest_cb}')
print(f'PL_IN          {swiftest_pl}')
print(f'TP_IN          {swiftest_tp}')
print(f'IN_TYPE        REAL8')
print(f'ISTEP_OUT      {iout:d}')
print(f'ISTEP_DUMP     {iout:d}')
print(f'BIN_OUT        {swiftest_bin}')
print(f'OUT_TYPE       REAL8')
print(f'OUT_FORM       XV')
print(f'OUT_STAT       REPLACE')
print(f'CHK_CLOSE      yes')
print(f'CHK_RMIN       {rmin}')
print(f'CHK_RMAX       {rmax}')
print(f'CHK_EJECT      {rmax}')
print(f'CHK_QMIN       {rmin}')
print(f'CHK_QMIN_COORD HELIO')
print(f'CHK_QMIN_RANGE {rmin} {rmax}')
print(f'ENC_OUT        {swiftest_enc}')
print(f'EXTRA_FORCE    no')
print(f'BIG_DISCARD    no')
print(f'ROTATION       no')
print(f'GR             no')
print(f'MU2KG          {MU2KG}')
print(f'DU2M           {DU2M}')
print(f'TU2S           {TU2S}')
print(f'RHILL_PRESENT  yes')
print(f'MTINY          1e-12')




