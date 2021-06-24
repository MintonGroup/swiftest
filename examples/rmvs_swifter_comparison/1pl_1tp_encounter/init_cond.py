"""
For testing RMVS, the code generates clones of test particles based on one that is fated to impact Mercury.
To use the script, modify the variables just after the  "if __name__ == '__main__':" line
"""
import numpy as np
from astroquery.jplhorizons import Horizons
import astropy.constants as const 
import swiftestio as swio
from scipy.io import FortranFile
import sys

#Values from JPL Horizons
AU2M = np.longdouble(const.au.value)
GMSunSI = np.longdouble(const.GM_sun.value)
Rsun = np.longdouble(const.R_sun.value)
GC = np.longdouble(const.G.value)
JD = 86400
year = np.longdouble(365.25 * JD)
c = np.longdouble(299792458.0)
MSun_over_Mpl = np.array([6023600.0,
                 408523.71,
                 328900.56,
                 3098708.,
                 1047.3486,
                 3497.898,
                 22902.98,
                 19412.24,
                 1.35e8], dtype=np.longdouble)

MU2KG    = np.longdouble(GMSunSI / GC) #Conversion from mass unit to kg
DU2M     = np.longdouble(AU2M)         #Conversion from radius unit to centimeters
TU2S     = np.longdouble(year)         #Conversion from time unit to seconds
GU       = np.longdouble(GC / (DU2M**3 / (MU2KG * TU2S**2)))
GMSun = np.longdouble(GMSunSI / (DU2M**3 / TU2S**2))

# Solar oblatenes values: From Mecheri et al. (2004), using Corbard (b) 2002 values (Table II)
J2 = np.longdouble(2.198e-7) * (Rsun / DU2M)**2
J4 = np.longdouble(-4.805e-9) * (Rsun / DU2M)**4

#Planet Msun/M ratio
MSun_over_Mpl = {
   'mercury'     : np.longdouble(6023600.0),
   'venus'       : np.longdouble(408523.71),
   'earthmoon'   : np.longdouble(328900.56),
   'mars'        : np.longdouble(3098708.),
   'jupiter'     : np.longdouble(1047.3486),
   'saturn'      : np.longdouble(3497.898),
   'uranus'      : np.longdouble(22902.98),
   'neptune'     : np.longdouble(19412.24),
   'plutocharon' : np.longdouble(1.35e8)
}

#Planet radii in meters
Rpl = {
   'mercury'     : np.longdouble(2439.4e3),
   'venus'       : np.longdouble(6051.8e3),
   'earthmoon'   : np.longdouble(6371.0084e3), # Earth only for radius
   'mars'        : np.longdouble(3389.50e3),
   'jupiter'     : np.longdouble(69911e3),
   'saturn'      : np.longdouble(58232.0e3),
   'uranus'      : np.longdouble(25362.e3),
   'neptune'     : np.longdouble(24622.e3),
   'plutocharon' : np.longdouble(1188.3e3)
}

THIRDLONG = np.longdouble(1.0) / np.longdouble(3.0)

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

# Simple initial conditions of a circular planet with one test particle in a close encounter state
# Simulation start, stop, and output cadence times
t_0	  = 0 # simulation start time
deltaT	= 0.25 * JD / TU2S   # simulation step size
end_sim = year / TU2S #10 * JD / TU2S # simulation end time
t_print = deltaT  #output interval to print results

iout = int(np.ceil(t_print / deltaT))
rmin = Rsun / DU2M
rmax = 1000.0

npl = 1
plid = 2
tpid = 100

radius = np.double(Rpl['earthmoon'] / DU2M)
mass = np.double(GMSun * MSun_over_Mpl['earthmoon']**-1)
apl = np.longdouble(1.0)
atp = np.longdouble(1.01)
vpl = np.longdouble(2 * np.pi)
vtp = np.longdouble(2 * np.pi / np.sqrt(atp))

p_pl = np.array([apl, 0.0, 0.0], dtype=np.double)
v_pl = np.array([0.0, vpl, 0.0], dtype=np.double)

p_tp = np.array([atp, 0.0, 0.0], dtype=np.double)
v_tp = np.array([0.0, vtp, 0.0], dtype=np.double)

Rhill = apl * ((3 * MSun_over_Mpl['earthmoon'])**(-THIRDLONG))

#Make Swifter files
plfile = open(swifter_pl, 'w')
print(npl+1, f'! Planet input file generated using init_cond.py',file=plfile)
print(1,GMSun,file=plfile)
print('0.0 0.0 0.0',file=plfile)
print('0.0 0.0 0.0',file=plfile)
print(plid,"{:.23g}".format(mass),Rhill, file=plfile)
print(radius, file=plfile)
print(*p_pl, file=plfile)
print(*v_pl, file=plfile)
plfile.close()

tpfile = open(swifter_tp, 'w')
print(1,file=tpfile)
print(tpid, file=tpfile)
print(*p_tp, file=tpfile)
print(*v_tp, file=tpfile)
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
print(f'OUT_STAT      NEW')
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
cbfile = FortranFile(swiftest_cb, 'w')
Msun = np.double(1.0)
cbfile.write_record(np.double(GMSun))
cbfile.write_record(np.double(rmin))
cbfile.write_record(np.double(J2))
cbfile.write_record(np.double(J4))
cbfile.close()

plfile = FortranFile(swiftest_pl, 'w')
plfile.write_record(npl)

plfile.write_record(plid)
plfile.write_record(p_pl[0])
plfile.write_record(p_pl[1])
plfile.write_record(p_pl[2])
plfile.write_record(v_pl[0])
plfile.write_record(v_pl[1])
plfile.write_record(v_pl[2])
plfile.write_record(mass)
plfile.write_record(radius)
plfile.close()
tpfile = FortranFile(swiftest_tp, 'w')
ntp = 1
tpfile.write_record(ntp)
tpfile.write_record(tpid)
tpfile.write_record(p_tp[0])
tpfile.write_record(p_tp[1])
tpfile.write_record(p_tp[2])
tpfile.write_record(v_tp[0])
tpfile.write_record(v_tp[1])
tpfile.write_record(v_tp[2])

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





