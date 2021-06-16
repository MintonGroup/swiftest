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

MU2KG    = np.longdouble(GMSunSI / GC) #Conversion from mass unit to kg
DU2M     = np.longdouble(AU2M)         #Conversion from radius unit to centimeters
TU2S     = np.longdouble(year)         #Conversion from time unit to seconds
GU       = np.longdouble(GC / (DU2M**3 / (MU2KG * TU2S**2)))
GMSun = np.longdouble(GMSunSI / (DU2M**3 / TU2S**2))

# Solar oblatenes values: From Mecheri et al. (2004), using Corbard (b) 2002 values (Table II)
J2 = 0.0 #np.longdouble(2.198e-7) * (Rsun / DU2M)**2
J4 = 0.0 #np.longdouble(-4.805e-9) * (Rsun / DU2M)**4

npl = 9
ntp = 2 * npl

# Planet ids
plid = {
   'mercury'     : '1',
   'venus'       : '2',
   'earthmoon'   : '3',
   'mars'        : '4',
   'jupiter'     : '5',
   'saturn'      : '6',
   'uranus'      : '7',
   'neptune'     : '8',
   'plutocharon' : '9'
}

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

# Both codes use the same tp input file
tpin = "tp.in"

swifter_input  = "param.swifter.in"
swifter_pl     = "pl.swifter.in"
swifter_bin    = "bin.swifter.dat"
swifter_enc    = "enc.swifter.dat"

swiftest_input = "config.swiftest.in"
swiftest_pl    = "pl.swiftest.in"
swiftest_cb    = "cb.swiftest.in"
swiftest_bin   = "bin.swiftest.dat"
swiftest_enc   = "enc.swiftest.dat"

# Simple initial conditions of a circular planet with one test particle in a close encounter state
# Simulation start, stop, and output cadence times
t_0	  = 0 # simulation start time
deltaT	= 1.00 * JD / TU2S   # simulation step size
end_sim = year / TU2S # simulation end time
t_print = deltaT  #output interval to print results

iout = int(np.ceil(t_print / deltaT))
rmin = Rsun / DU2M
rmax = 1000.0

# Dates to fetch planet ephemerides from JPL Horizons
tstart = '2021-06-15'
tend = '2021-06-16'
tstep = '1d'

#######################################################
# Start generating initial conditions for the planets #
#######################################################
pldata = {}
tpdata = {}
plvec = {}
tpvec = {}
Rhill = {}

for key,val in plid.items():
   pldata[key] = Horizons(id=val, id_type='majorbody',location='@sun',
            epochs={'start': tstart, 'stop': tend,
            'step': tstep})
   plvec[key] = np.array([
         pldata[key].vectors()['x'][0],
         pldata[key].vectors()['y'][0],
         pldata[key].vectors()['z'][0],
         pldata[key].vectors()['vx'][0],
         pldata[key].vectors()['vy'][0],
         pldata[key].vectors()['vz'][0]
         ])
   Rhill[key] = pldata[key].elements()['a'][0] * (3 * MSun_over_Mpl[key])**(-1.0 / 3.0)

# For each planet, we will initialize a pair of test particles. One on its way in, and one on its way out. We will also initialize two additional particles that don't encounter anything
tpnames = range(101, 101 + ntp)
tpxv1 = np.empty((6))
tpxv2 = np.empty((6))
for idx,key in enumerate(plid):
    rstart = 2 * Rpl[key] / DU2M # Start the test particles at a multiple of the planet radius away
    vstart = 1.5 * np.sqrt(2 * GMSun / MSun_over_Mpl[key] / rstart)  # Start the test particle velocities at a multiple of the escape speed
    xvstart = np.array([rstart / np.sqrt(2.0), rstart / np.sqrt(2.0), 0.0, vstart, 0.0, 0.0])
    # The positions and velocities of each pair of test particles will be in reference to a planet
    tpxv1 = plvec[key] + xvstart
    tpxv2 = plvec[key] - xvstart
    tpvec[tpnames[idx]] = tpxv1
    tpvec[tpnames[idx + npl]] = tpxv2

# TP file
tpfile = open(tpin, 'w')
print(ntp, file=tpfile)
for tpid, tp in tpvec.items():
    print(tpid, file=tpfile)
    print(f'{tp[0]} {tp[1]} {tp[2]}', file=tpfile)
    print(f'{tp[3]} {tp[4]} {tp[5]}', file=tpfile)
tpfile.close()

# Swifter PL file
plfile = open(swifter_pl, 'w')
print(npl + 1, file=plfile)
print(1,GMSun,file=plfile)
print('0.0 0.0 0.0',file=plfile)
print('0.0 0.0 0.0',file=plfile)
for key, pl in plvec.items():
    print(f'{int(plid[key])+1} {GMSun / MSun_over_Mpl[key]} {Rhill[key]} ! {key}', file=plfile)
    print(f'{Rpl[key] / DU2M}', file=plfile)
    print(f'{pl[0]} {pl[1]} {pl[2]}', file=plfile)
    print(f'{pl[3]} {pl[4]} {pl[5]}', file=plfile)
plfile.close()

# Swiftest Central body and pl file
cbfile = open(swiftest_cb, 'w')
print(GMSun, file=cbfile)
print(rmin, file=cbfile)
print(J2, file=cbfile)
print(J4, file=cbfile)
cbfile.close()
plfile = open(swiftest_pl, 'w')
print(npl, file=plfile)
for key, pl in plvec.items():
    print(f'{int(plid[key]) + 1} {GMSun / MSun_over_Mpl[key]} ! {key}', file=plfile)
    print(f'{Rpl[key] / DU2M}', file=plfile)
    print(f'{pl[0]} {pl[1]} {pl[2]}', file=plfile)
    print(f'{pl[3]} {pl[4]} {pl[5]}', file=plfile)
plfile.close()

# Swifter parameter file
sys.stdout = open(swifter_input, "w")
print(f'! Swifter input file generated using init_cond.py')
print(f'T0            {t_0} ')
print(f'TSTOP         {end_sim}')
print(f'DT            {deltaT}')
print(f'PL_IN         {swifter_pl}')
print(f'TP_IN         {tpin}')
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

# Swiftest configuration file
sys.stdout = open(swiftest_input, "w")
print(f'! Swiftest input file generated using init_cond.py')
print(f'T0             {t_0} ')
print(f'TSTOP          {end_sim}')
print(f'DT             {deltaT}')
print(f'CB_IN          {swiftest_cb}')
print(f'PL_IN          {swiftest_pl}')
print(f'TP_IN          {tpin}')
print(f'IN_TYPE        ASCII')
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





