"""
This script generates initial conditions for the solar using JPL Horizons data.

For testing RMVS, the code generates clones of test particles based on one that is fated to impact Mercury.
To use the script, modify the variables just after the  "if __name__ == '__main__':" line
"""
import numpy as np
import sys
from astroquery.jplhorizons import Horizons
import astropy.constants as const 
import swiftestio as swio
from scipy.io import FortranFile

from numpy.random import default_rng

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

# Simulation start, stop, and output cadence times
t_0	  = 0 # simulation start time
deltaT	= 0.25 * JD / TU2S   # simulation step size
end_sim = 1000 * year  / TU2S # simulation end time
t_print = year / TU2S #output interval to print results

# Solar oblatenes values: From Mecheri et al. (2004), using Corbard (b) 2002 values (Table II)
J2 = np.longdouble(2.198e-7) * (Rsun / DU2M)**2
J4 = np.longdouble(-4.805e-9) * (Rsun / DU2M)**4

tstart = '2021-01-28'
tend = '2021-01-29'
tstep = '1d'
planetid = {
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
npl = 9

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

pdata = {}
plvec = {}
Rhill = {}
THIRDLONG = np.longdouble(1.0) / np.longdouble(3.0)

for key,val in planetid.items():
   pdata[key] = Horizons(id=val, id_type='majorbody',location='@sun',
            epochs={'start': tstart, 'stop': tend,
            'step': tstep})
   plvec[key] = np.array([pdata[key].vectors()['x'][0],
         pdata[key].vectors()['y'][0], 
         pdata[key].vectors()['z'][0], 
         pdata[key].vectors()['vx'][0], 
         pdata[key].vectors()['vy'][0], 
         pdata[key].vectors()['vz'][0] 
         ])
   Rhill[key] = np.longdouble(pdata[key].elements()['a'][0]) * (3 * MSun_over_Mpl[key])**(-THIRDLONG)

if __name__ == '__main__':

   nclones = 10
   xv_dispersion_factor = 0.01 # randomly alter x and v vectors by this fraction of the original

   # Get tp data from collidor simulation to produce the parent of the clones
   inparfile  = "param.tpcollider.in"
   paramfile = swio.read_swifter_param(inparfile)
   swifterdat = swio.swifter2xr(paramfile)
   px = swifterdat.isel(time=-1).sel(id=100)['px'].values.item()
   py = swifterdat.isel(time=-1).sel(id=100)['py'].values.item()
   pz = swifterdat.isel(time=-1).sel(id=100)['pz'].values.item()

   vx = swifterdat.isel(time=-1).sel(id=100)['vx'].values.item()
   vy = swifterdat.isel(time=-1).sel(id=100)['vy'].values.item()
   vz = swifterdat.isel(time=-1).sel(id=100)['vz'].values.item()

   jangofett = np.array([px, py, pz, -vx, -vy, -vz])

   # generate random clones
   rng = default_rng()
   clone_xv = (rng.standard_normal((nclones, 6)) - 0.5) * xv_dispersion_factor + 1.0
   clonetroops = jangofett * clone_xv
   clonenames = range(100, 100 + nclones)
   tpvec = dict(zip(clonenames,clonetroops))

   # Convert from AU-day to AU-year just because I find it easier to keep track of the sim progress
   for plid in plvec:
      plvec[plid][3:] *= year / JD

   for tpid in tpvec:
      tpvec[tpid][3:] *= year / JD

   # Names of all output files
   swifter_input  = "param.swifter.in"
   swifter_pl     = "pl.swifter.in"
   swifter_tp     = "tp.swifter.in"
   swifter_bin    = "bin.swifter.dat"
   swifter_enc    = "enc.swifter.dat"

   swiftest_input = "config.swiftest.in"
   swiftest_pl    = "pl.swiftest.in"
   swiftest_tp    = "tp.swiftest.in"
   swiftest_cb    = "cb.swiftest.in"
   swiftest_bin   = "bin.swiftest.dat"
   swiftest_enc   = "enc.swiftest.dat"

   iout = int(np.ceil(t_print / deltaT))
   rmin = Rsun / DU2M
   rmax = 1000.0

   #Make Swifter files
   plfile = open(swifter_pl, 'w')
   print(npl+1, f'! Planet input file generated using init_cond.py using JPL Horizons data for the major planets (and Pluto) for epoch {tstart}' ,file=plfile)
   print(1,GMSun,file=plfile)
   print('0.0 0.0 0.0',file=plfile)
   print('0.0 0.0 0.0',file=plfile)
   for i, plid in enumerate(plvec): 
      print(i + 2,"{:.23g}".format(GMSun * MSun_over_Mpl[plid]**-1),Rhill[plid], file=plfile)
      print(Rpl[plid] / DU2M, file=plfile)
      print(plvec[plid][0],plvec[plid][1],plvec[plid][2], file=plfile)
      print(plvec[plid][3],plvec[plid][4],plvec[plid][5], file=plfile)
   plfile.close()

   tpfile = open(swifter_tp, 'w')
   print(nclones,file=tpfile)
   for tpid, tp in tpvec.items():
      print(tpid, file=tpfile)
      print(f'{tp[0]} {tp[1]} {tp[2]}', file=tpfile)
      print(f'{tp[3]} {tp[4]} {tp[5]}', file=tpfile)
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

   name = np.empty(npl, dtype=np.int32)
   px = np.empty(npl, dtype=np.double)
   py = np.empty(npl, dtype=np.double)
   pz = np.empty(npl, dtype=np.double)
   vx = np.empty(npl, dtype=np.double)
   vy = np.empty(npl, dtype=np.double)
   vz = np.empty(npl, dtype=np.double)
   mass = np.empty(npl, dtype=np.double)
   Gmass = np.empty(npl, dtype=np.double)
   radius = np.empty(npl, dtype=np.double)
   for i, plid in enumerate(plvec):
      name[i] = i + 2
      px[i] = plvec[plid][0]
      py[i] = plvec[plid][1]
      pz[i] = plvec[plid][2]
      vx[i] = plvec[plid][3]
      vy[i] = plvec[plid][4]
      vz[i] = plvec[plid][5]
      Gmass[i] = GMSun * MSun_over_Mpl[plid]**-1
      radius[i] = Rpl[plid] / DU2M
   plfile.write_record(name.T)
   plfile.write_record(px.T)
   plfile.write_record(py.T)
   plfile.write_record(pz.T)
   plfile.write_record(vx.T)
   plfile.write_record(vy.T)
   plfile.write_record(vz.T)
   plfile.write_record(Gmass.T)
   plfile.write_record(radius.T)
   plfile.close()
   tpfile = FortranFile(swiftest_tp, 'w')
   ntp = nclones
   tpfile.write_record(ntp)
   name = np.empty(ntp, dtype=np.int32)
   px = np.empty(ntp, dtype=np.double)
   py = np.empty(ntp, dtype=np.double)
   pz = np.empty(ntp, dtype=np.double)
   vx = np.empty(ntp, dtype=np.double)
   vy = np.empty(ntp, dtype=np.double)
   vz = np.empty(ntp, dtype=np.double)
   for i, tpid in enumerate(tpvec):
      name[i] = int(tpid)
      px[i] = tpvec[tpid][0]
      py[i] = tpvec[tpid][1]
      pz[i] = tpvec[tpid][2]
      vx[i] = tpvec[tpid][3]
      vy[i] = tpvec[tpid][4]
      vz[i] = tpvec[tpid][5]
   tpfile.write_record(name.T)
   tpfile.write_record(px.T)
   tpfile.write_record(py.T)
   tpfile.write_record(pz.T)
   tpfile.write_record(vx.T)
   tpfile.write_record(vy.T)
   tpfile.write_record(vz.T)

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


   sys.stdout = sys.__stdout__



