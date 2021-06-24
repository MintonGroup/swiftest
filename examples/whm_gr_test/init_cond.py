import numpy as np
import sys
from astroquery.jplhorizons import Horizons
import astropy.constants as const 

#Values from JPL Horizons
AU2M = const.au.value
GMSunSI = const.GM_sun.value
Rsun = const.R_sun.value
GC = const.G.value
JD = 86400
year = 365.25 * JD
c = 299792458.0
MSun_over_Mpl = [6023600.0,
                 408523.71,
                 328900.56,
                 3098708.,
                 1047.3486,
                 3497.898,
                 22902.98,
                 19412.24,
                 1.35e8]

MU2KG    = GMSunSI / GC #Conversion from mass unit to kg
DU2M     = AU2M         #Conversion from radius unit to centimeters
TU2S     = year         #Conversion from time unit to seconds
GU       = GC / (DU2M**3 / (MU2KG * TU2S**2))

GMSun = GMSunSI / (DU2M**3 / TU2S**2)

t_print = 10.e0 * year / TU2S #output interval to print results
deltaT	= 0.25 * JD / TU2S    #timestep simulation
end_sim = 1.0e3 * year / TU2S + t_print #end time

# Solar oblatenes values: From Mecheri et al. (2004), using Corbard (b) 2002 values (Table II)
J2 = 2.198e-7 * (Rsun / DU2M)**2
J4 = -4.805e-9 * (Rsun / DU2M)**4

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
   'mercury'     : 6023600.0,
   'venus'       : 408523.71,
   'earthmoon'   : 328900.56,
   'mars'        : 3098708.,
   'jupiter'     : 1047.3486,
   'saturn'      : 3497.898,
   'uranus'      : 22902.98,
   'neptune'     : 19412.24,
   'plutocharon' :  1.35e8
}

#Planet radii in meters
Rpl = {
   'mercury'     : 2439.4e3,
   'venus'       : 6051.8e3,
   'earthmoon'   : 6371.0084e3, # Earth only for radius
   'mars'        : 3389.50e3,
   'jupiter'     : 69911e3,
   'saturn'      : 58232.0e3,
   'uranus'      : 25362.e3,
   'neptune'     : 24622.e3,
   'plutocharon' : 1188.3e3
}

pdata = {}
plvec = {}
Rhill = {}

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
   Rhill[key] = pdata[key].elements()['a'][0] * (3 * MSun_over_Mpl[key])**(-1.0 / 3.0)


if __name__ == '__main__':
   # Convert from AU-day to AU-year just because I find it easier to keep track of the sim progress
   for plid in plvec:
      plvec[plid][3:] *= year / JD

   # Names of all output files
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

   # Simulation start, stop, and output cadence times
   t_0	  = 0 # simulation start time
   end_sim = 1000.0e0 * year / TU2S # simulation end time
   deltaT	= 0.25 * JD / TU2S   # simulation step size
   t_print = 1.0 * year / TU2S #output interval to print results

   iout = int(np.ceil(t_print / deltaT))
   rmin = Rsun / DU2M
   rmax = 1000.0

   #Make Swifter files
   plfile = open(swifter_pl, 'w')
   print(f'{npl+1} ! Planet input file generated using init_cond.py using JPL Horizons data for the major planets (and Pluto) for epoch {tstart}' ,file=plfile)
   print(f'1 {GMSun}',file=plfile)
   print(f'0.0 0.0 0.0',file=plfile)
   print(f'0.0 0.0 0.0',file=plfile)
   for i, plid in enumerate(plvec): 
      print(f'{i + 2} {GMSun / MSun_over_Mpl[plid]} {Rhill[plid]}', file=plfile)
      print(f'{Rpl[plid] / DU2M}', file=plfile)
      print(f'{plvec[plid][0]} {plvec[plid][1]} {plvec[plid][2]}', file=plfile)
      print(f'{plvec[plid][3]} {plvec[plid][4]} {plvec[plid][5]}', file=plfile)
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
   print(f'OUT_FORM      EL')
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
   print(f'C              {c / (DU2M / TU2S)}')

   #Now make Swiftest files
   cbfile = open(swiftest_cb, 'w')
   print(f'{1.0}',file=cbfile)
   print(f'{rmin}',file=cbfile)
   print(f'{J2}',file=cbfile)
   print(f'{J4}',file=cbfile)

   plfile = open(swiftest_pl, 'w')
   print(npl,file=plfile)

   for i, plid in enumerate(plvec): 
      print(f'{i + 2} {1.0 / MSun_over_Mpl[plid]}', file=plfile)
      print(f'{Rpl[plid] / DU2M}', file=plfile)
      print(f'{plvec[plid][0]} {plvec[plid][1]} {plvec[plid][2]}', file=plfile)
      print(f'{plvec[plid][3]} {plvec[plid][4]} {plvec[plid][5]}', file=plfile)
   plfile.close()
   tpfile = open(swiftest_tp, 'w')
   print(0,file=tpfile)
   tpfile.close()

   sys.stdout = open(swiftest_input, "w")
   print(f'! Swiftest input file generated using init_cond.py')
   print(f'T0             {t_0} ')
   print(f'TSTOP          {end_sim}')
   print(f'DT             {deltaT}')
   print(f'CB_IN          {swiftest_cb}')
   print(f'PL_IN          {swiftest_pl}')
   print(f'TP_IN          {swiftest_tp}')
   print(f'IN_TYPE        ASCII')
   print(f'ISTEP_OUT      {iout:d}')
   print(f'ISTEP_DUMP     {iout:d}')
   print(f'BIN_OUT        {swiftest_bin}')
   print(f'OUT_TYPE       REAL8')
   print(f'OUT_FORM       EL')
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
   print(f'GR             yes')
   print(f'MU2KG          {MU2KG}')
   print(f'DU2M           {DU2M}')
   print(f'TU2S           {TU2S}')


   sys.stdout = sys.__stdout__
