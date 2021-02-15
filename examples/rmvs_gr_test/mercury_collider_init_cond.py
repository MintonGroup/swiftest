"""
This script generates initial conditions for the solar using JPL Horizons data. 
Makes a test particle located at Mercury's position, 500 years into the future so that I can model an impact onto Mercury.

To use the script, modify the variables just after the  "if __name__ == '__main__':" line
"""
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



MU2KG   = GMSunSI / GC #Conversion from mass unit (G * Msun) to kg
DU2M   = AU2M      #Conversion from distance unit (AU) to meters
TU2S   = JD       #Conversion from time unit (Julian Day) to seconds
GU     = GC / (DU2M**3 / (MU2KG * TU2S**2))

GMSun = GMSunSI / (DU2M**3 / TU2S**2)



# Solar oblateness values: From Mecheri et al. (2004), using Corbard (b) 2002 values (Table II)
J2 = 2.198e-7 * (Rsun / DU2M)**2
J4 = -4.805e-9 * (Rsun / DU2M)**4

npl = 9

tstart = '2421-01-31'
tend = '2421-02-01'
tstep = '1d'

# All planets except for Mercury
planetid = {
   'venus'      : '2',
   'earthmoon'   : '3',
   'mars'      : '4',
   'jupiter'    : '5',
   'saturn'     : '6',
   'uranus'     : '7',
   'neptune'    : '8',
   'plutocharon' : '9'
}

tpid = {'mercury_impactor' : '1'}

#Planet Msun/M ratio
MSun_over_Mpl = {
   'mercury'    : 6023600.0,
   'venus'      : 408523.71,
   'earthmoon'   : 328900.56,
   'mars'      : 3098708.,
   'jupiter'    : 1047.3486,
   'saturn'     : 3497.898,
   'uranus'     : 22902.98,
   'neptune'    : 19412.24,
   'plutocharon' :  1.35e8
}

#Planet radii in meters
Rpl = {
   'mercury'    : 2439.4e3,
   'venus'      : 6051.8e3,
   'earthmoon'   : 6371.0084e3, # Earth only for radius
   'mars'      : 3389.50e3,
   'jupiter'    : 69911e3,
   'saturn'     : 58232.0e3,
   'uranus'     : 25362.e3,
   'neptune'    : 24622.e3,
   'plutocharon' : 1188.3e3
}

pdata = {}
vec = {}
Rhill = {}

for key,val in planetid.items():
   pdata[key] = Horizons(id=val, id_type='majorbody',location='@sun',
            epochs={'start': tstart, 'stop': tend,
            'step': tstep})
   vec[key] = np.array([pdata[key].vectors()['x'][0],
               pdata[key].vectors()['y'][0], 
               pdata[key].vectors()['z'][0], 
               pdata[key].vectors()['vx'][0], 
               pdata[key].vectors()['vy'][0], 
               pdata[key].vectors()['vz'][0]
         ])
   Rhill[key] = pdata[key].elements()['a'][0] * (3 * MSun_over_Mpl[key])**(-1.0 / 3.0)

# Make a test particle initially at Mercury's position but with 10% higher velocity
vfactor = 1.00
tpdata = Horizons(id='1', id_type='majorbody',location='@sun',
            epochs={'start': tstart, 'stop': tend,
            'step': tstep}) 
tpvec = [tpdata.vectors()['x'][0], 
         tpdata.vectors()['y'][0], 
         tpdata.vectors()['z'][0],
         vfactor * tpdata.vectors()['vx'][0], 
         vfactor * tpdata.vectors()['vy'][0],
         vfactor * tpdata.vectors()['vz'][0]
      ]

if __name__ == '__main__':
   # Names of all output files
   swifter_input  = "param.tpcollider.in"
   swifter_pl     = "pl.tpcollider.in"
   swifter_tp     = "tp.tpcollider.in"
   swifter_bin    = "bin.tpcollider.dat"
   swifter_enc    = "enc.tpcollider.dat"

   # Simulation start, stop, and output cadence times
   t_0	  = 0 * year / TU2S # simulation start time
   deltaT	= 0.1 * JD / TU2S   # simulation step size
   t_print = 4.e0 * year / TU2S #output interval to print results
   end_sim = 400.0 * year / TU2S # simulation end time


   iout = int(np.ceil(t_print / deltaT))
   rmin = Rsun / DU2M
   rmax = 10000.0

   #Make Swifter files
   # Reverse all the velocity signs to fake going backward in time
   plfile = open(swifter_pl, 'w')
   print(f'{npl} ! Planet input file generated using init_cond.py using JPL Horizons data for the major planets (and Pluto) for epoch {tstart}' ,file=plfile)
   print(f'1 {GMSun}',file=plfile)
   print(f'0.0 0.0 0.0',file=plfile)
   print(f'0.0 0.0 0.0',file=plfile)
   for i, key in enumerate(planetid): 
      print(f'{i + 3} {GMSun / MSun_over_Mpl[key]} {Rhill[key]}', file=plfile)
      print(f'{Rpl[key] / DU2M}', file=plfile)
      print(f'{vec[key][0]} {vec[key][1]} {vec[key][2]}', file=plfile)
      print(f'{-vec[key][3]} {-vec[key][4]} {-vec[key][5]}', file=plfile)
   plfile.close()
   tpfile = open(swifter_tp, 'w')
   print(1,file=tpfile)
   print(100,file=tpfile)
   print(f'{tpvec[0]} {tpvec[1]} {tpvec[2]}', file=tpfile)
   print(f'{-tpvec[3]} {-tpvec[4]} {-tpvec[5]}', file=tpfile)
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



