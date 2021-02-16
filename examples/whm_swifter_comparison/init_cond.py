import numpy as np
import sys
from astroquery.jplhorizons import Horizons
import astropy.constants as const 
from scipy.io import FortranFile

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
end_sim = 1 * year  / TU2S # simulation end time
t_print = deltaT #year / TU2S #output interval to print results


# Solar oblatenes values: From Mecheri et al. (2004), using Corbard (b) 2002 values (Table II)
J2 = 0.0 #np.longdouble(2.198e-7) * (Rsun / DU2M)**2
J4 = 0.0 #np.longdouble(-4.805e-9) * (Rsun / DU2M)**4

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

asteroidid = {
   '100001'  : 'Ceres',
   '100002'  : 'Pallas',
   '100003'  : 'Juno',
   '100004'  : 'Vesta'
}
ntp = 4

tdata = {}
tpvec = {}
for key,val in asteroidid.items():
   tdata[key] = Horizons(id=val, id_type='smallbody', location='@sun',
            epochs={'start': tstart, 'stop': tend,
            'step': tstep})
   tpvec[key] = np.array([tdata[key].vectors()['x'][0],
         tdata[key].vectors()['y'][0], 
         tdata[key].vectors()['z'][0], 
         tdata[key].vectors()['vx'][0], 
         tdata[key].vectors()['vy'][0], 
         tdata[key].vectors()['vz'][0] 
         ])


if __name__ == '__main__':
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
   rmax = np.longdouble(1000.0)
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
   print(ntp,file=tpfile)
   for tpid, tp in tpvec.items():
      print(tpid, file=tpfile)
      print(tp[0],tp[1],tp[2], file=tpfile)
      print(tp[3],tp[4],tp[5], file=tpfile)
   tpfile.close()

   sys.stdout = open(swifter_input, "w")
   print('! Swifter input file generated using init_cond.py')
   print('T0            ',t_0)
   print('TSTOP         ',end_sim)
   print('DT            ',deltaT)
   print('PL_IN         ',swifter_pl)
   print('TP_IN         ',swifter_tp)
   print('IN_TYPE       ASCII')
   print('ISTEP_OUT     ',iout)
   print('ISTEP_DUMP    ',iout)
   print('BIN_OUT       ',swifter_bin)
   print('OUT_TYPE      REAL8')
   print('OUT_FORM      XV')
   print('OUT_STAT      NEW')
   print('J2            ',J2)
   print('J4            ',J4)
   print('CHK_CLOSE     yes')
   print('CHK_RMIN      ',rmin)
   print('CHK_RMAX      ',rmax)
   print('CHK_EJECT     ',rmax)
   print('CHK_QMIN      ',rmin)
   print('CHK_QMIN_COORD HELIO')
   print('CHK_QMIN_RANGE ',rmin,rmax)
   print('ENC_OUT        ',swifter_enc)
   print('EXTRA_FORCE    no')
   print('BIG_DISCARD    no')
   print('RHILL_PRESENT  yes')

   sys.stdout = sys.__stdout__
   #Now make Swiftest files
   #cbfile = open(swiftest_cb, 'w')
   cbfile = FortranFile(swiftest_cb, 'w')
   #print(1.0,file=cbfile)
   #print(rmin,file=cbfile)
   #print(J2,file=cbfile)
   #print(J4,file=cbfile)
   Msun = np.double(1.0)
   cbfile.write_record(np.double(GMSun))
   cbfile.write_record(np.double(rmin))
   cbfile.write_record(np.double(J2))
   cbfile.write_record(np.double(J4))
   cbfile.close()

   #plfile = open(swiftest_pl, 'w')
   plfile = FortranFile(swiftest_pl, 'w')
   #print(npl,file=plfile)
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
   #for i, plid in enumerate(plvec): 
   #  print(i + 2,"{:.23g}".format(np.longdouble(MSun_over_Mpl[plid]**-1)), file=plfile)
   #  print(Rpl[plid] / DU2M, file=plfile)
   #  print(plvec[plid][0], plvec[plid][1], plvec[plid][2], file=plfile)
   #  print(plvec[plid][3], plvec[plid][4], plvec[plid][5], file=plfile)
   plfile.close()
   #tpfile = open(swiftest_tp, 'w')
   tpfile = FortranFile(swiftest_tp, 'w')
   #print(ntp,file=tpfile)
   tpfile.write_record(ntp)
   #for tpid, tp in tpvec.items():
   #   print(tpid, file=tpfile)
   #   print(tp[0],tp[1],tp[2], file=tpfile)
   #   print(tp[3],tp[4],tp[5], file=tpfile)

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
   print('! Swiftest input file generated using init_cond.py')
   print('T0             ',t_0)
   print('TSTOP          ',end_sim)
   print('DT             ',deltaT)
   print('CB_IN          ',swiftest_cb)
   print('PL_IN          ',swiftest_pl)
   print('TP_IN          ',swiftest_tp)
   print('IN_TYPE        REAL8')
   print('ISTEP_OUT      ',iout)
   print('ISTEP_DUMP     ',iout)
   print('BIN_OUT        ',swiftest_bin)
   print('OUT_TYPE       REAL8')
   print('OUT_FORM       XV')
   print('OUT_STAT       REPLACE')
   print('CHK_CLOSE      yes')
   print('CHK_RMIN       ',rmin)
   print('CHK_RMAX       ',rmax)
   print('CHK_EJECT      ',rmax)
   print('CHK_QMIN       ',rmin)
   print('CHK_QMIN_COORD HELIO')
   print('CHK_QMIN_RANGE ',rmin,rmax)
   print('ENC_OUT        ',swiftest_enc)
   print('EXTRA_FORCE    no')
   print('BIG_DISCARD    no')
   print('ROTATION       no')
   print('GR             no')
   print('MU2KG          ',MU2KG)
   print('DU2M           ',DU2M)
   print('TU2S           ',TU2S)


   sys.stdout = sys.__stdout__
