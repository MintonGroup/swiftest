from typing import List, Any

import numpy as np
import sys



###***Define initial conditions***###

#Units will be in terms of planet mass, planet radius, and years
G	    = 6.674e-8                          #Gravitational constant (cgs)
year    = 3600*24*365.25                #seconds in 1 year
AU      = 1.4960e13
M_Sun   = 1.9891e33
N   = 199           #number of bins in disk

M_Saturn    = 5.6834e29
R_Saturn    = 60268.0e5
R_Planet = R_Saturn
M_Planet = M_Saturn
a_Saturn    = 9.582172*AU
T_Saturn    = 10.57*3600
W_Saturn    = 2.0*np.pi/T_Saturn
L_Saturn    = 2.0*M_Saturn*R_Saturn**2.0/5.0*W_Saturn
Rhill_Saturn = a_Saturn * (M_Saturn / (3 * M_Sun))**(1.0 / 3.0)

k_2         = 0.104 #tidal love number for primary
Q           = 3000. #tidal dissipation factor for primary
Q_s         = 1.0e-5    #tidal dissipation factor for satellites

rho_sat = 900.0 # Satellite/ring particle mass density in gm/cm**3

J2 = 0.0 #Add these in later
J4 = 0.0

#The following are Unit conversion factors
MU2GM    =     1000.0 #M_Saturn          #Conversion from mass unit to grams
DU2CM    =     100.0 #R_Saturn                       #Conversion from radius unit to centimeters
TU2S     =     1.0 #year                           #Conversion from time unit to seconds
GU       = G / (DU2CM**3 / (MU2GM * TU2S**2))

#Primary body definitions
RP    = R_Saturn / DU2CM #1.0
MP    = M_Saturn / MU2GM #1.0
rhoP  = 3.0 * MP / (4.0 * np.pi * RP**3) #Density of primary
TP    =  T_Saturn / TU2S
IP = 2.0 / 5.0 * MP * RP**2
wP = np.array([0.0,0.0,1.0]) * 2.0 * np.pi / TP
LP =  IP * wP


###For the disk:

r_pdisk = 1.0  #disk particle size
m_pdisk = (4.0 * np.pi*r_pdisk**3)/3.0 * rho_sat   #disk particle size (mass in g)

# gamma	= 0.3	    #ang momentum efficiency factor
inside = 0  #bin id of innermost ring bin (can increase if primary accretes a lot mass through 'Update.py'
r_I	= 100e6      #inside radius of disk is at the embryo's surface
r_F	= 120.47646073197e6  #outside radius of disk
deltar = (r_F - r_I) / N	#width of a bin
deltaX = (2 * np.sqrt(r_F) - 2 * np.sqrt(r_I)) / N  #variable changed bin width used for viscosity calculations
r = []
deltaA = []
m = []
sigma = []
X = []
w = []
I = []
R = []
Torque_to_disk = []

def f(x):
    #Creates initial values for the disk and prints them out
    for a in range(int(N)):
        X.append(2 * np.sqrt(r_I) + deltaX * (a + 0.5))
        r.append((0.5 * X[a])**2)
        deltar = (0.5 * (X[a] + deltaX))**2 - (0.5 * X[a])**2
        deltaA.append(2*np.pi*r[a]*deltar)

        #Power law surface mass density profile
        if x == 0:
            m.append(2*np.pi*alpha/(2-p)*((r[a] + deltar/2)**(2-p) - (r[a] - deltar/2)**(2-p)))
            sigma.append(m[a]/deltaA[a])
        #Gaussian surface mass density profile
        elif x == 1:
            centroid = 110e6    #radius of the center of the gaussian
            spread = 500e3      #width of the gaussian
            sigma_peak = 6e4 #scale factor to get a given mass
            sigma.append(sigma_peak * np.exp(-(r[a]-centroid)**2/(2*spread**2)))
            m.append(sigma[a]*deltaA[a])
        else:
            print('You have not chosen a valid disk model')
        R.append(r[a]**2 + deltar**2/4)
        I.append(m[a]*R[a])
        w.append((GU*MP/r[a]**3)**0.5)
        Torque_to_disk.append(0.0)

f(1) #Make a Gaussian ring


outfile = open('ring.in', 'w')
print(N, file=outfile)
print(r_I, r_F, file=outfile)
print(r_pdisk, GU * m_pdisk, file=outfile)


for a in range(int(N)):
    print(GU * sigma[a],file=outfile)


t_0	= 0
end_sim = 1.1e5 * year / TU2S  #end time
t_print = 1.e2 * year / TU2S #output interval to print results
deltaT	= 1.e2 * year / TU2S  #timestep simulation




plfile = open('pl.in', 'w')
print(1,file=plfile)
print(f'1 {GU*MP}',file=plfile)
print(f'0.0 0.0 0.0',file=plfile)
print(f'0.0 0.0 0.0',file=plfile)
plfile.close()
tpfile = open('tp.in', 'w')
print(0,file=tpfile)
tpfile.close()


iout = int(np.ceil(end_sim / (deltaT * t_print)))
rmin = RP
rmax = Rhill_Saturn / DU2CM

mtiny = 1e-8 #roughly 1/3 the mass of Puck



sys.stdout = open("param.in", "w")
print(f'!Parameter file for the SyMBA-RINGMOONS test')
print(f'!NPLMAX         -1 ')
print(f'!NTPMAX         -1')
print(f'T0             {t_0} ')
print(f'TSTOP          {end_sim}')
print(f'DT             {deltaT}')
print(f'PL_IN          pl.in')
print(f'TP_IN          tp.in')
print(f'IN_TYPE        ASCII')
print(f'ISTEP_OUT      {iout:d}')
print(f'BIN_OUT        bin.dat')
print(f'OUT_TYPE       REAL8')
print(f'OUT_FORM       EL')
print(f'OUT_STAT       NEW')
print(f'ISTEP_DUMP     {iout:d}')
print(f'J2             {J2}')
print(f'J4             {J4}')
print(f'CHK_CLOSE      yes')
print(f'CHK_RMIN       {rmin}')
print(f'CHK_RMAX       {rmax}')
print(f'CHK_EJECT      {rmax}')
print(f'CHK_QMIN       {rmin}')
print(f'CHK_QMIN_COORD HELIO')
print(f'CHK_QMIN_RANGE {rmin} {rmax}')
print(f'ENC_OUT        enc.dat')
print(f'EXTRA_FORCE    no')
print(f'BIG_DISCARD    no')
print(f'RHILL_PRESENT  yes')
print(f'MU2GM          {MU2GM}')
print(f'DU2CM          {DU2CM}')
print(f'TU2S           {TU2S}')
print(f'MTINY          {mtiny}')


sys.stdout = sys.__stdout__

