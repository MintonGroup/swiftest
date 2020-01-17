import numpy as np
import sys
from uranus_system import *



#The following are Unit conversion factors
MU2GM    =     M_Uranus          #Conversion from mass unit to grams
DU2CM    =     R_Uranus                       #Conversion from radius unit to centimeters
TU2S     =     year                           #Conversion from time unit to seconds
GU       = G / (DU2CM**3 / (MU2GM * TU2S**2))

r_pdisk = 100.0e2 / DU2CM #disk particle size
rho_pdisk = 1.2 * DU2CM**3 / MU2GM # Satellite/ring particle mass density in gm/cm**3
rho_sat   = rho_pdisk # Satellite/ring particle mass density in gm/cm**3


t_0	= 0
t_print = 1e3 * year / TU2S #output interval to print results
deltaT	= 1e3 * year / TU2S  #timestep simulation
end_sim = 1.0e6 * year / TU2S + t_print #end time
feeding_zone_factor = 2.2
rkf_tol = 1e-8

Nbins    = 256      #number of bins in disk



sigma_peak = 0.6e4 * DU2CM ** 2 / MU2GM  # scale factor to get a given mass


k_2         = k2_Uranus #tidal love number for primary
Q           = Q_Uranus_TWlo #tidal dissipation factor for primary
Q_s         = 1.0e-5    #tidal dissipation factor for satellites

J2 = J2_Uranus
J4 = J4_Uranus

#Primary body definitions
RP    = R_Uranus / DU2CM  # Radius of primary
RP_eq = R_Uranus_eq / DU2CM
MP    = M_Uranus / MU2GM  # Mass of primary
rhoP  = 3.0 * MP / (4.0 * np.pi * RP**3) #Density of primary
TP    =  T_Uranus / TU2S # Rotation rate of primary
IPp = Ipolar_Uranus # Polar moment of inertia of primary
IPe = J2 + IPp # equatorial moment of inertia of primary
FRL = 2.456 * RP * (rhoP / rho_pdisk)**(1./3.)
RRL = 1.44 * RP * (rhoP / rho_sat)**(1./3.)
Rsync = (GU * MP * TP**2 / (4 * np.pi**2))**(1./3.)

r_I	= 0.99 * RP      #inside radius of disk is at the embryo's surface
r_F	= 2.0 * FRL  #outside radius of disk

wP = np.array([0.0,0.0,1.0]) * 2.0 * np.pi / TP # rotation vector of primary
IP = np.array([IPe, IPe, IPp]) # Principal moments of inertia

###For the disk:
m_pdisk = (4.0 * np.pi*r_pdisk**3)/3.0 * rho_pdisk   #disk particle size (mass in g)

# gamma	= 0.3	    #ang momentum efficiency factor
inside = 0  #bin id of innermost ring bin (can increase if primary accretes a lot mass through 'Update.py'

deltar = (r_F - r_I) / Nbins	#width of a bin
deltaX = (2 * np.sqrt(r_F) - 2 * np.sqrt(r_I)) / Nbins  #variable changed bin width used for viscosity calculations
r = []
deltaA = []
m = []
sigma = []
X = []
w = []
I = []
R = []
Torque_to_disk = []

def f():

    #Creates initial values for the disk and prints them out
    for a in range(int(Nbins)):
        X.append(2 * np.sqrt(r_I) + deltaX * (a + 0.5))
        r.append((0.5 * X[a])**2)
        deltar = (0.5 * (X[a] + deltaX))**2 - (0.5 * X[a])**2
        deltaA.append(2*np.pi*r[a]*deltar)

        #if a >= Nbins - 45:
        #if r[a] > FRL:
        sigma.append(0.0)
        #else:
        #    sigma.append(sigma_peak * (r[a] / RP) ** (-3))
        m.append(sigma[a] * deltaA[a])
        R.append(r[a]**2 + deltar**2/4)
        I.append(m[a]*R[a])
        w.append((GU*MP/r[a]**3)**0.5)
        Torque_to_disk.append(0.0)

if __name__ == '__main__':
    f() #Make a power law ring

    Nseeds   = 1
    Gmseed = [4 * M_Mira * GU / MU2GM]
    aseed= [1.1* RRL]


    outfile = open('ring.in', 'w')
    print(Nbins, Nseeds, feeding_zone_factor, rkf_tol, file=outfile)
    print(r_I, r_F, file=outfile)
    print(r_pdisk, GU * m_pdisk, file=outfile)

    for a in range(int(Nbins)):
        print(GU * sigma[a],file=outfile)

    for a in range(int(Nseeds)):
        print(aseed[a], Gmseed[a], file=outfile)

    plfile = open('pl.in', 'w')
    print(1,file=plfile)
    print(f'1 {GU*MP}',file=plfile)
    print(f'{RP} {k_2} {Q}', file=plfile)
    np.savetxt(plfile, IP, newline=' ')
    print('',file=plfile)
    np.savetxt(plfile, wP, newline=' ')
    print('',file=plfile)
    print(f'0.0 0.0 0.0',file=plfile)
    print(f'0.0 0.0 0.0',file=plfile)
    plfile.close()
    tpfile = open('tp.in', 'w')
    print(0,file=tpfile)
    tpfile.close()


    iout = int(np.ceil(t_print / deltaT))
    rmin = RP
    rmax = Rhill_Uranus / DU2CM

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
    print(f'RING_OUTFILE   ring.dat')
    print(f'ROTATION       yes')
    print(f'PREDPREY       no')


    sys.stdout = sys.__stdout__

