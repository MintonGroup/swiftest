import numpy as np
import sys

###***Define initial conditions***###

#Units will be in terms of planet mass, planet radius, and years
G	    = 6.674e-8                          #Gravitational constant (cgs)
year    = 3600*24*365.25                #seconds in 1 year
AU      = 1.4960e13
M_Sun   = 1.9891e33


N   = 175            #number of bins in disk

#M_Uranus =  8.6810e28 # mass of Uranus in gm
#R_Uranus = 25362e5 # Radius of Uranus in cm
#T_Uranus = 0.71833 * 24 * 60 * 60 # Rotation period of Uranus in seconds
#W_Uranus    = 2.0 * np.pi / T_Uranus # Angular rotation rate of Uranus
#L_Uranus    = 2.0 * M_Uranus * R_Uranus**2 / 5.0 * W_Uranus  #spin angular momentum of Uranusa
#a_Uranus    = 19.2184 * AU
#Rhill_Uranus = a_Uranus * (M_Uranus / (3 * M_Sun))**(1.0 / 3.0)

M_Saturn    = 5.6846e29
R_Saturn    = 60268.0e5
a_Saturn    = 9.582172*AU
T_Saturn    = 10.57*3600
W_Saturn    = 2.0*pi/T_Saturn
L_Saturn    = 2.0*M_Saturn*R_Saturn**2.0/5.0*W_Saturn
Rhill_Saturn = a_Saturn * (M_Saturn / (3 * M_Sun))**(1.0 / 3.0)


k_2         = 0.104 #tidal love number for primary
Q           = 3000. #tidal dissipation factor for primary
Q_s         = 1.0e-5    #tidal dissipation factor for satellites

rho_sat = 1.2 # Satellite/ring particle mass density in gm/cm**3

J2 = 0.0 #Add these in later
J4 = 0.0

#The following are Unit conversion factors
MU2GM    =      M_Uranus          #Conversion from mass unit to grams
DU2CM    =      R_Uranus                       #Conversion from radius unit to centimeters
TU2S     =      year                           #Conversion from time unit to seconds
GU       = G / (DU2CM**3 / (MU2GM * TU2S**2))


#Primary body definitions
RP    = 1.0
MP    = 1.0
rhoP  = 3.0 * MP / (4.0 * np.pi * RP**3) #Density of primary
TP    =  T_Uranus / TU2S
IP = 2.0 / 5.0 * MP * RP**2
wP = np.array([0.0,0.0,1.0]) * 2.0 * np.pi / TP
LP =  IP * wP








###For the disk:

r_p100m = 100e2                        #radius of 100 m planetesimal (cm)
m_p100m = (4.0 * np.pi * r_p100m**3.0) / 3.0 * rho_sat     #mass of 100 m planetesimal (g)

r_pdisk = r_p100m    #disk particle size (radius)
m_pdisk = m_p100m  #disk particle size (mass)

# gamma	= 0.3	    #ang momentum efficiency factor
inside = 0  #bin id of innermost ring bin (can increase if primary accretes a lot mass through 'Update.py'
p	=   3   #Power law index that surface mass density falls off as (cannot be exactly 2)
alpha	=  4.12e32      #Arbitrary constant (sets initial mass of disk in "build.py"
r_I	= R_Uranus     #inside radius of disk is at the embryo's surface
r_F	= 5.0 * R_Uranus  #outside radius of disk
deltar = (r_F - r_I) / N	#width of a bin
deltaX = (2. * r_F**0.5 - 2.*r_I**0.5) / N  #variable changed bin width used for viscosity calculations
r = []
deltaA = []
m = []
sigma = []

def f(x):
    #Creates initial values for the disk and prints them out
    for a in range(int(N)):
        r.append(r_I + deltar*(a+0.5))
        deltaA.append(2*np.pi*deltar*(r_I + deltar*(a + 0.5)))

        #Power law surface mass density profile
        if x == 0:
            m.append(2*np.pi*alpha/(2-p)*((r[a] + deltar/2)**(2-p) - (r[a] - deltar/2)**(2-p)))
            sigma.append(m[a]/deltaA[a])
        #Gaussian surface mass density profile
        elif x == 1:
            centroid = 20       #bin id of the center of the gaussian
            spread = 1          #width of the gaussian
            mass_scale = 2.3e6  #scale factor to get a given mass
            sigma.append(mass_scale/(5*(2*np.pi)**(0.5))*np.exp(-(a-centroid)**2/(2*spread**2)))
            m.append(sigma[a]*deltaA[a])
        else:
            print('You have not chosen a valid disk model')

f(1) #Make a Gaussian ring

outfile = open('ring.in', 'w')
print(N, file=outfile)
print(deltar / DU2CM, file=outfile)
print(r_pdisk / DU2CM, G * (m_pdisk / MU2GM), file=outfile)


for a in range(int(N)):
    #print(sigma[a] * DU2CM**2 / MU2GM,file=outfile)
    print(sigma[a] ,file=outfile)







t_0	= 0
end_sim = 5.e9  #end time
t_print = 2.e4  #output interval to print results
deltaT	= 1.e2  #timestep simulation




plfile = open('pl.in', 'w')
print(1,file=plfile)
print(f'1 {G*MP}',file=plfile)
print(f'0.0 0.0 0.0',file=plfile)
print(f'0.0 0.0 0.0',file=plfile)
plfile.close()
tpfile = open('tp.in', 'w')
print(0,file=tpfile)
tpfile.close()


iout = int(np.ceil(end_sim / (deltaT * t_print)))
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




