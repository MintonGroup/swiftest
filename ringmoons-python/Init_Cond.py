import numpy as np
###***Define initial conditions***###
N   = 175            #number of bins in disk

M_Uranus =  8.6810e28 # mass of Uranus in gm
R_Uranus = 25362e5 # Radius of Uranus in cm
T_Uranus = 0.71833 * 24 * 60 * 60 # Rotation period of Uranus in seconds
W_Uranus    = 2.0 * np.pi / T_Uranus # Angular rotation rate of Uranus
L_Uranus    = 2.0 * M_Uranus * R_Uranus**2.0 / 5.0 * W_Uranus  #spin angular momentum of Uranus

k_2         = 0.104 #tidal love number for primary
Q           = 3000. #tidal dissipation factor for primary
Q_s         = 1.0e-5    #tidal dissipation factor for satellites
Sat_e_init  = 1.0e-5    #Initial eccentricity of satellites

rho_sat = 1.2 # Satellite/ring particle mass density in gm/cm**3


#Primary body definitions
RP    =  R_Uranus
MP    =  M_Uranus
rhoP     = 3.0 * MP / (4.0 * np.pi * RP**3) #Density of primary
LP = L_Uranus

p	=   3   #Power law index that surface mass density falls off as (cannot be exactly 2)
alpha	=  4.12e32      #Arbitrary constant (sets initial mass of disk in "build.py"
r_F	= 5.0 * RP  #outside radius of disk

#Units will be in terms of planet mass, planet radius, and years
G	= 6.674e-8                          #Gravitational constant (cgs)
year    = 3600*24*365.25                #seconds in 1 year
#The following should match the Swifter param.in file
MU2GM    =      M_Uranus          #Conversion from mass unit to grams
DU2CM    =      R_Uranus                       #Conversion from radius unit to centimeters
TU2S     =      year                           #Conversion from time unit to seconds




###For the disk:

r_p100m = 100e2                         #radius of 100 m planetesimal (cm)
m_p100m = (4.0 * np.pi * r_p100m**3.0) / 3.0 * rho_sat     #mass of 100 m planetesimal (g)

r_pdisk = r_p100m    #disk particle size (radius)
m_pdisk = m_p100m  #disk particle size (mass)

# gamma	= 0.3	    #ang momentum efficiency factor
inside = 0  #bin id of innermost ring bin (can increase if primary accretes a lot mass through 'Update.py'


# Ntotal = (M_Mars - M)/m_p #Total number of collisions in Mars growth
r_I	= RP     #inside radius of disk is at the embryo's surface
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

outfile = open('ring.in', 'w')
print(N, deltar / RP, file=outfile)
print(r_pdisk, m_pdisk, file=outfile)

f(1)

for a in range(int(N)):
    print(sigma[a],file=outfile)
