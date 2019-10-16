from Constants import *
from Lists import *
import csv
###***Define initial conditions***###

R_Planet    = R_Uranus  #radius of primary
M_Planet    = M_Uranus  #mass of primary
rho     = 3.0*M_Planet/(4.0*pi*R_Planet**3) #final density of Mars
L           = L_Uranus  #spin angular momentum of primary

k_2         = 0.104 #tidal love number for primary
Q           = 3000. #tidal dissipation factor for primary
Q_s         = 1.0e-5    #tidal dissipation factor for satellites
Sat_e_init  = 1.0e-5    #Initial eccentricity of satellites
###For the disk:
# r_p = r_p100k       #Set impactor size (radius)
# m_p = m_p100k       #Set impactor size (mass)
# A_p = pi*r_p**2.0   #Cross sectional area of impactor
r_pdisk = r_p100m    #disk particle size (radius)
m_pdisk = m_p100m  #disk particle size (mass)
deltaT	= 1.e2*year       #timestep simulation
#interval = 1e5       #number of iterations before Update and Restructure are run
p	=   3   #Order that surface mass density falls off as (cannot be exactly 2)
alpha	=  4.12e32      #Arbitrary constant (sets initial mass of disk in "build.py"
t_0	= 0	            #Set initial time to zero
N   = 175            #number of bins in disk
r_F	= 5.0*R_Planet  #outside radius of disk

# gamma	= 0.3	    #ang momentum efficiency factor
Mars_Accrete = 0.0  #Mass that passes inside the disk to Mars, initialize at zero
inside = 0  #bin id of innermost ring bin (can increase if primary accretes a lot mass through 'Update.py'


###For the embryo:
# Sig = 8*(1.5**(-1.5)) #g/cm**2
# pm = 1.5 #g/cm**3
# pM = 4 #g/cm**3
# C = 6*pi*(.75*(1/pM))**(1./3.)*(G/M_Sun)**(.5)
# b = 9.8
# t=0     #initial time
# M_init = 1.68e25 #g Initial mass of the embryo
M = M_Planet #M_Uranus
R_Embryo = R_Planet #(3.0*M/(4.0*pi*rho))**(1.0/3.0)



# Ntotal = (M_Mars - M)/m_p #Total number of collisions in Mars growth
r_I	= R_Embryo          #inside radius of disk is at the embryo's surface
deltar = (r_F - r_I)/N	#width of a bin
deltaX = (2.*r_F**0.5 - 2.*r_I**0.5)/N  #variable changed bin width used for viscosity calculations
t_print = 2.e4       #output interval (in years) to print results