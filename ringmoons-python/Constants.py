from __future__ import division
from math import *
###***Physical constants***###

G	= 6.674e-8                          #Gravitational constant (cgs)
M_Sun   = 1.9891e33
AU      = 1.4960e13
M_Earth = 5.9729e27
a_Mars  = 1.523679*AU
M_Mars	= 6.4185e26                     #(final) grams
R_Mars	= 3.394e8                       #(final) cm
# rho_sat = (1.471 + 1.876)/2.0
rho_sat = 1.2
T_Mars  = 24.622962*3600                #Rotation period of Mars(s)
W_Mars  = 2*pi/T_Mars                   #Angular rotation of Mars
L_Mars  = 2*M_Mars*R_Mars**2/5*W_Mars   #Angular Momentum of Mars
a_Phobos = 9.376e8
M_Phobos = 1.0659e19
R_Phobos = 1.112667e6
a_Deimos = 2.3463e9
M_Deimos = 1.4762e18
R_Deimos = 6.2e5

M_Saturn    = 5.6846e29
R_Saturn    = 60268.0e5
a_Saturn    = 9.582172*AU
T_Saturn    = 10.57*3600
W_Saturn    = 2.0*pi/T_Saturn
L_Saturn    = 2.0*M_Saturn*R_Saturn**2.0/5.0*W_Saturn

M_Uranus    = 8.6810e28
R_Uranus    = 25362e5
a_Uranus    = 19.2184*AU
T_Uranus    = 0.71833*24*3600
W_Uranus    = 2.0*pi/T_Uranus
L_Uranus    = 2.0*M_Uranus*R_Uranus**2.0/5.0*W_Uranus

M_Neptune    = 1.0243e29
R_Neptune    = 24622e5
a_Neptune    = 30.110387*AU
T_Neptune    = 0.6713*24*3600
W_Neptune    = 2.0*pi/T_Neptune
L_Neptune    = 2.0*M_Neptune*R_Neptune**2.0/5.0*W_Neptune

# v_K     = (G*M_Sun/a_Mars)**(1/2)
year    = 3600*24*365.25                #seconds in 1 year
r_p1    = 1e2                           #radius of 1 m planetesimal (cm)
m_p1    = (4.0*pi*r_p1**3.0)/3.0*rho_sat        #mass of 1 m planetesimal
r_p10   = 1.e3                          #radius of 10m particle (cm)
m_p10    = (4.0*pi*r_p10**3.0)/3.0*rho_sat        #mass of 10 m planetesimal
r_p500m = 500e2                         #radius of 500 m planetesimal (cm)
m_p500m = (4.0*pi*r_p500m**3.0)/3.0*rho_sat     #mass of 500 m planetesimal (g)
r_p100m = 100e2                         #radius of 100 m planetesimal (cm)
m_p100m = (4.0*pi*r_p100m**3.0)/3.0*rho_sat     #mass of 100 m planetesimal (g)
r_p1k    = 1e5                          #radius of 1 km planetesimal (cm)
m_p1k    = (4.0*pi*r_p1k**3.0)/3.0*rho_sat      #mass of 1 km planetesimal (g)
r_p10k    = 1e6                         #radius of 10 km planetesimal (cm)
m_p10k    = (4.0*pi*r_p10k**3.0)/3.0*rho_sat    #mass of 10 km planetesimal (cm)
r_p100k     = 1e7                       #radius of 100 km planetesimal (cm)
m_p100k     = (4.0*pi*r_p100k**3.0)/3.0*rho_sat #mass of 100 km planetesimal (g)
r_p1000k    = 1e8                       #mass of 1000 km planetesimal (g)
m_p1000k    = (4.0*pi*r_p1000k**3.0)/3.0*rho_sat #radius of 1000 km planetesimal (g)
