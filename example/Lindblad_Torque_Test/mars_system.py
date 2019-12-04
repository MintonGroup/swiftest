#Mars and its satellites in cgs units
import numpy as np

#Units will be in terms of planet mass, planet radius, and years
G	    = 6.674e-8                      #Gravitational constant (cgs)
year    = 3600*24*365.25                #seconds in 1 year
AU      = 1.4960e13
M_Sun   = 1.9891e33

a_Mars    = 1.523679 * AU
M_Mars    = 6.41714e26
Rhill_Mars = a_Mars * (M_Mars / (3 * M_Sun))**(1.0 / 3.0)
R_Mars    = 3389.5e5
R_Mars_eq = 3396.2e5
T_Mars    = 1.025957 * 24 * 60 * 60
W_Mars    = 2.0 * np.pi / T_Mars
Ipolar_Mars = 0.365

J2_Mars = 1960e-6
J4_Mars = -19e-6

k2_Mars = 0.164
Q_Mars  = 99.5

a_Phobos    = 9376e5
M_Phobos    = 1.0659e19

a_Deimos    = 23463.2e5
M_Deimos    = 0.14762e19

Sat_r_RM = np.array([a_Phobos, a_Deimos])
Sat_M_Mass = np.array([M_Phobos, M_Deimos])
