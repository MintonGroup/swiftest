#Saturn and its satellites in cgs units
import numpy as np

#Units will be in terms of planet mass, planet radius, and years
G	    = 6.674e-8                      #Gravitational constant (cgs)
year    = 3600*24*365.25                #seconds in 1 year
AU      = 1.4960e13
M_Sun   = 1.9891e33

a_Saturn    = 9.582172 * AU
M_Saturn    = 5.6834e29
Rhill_Saturn = a_Saturn * (M_Saturn / (3 * M_Sun))**(1.0 / 3.0)
R_Saturn    = 58232e5
R_Saturn_eq = 60268e5
T_Saturn    = 10.57 * 3600
W_Saturn    = 2.0 * np.pi / T_Saturn
Ipolar_Saturn = 0.210
#Anderson & Schubert post Cassini values
J2_Saturn = 16290.71e-6
J4_Saturn = -935.5e-6

#Lainey et al. (2017)
k2_Saturn = 0.390
Q_Saturn  = 2450.0

a_Pan        = 133584e5
M_Pan        = 4.95e18

a_Daphnis    = 136505e5
M_Daphnis    = 0.08e18

a_Atlas      = 137670e5
M_Atlas      = 6.6e18

a_Prometheus = 139380e5
M_Prometheus = 159.5e18

a_Pandora    = 141720e5
M_Pandora    = 137.1e18

a_Epimetheus = 151422e5
M_Epimetheus = 526.0e18

a_Janus      = 151472e5
M_Janus      = 1897.0e18

a_Mimas      = 185404e5
M_Mimas      = 37493.0e18

a_Methone    = 194440e5
M_Methone    = 0.02e18

a_Anthe      = 197700e5
M_Anthe      = 0.0015e18

a_Pallene    = 212280e5
M_Pallene    = 0.05e18

a_Enceladus  = 237950e5
M_Enceladus  = 108022e18

a_Tethys     = 294619e5
M_Tethys     = 61744e18

Sat_r_RM = np.array([a_Pan, a_Daphnis, a_Atlas, a_Prometheus, a_Pandora, a_Epimetheus, a_Janus, a_Mimas, a_Methone, a_Anthe, a_Pallene, a_Enceladus, a_Tethys])
Sat_M_Mass = np.array([M_Pan, M_Daphnis, M_Atlas, M_Prometheus, M_Pandora, M_Epimetheus, M_Janus, M_Mimas, M_Methone, M_Anthe, M_Pallene, M_Enceladus, M_Tethys])
