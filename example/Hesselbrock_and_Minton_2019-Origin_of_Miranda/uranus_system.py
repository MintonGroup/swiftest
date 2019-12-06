#satellites of Uranus
#inner sats from French, Dawson and Showalter (2015)
import numpy as np

#Units will be in terms of planet mass, planet radius, and years
G	    = 6.674e-8                       #Gravitational constant (cgs)
year    = 3600*24*365.25                #seconds in 1 year
AU      = 1.4960e13
M_Sun   = 1.9891e33


###***Define initial conditions***###
M_Uranus    = 8.6127e28
R_Uranus    = 25362e5
R_Uranus_eq = 25559e5
a_Uranus    = 19.2184 * AU
T_Uranus    = 0.71833 * 24 * 3600
W_Uranus    = 2.0 * np.pi / T_Uranus
Rhill_Uranus = a_Uranus * (M_Uranus / (3 * M_Sun))**(1.0 / 3.0)
Ipolar_Uranus = 0.230
J2_Uranus = 3341.29e-6
J4_Uranus = 33.61e-6


k2_Uranus = 0.104 #Gavrilov & Zharkov (1977)
Q_Uranus_TWlo = 11000.0
Q_Uranus_TWhi = 39000.0



a_Cord      = 4.98e9
M_Cord      = 3.88e19
R_Cord      = 2.01e6

a_Ophe      = 5.38e9
M_Ophe      = 5.10e19
R_Ophe      = 2.14e6

a_Bian      = 5.92e9
M_Bian      = 8.24e19
R_Bian      = 2.7e6

a_Cres      = 6.18e9
M_Cres      = 2.89e20
R_Cres      = 4.1e6

a_Desd      = 6.27e9
M_Desd      = 1.80e20
R_Desd      = 3.5e6

a_Juli      = 6.44e9
M_Juli      = 6.24e20
R_Juli      = 5.3e6

a_Port      = 6.61e9
M_Port      = 1.44e21
R_Port      = 7.0e6

a_Rosa      = 6.99e9
M_Rosa      = 1.95e20
R_Rosa      = 3.6e6

a_Cupi      = 7.44e9
M_Cupi      = 2.95e18
R_Cupi      = 8.9e5

a_Beli      = 7.53e9
M_Beli      = 3.82e20
R_Beli      = 4.5e6

a_Perd      = 7.64e9
M_Perd      = 9.85e18
R_Perd      = 1.33e6

a_Puck      = 8.6e9
M_Puck      = 2.23e21
R_Puck      = 8.1e6

a_Mab       = 9.77e9
M_Mab       = 7.99e17
R_Mab       = 1.24e6

a_Mira      = 1.30e10
M_Mira      = 6.59e22
R_Mira      = 2.35e7

a_Ariel     = 191020e5
M_Ariel     = 1353e21
R_Ariel     = 5.789e7

a_Umbriel   = 266300e5
M_Umbriel   = 1172e21
R_Umbriel   = 5.847e7

a_Titania   = 435910e5
M_Titania   = 3527e21
R_Titania   = 7.884e7

a_Oberon    = 583520e5
M_Oberon    = 3014e21
R_Oberon    = 7.614e7

eps_ring_width = 96.4e5
eps_ring_rad = 51149e5
eps_ring_r = np.linspace(eps_ring_rad - 0.5 * eps_ring_width, eps_ring_rad + 0.5 * eps_ring_width, num=100)
eps_ring_s = np.full_like(eps_ring_r,25.0)
eps_ring_s[0] = 0
eps_ring_s[-1] = 0


Sat_semimajor = np.array([a_Cord, a_Ophe, a_Bian, a_Cres, a_Desd, a_Juli, a_Port, a_Rosa, a_Cupi, a_Beli, a_Perd, a_Puck, a_Mab, a_Mira, a_Ariel, a_Umbriel, a_Titania, a_Oberon])
Sat_mass = np.array([M_Cord, M_Ophe, M_Bian, M_Cres, M_Desd, M_Juli, M_Port, M_Rosa, M_Cupi, M_Beli, M_Perd, M_Puck, M_Mab, M_Mira, M_Ariel, M_Umbriel, M_Titania, M_Oberon])
Sat_radius = np.array([R_Cord, R_Ophe, R_Bian, R_Cres, R_Desd, R_Juli, R_Port, R_Rosa, R_Cupi, R_Beli, R_Perd, R_Puck, R_Mab, R_Mira, R_Ariel, R_Umbriel, R_Titania, R_Oberon])
#Sat_name = np.array(["Cordelia", "Ophelia", "Bianca", "Cressida", "Desdemona", "Juliet", "Portia", "Rosalind", "Cupid", "Belinda", "Perdita", "Puck", "Mab", "Miranda", "Ariel", "Umbriel", "Titania", "Oberon"])
Sat_name = np.array(["Cor", "Oph", "Bia", "Cres", "Desd", "Juliet", "Por", "Ros", "Cu", "Bel", "Perd", "Puck", "Mab", "Miranda", "Ariel", "Umbriel", "Titania", "Oberon"])


