#satellites of Uranus
#inner sats from French, Dawson and Showalter (2015)
import numpy as np

M_Uranus    = 8.6127e28
R_Uranus    = 25362e5

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


Sat_r_RM = []
Sat_M_Mass = []
area = []
Sat_colors = []


# plot existing satellite system
Sat_r_RM.append(a_Cord / R_Uranus)
Sat_M_Mass.append(M_Cord)
area.append(np.pi * (75. * (10. * R_Cord / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Ophe / R_Uranus)
Sat_M_Mass.append(M_Ophe)
area.append(np.pi * (75. * (10. * R_Ophe / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Bian / R_Uranus)
Sat_M_Mass.append(M_Bian)
area.append(np.pi * (75. * (10. * R_Bian / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Cres / R_Uranus)
Sat_M_Mass.append(M_Cres)
area.append(np.pi * (75. * (10. * R_Cres / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Desd / R_Uranus)
Sat_M_Mass.append(M_Desd)
area.append(np.pi * (75. * (10. * R_Desd / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Juli / R_Uranus)
Sat_M_Mass.append(M_Juli)
area.append(np.pi * (75. * (10. * R_Juli / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Port / R_Uranus)
Sat_M_Mass.append(M_Port)
area.append(np.pi * (75. * (10. * R_Port / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Rosa / R_Uranus)
Sat_M_Mass.append(M_Rosa)
area.append(np.pi * (75. * (10. * R_Rosa / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Cupi / R_Uranus)
Sat_M_Mass.append(M_Cupi)
area.append(np.pi * (75. * (10. * R_Cupi / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Beli / R_Uranus)
Sat_M_Mass.append(M_Beli)
area.append(np.pi * (75. * (10. * R_Beli / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Perd / R_Uranus)
Sat_M_Mass.append(M_Perd)
area.append(np.pi * (75. * (10. * R_Perd / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Puck / R_Uranus)
Sat_M_Mass.append(M_Puck)
area.append(np.pi * (75. * (10. * R_Puck / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Mab / R_Uranus)
Sat_M_Mass.append(M_Mab)
area.append(np.pi * (75. * (10. * R_Mab / R_Uranus)) ** 2)
Sat_colors.append('k')

Sat_r_RM.append(a_Mira / R_Uranus)
Sat_M_Mass.append(M_Mira)
area.append(np.pi * (75. * (10. * R_Mira / R_Uranus)) ** 2)
Sat_colors.append('k')