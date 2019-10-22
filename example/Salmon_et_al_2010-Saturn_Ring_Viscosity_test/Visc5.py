###***Model for viscous spreading.  Can calculate evolution of the surface mass density for either constant or variable viscosity.  Inputs
###***are x = switch for constant or variable (0 cont, 1 var), y is the timestep, dt the time between impacts, and z is the location of the interior bin
from __future__ import division
from math import *
from Init_Cond import *

nu = []
S = []
for a in range(int(N)):
    S.append(0.)
    nu.append(0.)

def f(x, M, R, t):        #X is a switch between constant viscosity and variable

    for a in range(int(N)):
        if x == 0:
            nu[a] = 1000.0*(r[a]/R)**2.0
        else:
            if sigma[a] <= 1.0e-6:
                nu[a] = 0.
                S[a] = sigma[a]*X[a]
                continue

            r_hstar = r[a]/(2.0*r_pdisk)*(2 * m_pdisk/(3.0*M))**(1.0/3.0)
            tau = pi*r_pdisk**2*sigma[a]/m_pdisk


            if r_hstar < 1.0:
                sigma_r = 0.5*(1-tanh((2.*r_hstar-1)/(r_hstar*(r_hstar-1))))*sqrt(GU*m_pdisk/r_pdisk) + 0.5*(1+tanh((2.*r_hstar-1)/(r_hstar*(r_hstar-1))))*(2.0*r_pdisk*w[a])
            else:
                sigma_r = sqrt(GU*m_pdisk/r_pdisk)


            Q = w[a]*sigma_r/(3.36*GU*(sigma[a]+0.000001))

            if Q <= 4.0:

                nu_trans_stable = sigma_r**2.0/(2.0*w[a])*(0.46*tau/(1.0+tau**2.0))
                nu_grav_stable = 0.0
                nu_trans_unstable = 13.0*r_hstar**5.0*GU**2.0*sigma[a]**2.0/w[a]**3.0
                nu_grav_unstable = nu_trans_unstable

                y = Q/4.
                nu_trans = 0.5*(1-tanh((2.*y-1)/(y*(y-1))))*(nu_trans_stable) +  0.5*(1+tanh((2.*y-1)/(y*(y-1))))*(nu_trans_unstable)
                nu_grav = 0.5*(1-tanh((2.*y-1)/(y*(y-1))))*(nu_grav_stable) +  0.5*(1+tanh((2.*y-1)/(y*(y-1))))*(nu_grav_unstable)

            else:
                nu_trans = sigma_r**2.0/(2.0*w[a])*(0.46*tau/(1.0+tau**2.0))
                nu_grav = 0.0

            nu_coll = r_pdisk**2.0*w[a]*tau
            nu[a] = nu_trans + nu_grav + nu_coll

