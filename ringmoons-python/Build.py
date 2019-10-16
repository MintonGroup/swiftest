import numpy as np
###***My modules***###
from Init_Cond import *
from Lists import *

def f(x, M):
    Roche = 2.456 * RP * (rhoP / rho_sat)**(1./3.) #density of satellites
    T = 4.0 * np.pi * MP * MP**2.0 / (5.0 * LP)
    r_sync = (T * np.sqrt(G * MP) / (2.0 * np.pi))**(2.0 / 3.0)


    #Creates initial values for the disk and prints them out
    for a in range(int(N)):

        r.append(r_I + deltar*(a+0.5))
        X.append(2.*r[a]**0.5)
        R_in_RP.append(r[a]/MP)
        deltaA.append(2*np.pi*deltar*(r_I + deltar*(a + 0.5)))

        #Power law surface mass density profile
        if x == 0:
            m.append(2*np.pi*alpha/(2-p)*((r[a] + deltar/2)**(2-p) - (r[a] - deltar/2)**(2-p)))
            sigma.append(m[a]/deltaA[a])
            #build the threshold surface mass density for the system to push a satellite beyond synchronous
            if 0.5**(2./3.)*Roche < r[a] < 0.5**(2./3.)*r_sync: #bounds are set up from lindblad resonance locations at FRL and synch
                sigma_threshold.append(3.*k_2*MP**5.*M/(2.**7.*Q*r[a]**7.))   #see google doc May 15th, 2017
            else:
                sigma_threshold.append(0.)  #outside of relevant region of disk, doesn't matter

        #Gaussian surface mass density profile
        elif x == 1:
            centroid = 20       #bin id of the center of the gaussian
            spread = 1          #width of the gaussian
            mass_scale = 2.3e6  #scale factor to get a given mass
            sigma.append(mass_scale/(5*(2*np.pi)**(0.5))*exp(-(a-centroid)**2/(2*spread**2)))
            m.append(sigma[a]*deltaA[a])

            #build the threshold surface mass density for the system to push a satellite beyond synchronous
            if 0.5**(2./3.)*Roche < r[a] < 0.5**(2./3.)*r_sync: #bounds are set up from lindblad resonance locations at FRL and synch
                sigma_threshold.append(3.*k_2*MP**5.*M/(2.**7.*Q*r[a]**7.))   #see google doc May 15th, 2017
            else:
                sigma_threshold.append(0.)  #outside of relevant region of disk, doesn't matter

        else:
            print('You have not chosen a valid disk model')
        R.append(r[a]**2 + deltar**2/4)
        I.append(m[a]*R[a])
        w.append((G*M/r[a]**3)**0.5)
        Torque_to_disk.append(0.0)
