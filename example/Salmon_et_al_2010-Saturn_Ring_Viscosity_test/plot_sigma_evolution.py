import numpy as np
import matplotlib.pyplot as plt
import os
import Init_Cond as ic
import Visc5


figure = plt.figure(1, figsize=(8,6))    #this should be called later when we know the aspect ratio, but for now, I have it here
axes = {'a' : figure.add_subplot(221),
        'b' : figure.add_subplot(222),
        'c' : figure.add_subplot(223),
        'd' : figure.add_subplot(224)}
xmin = 100
xmax = 120
ymin = 0
ymax = 8

for key in axes:
    axes[key].set_xlim(xmin,xmax)
    axes[key].set_ylim(ymin, ymax)
    axes[key].set_xlabel('Distance to Saturn (1000 km)')
    axes[key].set_ylabel('$\Sigma$ ($10^4$ kg$\cdot$m$^{-2}$)')

S2010test_t0 = np.loadtxt('S2010test-t0.csv',delimiter=',')
S2010test_t1e3 = np.loadtxt('S2010test-t1e3.csv',delimiter=',')
S2010test_t1e4 = np.loadtxt('S2010test-t1e4.csv',delimiter=',')
S2010test_t1e5 = np.loadtxt('S2010test-t1e5.csv',delimiter=',')
axes['a'].plot(S2010test_t0[:,0]*1e-3, S2010test_t0[:,1]*1e-4, '-', color="red", linewidth=1.0, zorder=40)
axes['b'].plot(S2010test_t1e3[:,0]*1e-3, S2010test_t1e3[:,1]*1e-4, '-', color="red", linewidth=1.0, zorder=40)
axes['c'].plot(S2010test_t1e4[:,0]*1e-3, S2010test_t1e4[:,1]*1e-4, '-', color="red", linewidth=1.0, zorder=40)
axes['d'].plot(S2010test_t1e5[:,0]*1e-3, S2010test_t1e5[:,1]*1e-4, '-', color="red", linewidth=1.0, zorder=40)

axes['a'].title.set_text('$0$ years')
axes['b'].title.set_text('$10^3$ years')
axes['c'].title.set_text('$10^4$ years')
axes['d'].title.set_text('$10^5$ years')

r = np.asarray(ic.r) * 1e-2 * 1e-3 * 1e-3 #cm to 1000 km
sigma = np.asarray(ic.sigma) * 1e-3 * 1e4 * 1e-4 #g/cm**2 to 1e4 kg/m**2
axes['a'].plot(r, sigma, '-', color="black", linewidth=1.0, zorder = 50)

#ringfile = 'ring.in'
#with open(ringfile) as f:
   #N = int(f.readline())
  # sigma = np.empty(N,dtype=np.float64)
  # deltar=f.readline()
  # vals = [float(r) for r in f.readline().split()]
  # r_pdisk = vals[0]
  # Gm_pdisk = vals[1]
  # for i, line in enumerate(f):
  #     sigma[i] = line
#
#print(Gm_pdisk)

Visc5.f(1,ic.M_Saturn,t=0.0)
print(f'nu = {Visc5.nu[501]*1e-4} m^2 s^-1')

figure.tight_layout()
#plt.show()



#figname ="Salmon_et_al_2010-Saturn_Ring_Viscocity_evoloution.png"
#plt.savefig(figname,dpi=300 )
#os.system(f'open {figname}')


