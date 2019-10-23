import numpy as np
import matplotlib.pyplot as plt
import os
import Init_Cond as ic
from scipy.io import FortranFile

figure = plt.figure(1, figsize=(8,6))
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
axes['a'].plot(S2010test_t0[:,0]*1e-3, S2010test_t0[:,1]*1e-4, '-', color="red", linewidth=1.0, zorder=40, label = "Salmon et al. (2010)")
axes['b'].plot(S2010test_t1e3[:,0]*1e-3, S2010test_t1e3[:,1]*1e-4, '-', color="red", linewidth=1.0, zorder=40)
axes['c'].plot(S2010test_t1e4[:,0]*1e-3, S2010test_t1e4[:,1]*1e-4, '-', color="red", linewidth=1.0, zorder=40)
axes['d'].plot(S2010test_t1e5[:,0]*1e-3, S2010test_t1e5[:,1]*1e-4, '-', color="red", linewidth=1.0, zorder=40)

axes['a'].title.set_text('$0$ years')
axes['b'].title.set_text('$10^3$ years')
axes['c'].title.set_text('$10^4$ years')
axes['d'].title.set_text('$10^5$ years')
ring = {}

with FortranFile('ring.dat', 'r') as f:
    while True:
        try:
            t = f.read_reals(np.float64)
        except:
            break
        N = f.read_ints(np.int32)
        r = f.read_reals(np.float64)
        Gsigma = f.read_reals(np.float64)
        nu = f.read_reals(np.float64)
        kval = int(t / ic.t_print)
        ring[f'{kval}'] = [r, Gsigma, nu]

#convert the units
for key in ring:
    ring[key][0] *= ic.DU2CM * 1e-8 #convert radius to 1000 km
    ring[key][1] *= ic.MU2GM / ic.DU2CM**2 / ic.GU * 1e-3  # convert surface mass density to 1e4 kg/m^2
    ring[key][2] *= ic.DU2CM**2 / ic.TU2S * 1e-2 # convert viscosity to m^2/s

# These are the output times to plot
tout = np.array([0.0, 1e3, 1e4, 1e5]) * ic.year / ic.TU2S
nt = np.rint(tout / ic.t_print).astype(int)

axes['a'].plot(ring[f'{nt[0]}'][0], ring[f'{nt[0]}'][1], '-', color="black", linewidth=1.0, zorder = 50, label = "SyMBA-RINGMOONS")
axes['b'].plot(ring[f'{nt[1]}'][0], ring[f'{nt[1]}'][1], '-', color="black", linewidth=1.0, zorder = 50)
axes['c'].plot(ring[f'{nt[2]}'][0], ring[f'{nt[2]}'][1], '-', color="black", linewidth=1.0, zorder = 50)
axes['d'].plot(ring[f'{nt[3]}'][0], ring[f'{nt[3]}'][1], '-', color="black", linewidth=1.0, zorder = 50)
axes['a'].legend(loc='upper left',prop={'size': 8})
figure.tight_layout()
#plt.show()

figname ="Salmon_et_al_2010-Saturn_Ring_Viscocity_evoloution.png"
plt.savefig(figname,dpi=300 )
os.system(f'open {figname}')


