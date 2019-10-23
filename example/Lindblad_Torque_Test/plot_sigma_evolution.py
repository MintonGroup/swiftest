import numpy as np
import matplotlib.pyplot as plt
import os
import Init_Cond as ic
from scipy.io import FortranFile

figure = plt.figure(1, figsize=(8,9))
axes = {'a' : figure.add_subplot(321),
        'b' : figure.add_subplot(322),
        'c' : figure.add_subplot(323),
        'd' : figure.add_subplot(324),
        'e' : figure.add_subplot(325),
        'f' : figure.add_subplot(326)}
xmin = 1.0
xmax = 5.0
ymin = 1.0
ymax = 5e4

for key in axes:
    axes[key].set_xlim(xmin,xmax)
    axes[key].set_ylim(ymin, ymax)
    axes[key].set_xlabel('Distance to Uranus (RU)')
    axes[key].set_ylabel('$\Sigma$ (g$\cdot$cm$^{-2}$)')
    axes[key].set_yscale('log')

#S2010test_t0 = np.loadtxt('S2010test-t0.csv',delimiter=',')
#S2010test_t1e3 = np.loadtxt('S2010test-t1e3.csv',delimiter=',')
#S2010test_t1e4 = np.loadtxt('S2010test-t1e4.csv',delimiter=',')
#S2010test_t1e5 = np.loadtxt('S2010test-t1e5.csv',delimiter=',')
#axes['a'].plot(S2010test_t0[:,0]*1e-3, S2010test_t0[:,1]*1e-4, '-', color="red", linewidth=1.0, zorder=40, label = "Salmon et al. (2010)")
#axes['b'].plot(S2010test_t1e3[:,0]*1e-3, S2010test_t1e3[:,1]*1e-4, '-', color="red", linewidth=1.0, zorder=40)
#axes['c'].plot(S2010test_t1e4[:,0]*1e-3, S2010test_t1e4[:,1]*1e-4, '-', color="red", linewidth=1.0, zorder=40)
#axes['d'].plot(S2010test_t1e5[:,0]*1e-3, S2010test_t1e5[:,1]*1e-4, '-', color="red", linewidth=1.0, zorder=40)

axes['a'].title.set_text('$0$ My')
axes['b'].title.set_text('$17$ My')
axes['c'].title.set_text('$183$ Mys')
axes['d'].title.set_text('$545$ years')
axes['e'].title.set_text('$546$ years')
axes['f'].title.set_text('$720$ years')
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
    ring[key][0] /= ic.RP  #convert radius to planet radius
    ring[key][1] *= ic.MU2GM / ic.DU2CM**2 / ic.GU  # convert surface mass density to cgs
    ring[key][2] *= ic.DU2CM**2 / ic.TU2S # convert viscosity to cgs

# These are the output times to plot
tout = np.array([0.0, 1e4, 1e5, 1e7]) * ic.year / ic.TU2S
nt = np.rint(tout / ic.t_print).astype(int)

axes['a'].plot(ring[f'{nt[0]}'][0], ring[f'{nt[0]}'][1], '-', color="black", linewidth=1.0, zorder = 50, label = "SyMBA-RINGMOONS")
axes['b'].plot(ring[f'{nt[1]}'][0], ring[f'{nt[1]}'][1], '-', color="black", linewidth=1.0, zorder = 50)
axes['c'].plot(ring[f'{nt[2]}'][0], ring[f'{nt[2]}'][1], '-', color="black", linewidth=1.0, zorder = 50)
axes['d'].plot(ring[f'{nt[3]}'][0], ring[f'{nt[3]}'][1], '-', color="black", linewidth=1.0, zorder = 50)
#axes['a'].legend(loc='upper left',prop={'size': 8})
figure.tight_layout()
#plt.show()

figname ="Uranus_ring_satellite_evoloution.png"
plt.savefig(figname,dpi=300 )
os.system(f'open {figname}')


