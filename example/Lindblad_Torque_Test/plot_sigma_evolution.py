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

y2min = 1e13
y2max = 1e23
secaxes = {}

for key in axes:
    axes[key].set_xlim(xmin,xmax)
    axes[key].set_ylim(ymin, ymax)
    axes[key].set_xlabel('Distance to Uranus (RU)')
    axes[key].set_ylabel('$\Sigma$ (g$\cdot$cm$^{-2}$)')
    axes[key].set_yscale('log')
    secaxes[key] = axes[key].twinx()
    secaxes[key].set_yscale('log')
    secaxes[key].set_ylabel('Mass of satellite (g)')
    secaxes[key].set_ylim(y2min, y2max)

axes['a'].title.set_text('$0$ My')
axes['b'].title.set_text('$17$ My')
axes['c'].title.set_text('$183$ Mys')
axes['d'].title.set_text('$545$ years')
axes['e'].title.set_text('$546$ years')
axes['f'].title.set_text('$720$ years')
ring = {}
seeds = {}

with FortranFile('ring.dat', 'r') as f:
    while True:
        try:
            t = f.read_reals(np.float64)
        except:
            break
        Nbin = f.read_ints(np.int32)
        r = f.read_reals(np.float64)
        Gsigma = f.read_reals(np.float64)
        nu = f.read_reals(np.float64)
        kval = int(t / ic.t_print)
        ring[kval] = [r, Gsigma, nu]
        Nseeds = f.read_ints(np.int32)
        a = f.read_reals(np.float64)
        Gm = f.read_reals(np.float64)
        seeds[kval] = [a, Gm]


#convert the units
for key in ring:
    ring[key][0] /= ic.RP  #convert radius to planet radius
    ring[key][1] *= ic.MU2GM / ic.DU2CM**2 / ic.GU  # convert surface mass density to cgs
    ring[key][2] *= ic.DU2CM**2 / ic.TU2S # convert viscosity to cgs
    seeds[key][0] /= ic.RP
    seeds[key][1] *= ic.MU2GM / ic.GU

# These are the output times to plot

tout = np.array([0.0, 17e6, 183e6, 545e6, 546e6, 720e6]) * ic.year / ic.TU2S
nt = np.rint(tout / ic.t_print).astype(int)

tn = {  'a' : nt[0],
        'b' : nt[1],
        'c' : nt[2],
        'd' : nt[3],
        'e' : nt[4],
        'f' : nt[5]}

for key in axes:
    axes[key].plot(ring[tn[key]][0], ring[tn[key]][1], '-', color="black", linewidth=1.0, zorder = 50, label = "SyMBA-RINGMOONS")
    secaxes[key].scatter(seeds[tn[key]][0], seeds[tn[key]][1], marker='o', color="black", s=1)

#axes['b'].plot(ring[nt[1]][0], ring[nt[1]][1], '-', color="black", linewidth=1.0, zorder = 50)
#axes['c'].plot(ring[nt[2]][0], ring[nt[2]][1], '-', color="black", linewidth=1.0, zorder = 50)
#axes['d'].plot(ring[nt[3]][0], ring[nt[3]][1], '-', color="black", linewidth=1.0, zorder = 50)
#axes['e'].plot(ring[nt[4]][0], ring[nt[4]][1], '-', color="black", linewidth=1.0, zorder = 50)
#axes['f'].plot(ring[nt[5]][0], ring[nt[5]][1], '-', color="black", linewidth=1.0, zorder = 50)
axes['a'].legend(loc='upper left',prop={'size': 8})
figure.tight_layout()
#plt.show()

figname = "Uranus_ring_satellite_evoloution.png"
plt.savefig(figname,dpi=300 )
os.system(f'open {figname}')


