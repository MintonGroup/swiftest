import numpy as np
import matplotlib.pyplot as plt
import os


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

S2010test_t0 = np.loadtxt("S2010test-t0.csv", skiprows=1)
axes['a'].plot(S2010test_t0[:,0]*1e-3, S2010test_t0[:,1]*1e-4, '-', color="black", linewidth=2.0, zorder=40)

plt.show()



