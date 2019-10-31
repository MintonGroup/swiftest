import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import Init_Cond as ic
from scipy.io import FortranFile

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
xmin = 1.0
xmax = 5.0
ymin = 1.0
ymax = 5e4

y2min = 1e13
y2max = 1e26
ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin,ymax))

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin, ymax)
ax.set_xlabel('Distance to Uranus (RU)')
ax.set_ylabel('$\Sigma$ (g$\cdot$cm$^{-2}$)')
ax.set_yscale('log')
secax = ax.twinx()
secax.set_yscale('log')
secax.set_ylabel('Mass of satellite (g)')
secax.set_ylim(y2min, y2max)
line, = ax.plot([], [], lw=2)
scatter, = secax.scatter([], [], marker='o', color="black", s=1, zorder = 50)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

def init2():
    scatter.set_data([], [])
    return scatter,


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


# animation function.  This is called sequentially
def animate(i):
    line.set_data(ring[i][0], ring[i][1])
    scatter.set_data(seeds[i][0], seeds[i][1])
    return line,scatter,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=20, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('sigma_animate.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()