#!/usr/bin/env python3
import swiftest 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib.collections as clt
from scipy.spatial.transform import Rotation as R

xmin = -20.0
xmax = 20.0
ymin = -20.0
ymax = 20.0

#cases = ['supercat_head', 'supercat_off', 'disruption_head', 'disruption_off']
cases = ['disruption_off']

def scale_sim(ds, param):

    dsscale = ds

    GMtot = dsscale['GMass'].sum(skipna=True, dim="id").isel(time=0)
    rscale = sum(ds['Radius'].sel(id=[2, 3], time=0)).item()
    ds['Radius'] /= rscale

    dsscale['radmarker'] = dsscale['Radius'].fillna(0)

    dsscale['px'] /= rscale
    dsscale['py'] /= rscale
    dsscale['pz'] /= rscale

    mpx = dsscale['GMass'] * dsscale['px']
    mpy = dsscale['GMass'] * dsscale['py']
    mpz = dsscale['GMass'] * dsscale['pz']
    xbsys = mpx.sum(skipna=True, dim="id") / GMtot
    ybsys = mpy.sum(skipna=True, dim="id") / GMtot
    zbsys = mpz.sum(skipna=True, dim="id") / GMtot

    mvx = dsscale['GMass'] * dsscale['vx']
    mvy = dsscale['GMass'] * dsscale['vy']
    mvz = dsscale['GMass'] * dsscale['vz']
    vxbsys = mvx.sum(skipna=True, dim="id") / GMtot
    vybsys = mvy.sum(skipna=True, dim="id") / GMtot
    vzbsys = mvz.sum(skipna=True, dim="id") / GMtot

    dsscale['pxb'] = dsscale['px'] - xbsys
    dsscale['pyb'] = dsscale['py'] - ybsys
    dsscale['pzb'] = dsscale['pz'] - zbsys

    dsscale['vxb'] = dsscale['vx'] - vxbsys
    dsscale['vyb'] = dsscale['vy'] - vybsys
    dsscale['vzb'] = dsscale['vz'] - vzbsys

    return dsscale

class UpdatablePatchCollection(clt.PatchCollection):
    def __init__(self, patches, *args, **kwargs):
        self.patches = patches
        clt.PatchCollection.__init__(self, patches, *args, **kwargs)

    def get_paths(self):
        self.set_paths(self.patches)
        return self._paths

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""

    def __init__(self, ds, param):

        frame = 0
        nframes = ds['time'].size
        self.ds = scale_sim(ds, param)
        self.param = param
        self.rot_angle = {}

        self.clist = {'Initial conditions' : 'xkcd:windows blue',
                      'Disruption' : 'xkcd:baby poop',
                      'Supercatastrophic' : 'xkcd:shocking pink',
                      'Hit and run fragment' : 'xkcd:blue with a hint of purple',
                      'Central body'    : 'xkcd:almost black'}

        self.stream = self.data_stream(frame)
        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots(figsize=(8,8))
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1, frames=nframes,
                                          init_func=self.setup_plot, blit=False)
        self.ani.save(animfile, fps=60, dpi=300,
                      extra_args=['-vcodec', 'libx264'])

    def plot_pl_circles(self, pl, radmarker):
        patches = []
        for i in range(pl.shape[0]):
            s = plt.Circle((pl[i, 0], pl[i, 1]), radmarker[i])
            patches.append(s)
        return patches

    def vec_props(self, c):
        arrowprops = {
            'arrowstyle': '<|-',
            'mutation_scale': 20,
            'connectionstyle': 'arc3',
        }

        arrow_args = {
            'xycoords': 'data',
            'textcoords': 'data',
            'arrowprops': arrowprops,
            'annotation_clip': True,
            'zorder': 100,
            'animated' : True
        }
        aarg = arrow_args.copy()
        aprop = arrowprops.copy()
        aprop['color'] = c
        aarg['arrowprops'] = aprop
        aarg['color'] = c
        return aarg

    def plot_pl_vectors(self, pl, cval, r):
        varrowend, varrowtip = self.velocity_vectors(pl, r)
        arrows = []
        for i in range(pl.shape[0]):
            aarg = self.vec_props(cval[i])
            a = self.ax.annotate("",xy=varrowend[i],xytext=varrowtip[i], **aarg)
            arrows.append(a)
        return arrows

    def plot_pl_spins(self, pl, id, cval, len):
        sarrowend, sarrowtip = self.spin_arrows(pl, id, len)
        arrows = []
        for i in range(pl.shape[0]):
            aarg = self.vec_props(cval[i])
            aarg['arrowprops']['mutation_scale'] = 5
            aarg['arrowprops']['arrowstyle'] = "simple"
            a = self.ax.annotate("",xy=sarrowend[i],xytext=sarrowtip[i], **aarg)
            arrows.append(a)
        return arrows

    def origin_to_color(self, origin):
        cval = []
        for o in origin:
           c = self.clist[o]
           cval.append(c)

        return cval

    def velocity_vectors(self, pl, r):
        px = pl[:, 0]
        py = pl[:, 1]
        vx = pl[:, 2]
        vy = pl[:, 3]
        vmag = np.sqrt(vx ** 2 + vy ** 2)
        ux = np.zeros_like(vx)
        uy = np.zeros_like(vx)
        goodv = vmag > 0.0
        ux[goodv] = vx[goodv] / vmag[goodv]
        uy[goodv] = vy[goodv] / vmag[goodv]
        varrowend = []
        varrowtip = []
        for i in range(pl.shape[0]):
            vend = (px[i], py[i])
            vtip = (px[i] + vx[i] * self.v_length, py[i] + vy[i] * self.v_length)
            varrowend.append(vend)
            varrowtip.append(vtip)
        return varrowend, varrowtip

    def spin_arrows(self, pl, id, len):
        px = pl[:, 0]
        py = pl[:, 1]
        sarrowend = []
        sarrowtip = []
        for i in range(pl.shape[0]):
            endrel = np.array([0.0, len[i],  0.0])
            tiprel = np.array([0.0, -len[i], 0.0])
            r = R.from_rotvec(self.rot_angle[id[i]])
            endrel = r.apply(endrel)
            tiprel = r.apply(tiprel)
            send = (px[i] + endrel[0], py[i] + endrel[1])
            stip = (px[i] + tiprel[0], py[i] + tiprel[1])
            sarrowend.append(send)
            sarrowtip.append(stip)
        return sarrowend, sarrowtip

    def setup_plot(self):
        # First frame
        """Initial drawing of the scatter plot."""
        t, name, GMass, Radius, npl, pl, radmarker, origin = next(self.data_stream(0))

        cval = self.origin_to_color(origin)
        # set up the figure
        self.ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
        plt.axis('off')
        plt.tight_layout(pad=0)
        self.ax.set_aspect(1)
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)

        # Scale markers to the size of the system
        self.v_length = 0.50  # Length of arrow as fraction of velocity

        self.ax.margins(x=1, y=1)
        self.ax.set_xlabel('x distance / ($R_1 + R_2$)', fontsize='16', labelpad=1)
        self.ax.set_ylabel('y distance / ($R_1 + R_2$)', fontsize='16', labelpad=1)

        self.title = self.ax.text(0.50, 0.90, "", bbox={'facecolor': 'w', 'pad': 5}, transform=self.ax.transAxes,
                        ha="center", zorder=1000)

        self.title.set_text(titletext)
        self.patches = self.plot_pl_circles(pl, radmarker)

        self.collection = UpdatablePatchCollection(self.patches, color=cval, alpha=0.5, zorder=50)
        self.ax.add_collection(self.collection)
        #self.varrows = self.plot_pl_vectors(pl, cval, radmarker)
        self.sarrows = self.plot_pl_spins(pl, name, cval, radmarker)

        return self.collection, self.sarrows

    def update(self,frame):
        """Update the scatter plot."""
        t, name, GMass, Radius, npl, pl, radmarker, origin = next(self.data_stream(frame))
        cval = self.origin_to_color(origin)
        #varrowend, varrowtip = self.velocity_vectors(pl, radmarker)
        sarrowend, sarrowtip = self.spin_arrows(pl, name, radmarker)
        for i, p in enumerate(self.patches):
            p.set_center((pl[i, 0], pl[i,1]))
            p.set_radius(radmarker[i])
            p.set_color(cval[i])
            #self.varrows[i].set_position(varrowtip[i])
            #self.varrows[i].xy = varrowend[i]
            self.sarrows[i].set_position(sarrowtip[i])
            self.sarrows[i].xy = sarrowend[i]

        self.collection.set_paths(self.patches)
        return self.collection, self.sarrows

    def data_stream(self, frame=0):
        while True:
            d = self.ds.isel(time=frame)
            Radius = d['radmarker'].values
            GMass = d['GMass'].values
            x = d['pxb'].values
            y = d['pyb'].values
            vx = d['vxb'].values
            vy = d['vyb'].values
            name = d['id'].values
            npl = d.id.count().values
            id = d['id'].values
            rotx = d['rot_x'].values
            roty = d['rot_y'].values
            rotz = d['rot_z'].values

            radmarker = d['radmarker'].values
            origin = d['origin_type'].values

            t = self.ds.coords['time'].values[frame]
            self.mask = np.logical_not(np.isnan(x))

            x = np.nan_to_num(x, copy=False)
            y = np.nan_to_num(y, copy=False)
            vx = np.nan_to_num(vx, copy=False)
            vy = np.nan_to_num(vy, copy=False)
            radmarker = np.nan_to_num(radmarker, copy=False)
            GMass = np.nan_to_num(GMass, copy=False)
            Radius = np.nan_to_num(Radius, copy=False)
            rotx = np.nan_to_num(rotx, copy=False)
            roty = np.nan_to_num(roty, copy=False)
            rotz = np.nan_to_num(rotz, copy=False)
            rotvec = np.array([rotx, roty, rotz])
            self.rotvec = dict(zip(id, zip(*rotvec)))

            if frame == 0:
                tmp = np.zeros_like(rotvec)
                self.rot_angle = dict(zip(id, zip(*tmp)))
            else:
                t0 = self.ds.coords['time'].values[frame-1]
                dt = t - t0
                idxactive = np.arange(id.size)[self.mask]
                for i in id[idxactive]:
                    self.rot_angle[i] = self.rot_angle[i] + dt * np.array(self.rotvec[i])
            frame += 1
            yield t, name, GMass, Radius, npl, np.c_[x, y, vx, vy], radmarker, origin

for case in cases:
    if case == 'supercat_off':
        animfile = 'movies/supercat_off_axis.mp4'
        titletext = "Supercatastrophic - Off Axis"
        paramfile = 'param.supercatastrophic_off_axis.in'
    elif case == 'supercat_head':
        animfile = 'movies/supercat_headon.mp4'
        titletext = "Supercatastrophic - Head on"
        paramfile = 'param.supercatastrophic_headon.in'
    elif case == 'disruption_off':
        animfile = 'movies/disruption_off_axis.mp4'
        titletext = "Disruption - Off Axis"
        paramfile = 'param.disruption_off_axis.in'
    elif case == 'disruption_head':
        animfile = 'movies/disruption_headon.mp4'
        titletext = "Disruption- Head on"
        paramfile = 'param.disruption_headon.in'
    elif case == 'merger':
        animfile = 'movies/merger.mp4'
        titletext = "Merger"
        paramfile = 'param.merger.in'
    else:
        print(f'{case} is an unknown case')
        exit(-1)
    sim = swiftest.Simulation(param_file=paramfile)
    sim.bin2xr()
    ds = sim.ds
    print('Making animation')
    anim = AnimatedScatter(ds,sim.param)
    print('Animation finished')
    plt.close(fig='all')
