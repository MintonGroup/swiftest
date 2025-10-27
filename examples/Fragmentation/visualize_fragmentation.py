import numpy as np
import pyvista as pv
import xarray as xr

import swiftest


class SimulationVisualizer:
    """
    Visualizes a fragmentation simulation using PyVista.

    Parameters
    ----------
    simdir : str
        The directory containing the simulation data.
    nframes : int, optional
        The number of frames to visualize, by default 1000.
    """

    def __init__(self, simdir, nframes=1000):
        self.sim = swiftest.Simulation(simdir=simdir, read_data=True)

        self.nframes = nframes
        self.ds = self.encounter_combiner()
        # Set the initial frame to just after the collision
        tcollision = self.sim.collisions.sel(collision_id=1).time.item()
        tframe = self.ds.sel(time=tcollision, method="nearest").time.item()
        self.iframe = np.where(self.ds.time.values == tframe)[0].item()
        if tframe < tcollision:
            self.iframe += 1

        self.plotter = pv.Plotter()
        self.plotter.clear()
        self.plotter.enable_lightkit()
        self.sphere_actors = {}
        self.arrow_actors = {}
        self.arrow_vectors = {}
        self.arrow_base_points = {}
        self.arrow_meshes = {}
        self.rb, self.vb = self.barycenter(self.ds.isel(time=self.iframe))
        self.plotter_is_active = False

        self.update_time()
        self.plotter.add_key_event("n", self.update_time_plus)
        self.plotter.add_key_event("b", self.update_time_minus)
        self.plotter_is_active = True
        self.plotter.show()

    def get_state(self, name):
        """
        Get the state of a body at the current frame.

        Parameters
        ----------
        name : str
            The name of the body.

        Returns
        -------
        radius : float
            The radius of the body.
        rcenter : np.ndarray
            The position of the body relative to the barycenter.
        vcenter : np.ndarray
            The velocity of the body relative to the barycenter.
        vmag : float
            The magnitude of the velocity of the body.
        """
        ds = self.ds.isel(time=self.iframe).sel(name=name)
        rcenter = ds.rh.values - self.rb
        radius = ds.radius.values
        vcenter = ds.vh.values - self.vb
        vmag = ds.vh.magnitude().item()
        return radius, rcenter, vcenter, vmag

    def add_actors(self, name):
        """
        Add the sphere and arrow actors for a body.

        Parameters
        ----------
        name : str
            The name of the body.
        """
        radius, rcenter, vcenter, vmag = self.get_state(name)
        if np.isnan(vmag):
            return
        sphere = pv.Sphere(radius=radius, center=rcenter)
        self.plotter.add_mesh(sphere, name=name + "_body")
        self.sphere_actors[name] = sphere

        # Canonical arrow geometry: unit arrow along +X at origin
        base_arrow = pv.Arrow(
            start=(0.0, 0.0, 0.0),
            direction=(1.0, 0.0, 0.0),
            scale=1.0,
            shaft_radius=0.01,
            tip_radius=0.02,
            tip_length=0.05,
            tip_resolution=10,
        )

        # Make a working copy that will be shown and updated in-place
        arrow = base_arrow.copy(deep=True)

        # Build a 4x4 transform that maps +X to vdir, scales by vmag*1e-6, and translates to rcenter
        def _rodrigues_from_x(to_dir: np.ndarray) -> np.ndarray:
            x = np.array([1.0, 0.0, 0.0])
            v = to_dir / np.linalg.norm(to_dir)
            c = float(np.dot(x, v))
            if c > 0.999999:
                return np.eye(3)
            if c < -0.999999:
                # 180° around any axis perpendicular to X; choose Z unless v is ±X, then use Y
                k = np.array([0.0, 0.0, 1.0])
                if abs(v[2]) > 0.9:
                    k = np.array([0.0, 1.0, 0.0])
            else:
                k = np.cross(x, v)
                k /= np.linalg.norm(k)
            K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]], dtype=float)
            # angle via acos, but use Rodrigues with sin/cos from cross/dot for stability
            s = np.linalg.norm(np.cross(x, v))
            R = np.eye(3) + K * s + (K @ K) * (1.0 - c)
            return R

        vdir = vcenter / vmag
        R = _rodrigues_from_x(vdir)
        s = vmag * 1e-6
        M = np.eye(4)
        M[:3, :3] = R * s
        M[:3, 3] = rcenter
        arrow.transform(M, inplace=True)

        actor = self.plotter.add_mesh(arrow, name=name + "_vel")
        self.arrow_actors[name] = actor  # store the pyvista Actor for fast transforms if needed
        self.arrow_meshes[name] = arrow  # store the PolyData being displayed
        self.arrow_base_points[name] = base_arrow.points.copy()  # canonical points for rebuild

        return

    def update_actors(self, name):
        """
        Update the sphere and arrow actors for a body.

        Parameters
        ----------
        name : str
            The name of the body.
        """
        radius, rcenter, vcenter, vmag = self.get_state(name)
        if np.isnan(vmag):
            if name in self.sphere_actors:
                self.plotter.remove_actor(name + "_body")
                del self.sphere_actors[name]
            if name in self.arrow_actors:
                self.plotter.remove_actor(name + "_vel")
                del self.arrow_actors[name]
                del self.arrow_meshes[name]
                del self.arrow_base_points[name]
            return
        if name in self.sphere_actors and name in self.arrow_actors:
            delta_r = rcenter - self.sphere_actors[name].center
            self.sphere_actors[name].translate(delta_r, inplace=True)

            # Update the existing arrow geometry in-place using a single transform of the canonical points
            vdir = vcenter / vmag

            def _rodrigues_from_x(to_dir: np.ndarray) -> np.ndarray:
                x = np.array([1.0, 0.0, 0.0])
                v = to_dir / np.linalg.norm(to_dir)
                c = float(np.dot(x, v))
                if c > 0.999999:
                    return np.eye(3)
                if c < -0.999999:
                    k = np.array([0.0, 0.0, 1.0])
                    if abs(v[2]) > 0.9:
                        k = np.array([0.0, 1.0, 0.0])
                else:
                    k = np.cross(x, v)
                    k /= np.linalg.norm(k)
                K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]], dtype=float)
                s = np.linalg.norm(np.cross(x, v))
                R = np.eye(3) + K * s + (K @ K) * (1.0 - c)
                return R

            R = _rodrigues_from_x(vdir)
            s = vmag * 1e-6
            # Build 4x4 matrix once and apply directly to points
            M = np.eye(4)
            M[:3, :3] = R * s
            M[:3, 3] = rcenter

            base_pts = self.arrow_base_points[name]
            ones = np.ones((base_pts.shape[0], 1))
            hom = np.hstack([base_pts, ones])
            new_pts = (hom @ M.T)[:, :3]
            self.arrow_meshes[name].points = new_pts
        else:
            self.add_actors(name)
        return

    def encounter_combiner(self):
        """
        Combines simulation data with encounter data to produce a dataset that contains the position, mass, radius, etc of both.

        It will interpolate over empty time values to fill in gaps.
        """
        names = ["Target", "Projectile"]
        # Only keep a minimal subset of necessary data from the simulation and encounter datasets
        keep_vars = ["name", "rh", "vh", "Gmass", "radius", "rot"]
        keep_names = self.sim.data["name"].values[1:]
        data = self.sim.data[keep_vars].sel(name=keep_names)
        enc = self.sim.encounters[keep_vars].sel(name=keep_names).load()

        # Remove any encounter data at the same time steps that appear in the data to prevent duplicates
        t_not_duplicate = ~data["time"].isin(enc["time"]).values
        data = data.sel(time=t_not_duplicate)
        tgood = enc.time.where(~np.isnan(enc.time), drop=True)
        enc = enc.sel(time=tgood)

        # The following will combine the two datasets along the time dimension, sort the time dimension, and then fill in any time gaps with interpolation
        ds = xr.combine_nested([data, enc], concat_dim="time").sortby("time").interpolate_na(dim="time")

        # Rename the merged Target body so that their data can be combined
        tname = [n for n in ds["name"].data if names[0] in n]
        nottname = [n for n in ds["name"].data if names[0] not in n]
        dslist = []
        for n in tname:
            dsnew = ds.sel(name=n)
            dsnew["name"] = names[0]
            dslist.append(dsnew)

        newds = xr.merge(dslist, compat="no_conflicts")
        ds = xr.combine_nested([ds.sel(name=nottname), newds], concat_dim="name")

        # Interpolate in time to make a smooth, constant time step dataset
        # Add a bit of padding to the time, otherwise there are some issues with the interpolation in the last few frames.
        smooth_time = np.linspace(start=ds.time.values[0], stop=ds.time.values[-1], num=int(1.2 * self.nframes))
        ds = ds.interp(time=smooth_time)
        ds["rotangle"] = xr.zeros_like(ds["rot"])
        ds["rot"] = ds["rot"].fillna(0.0)

        return ds

    @staticmethod
    def barycenter(ds):
        """
        Compute the barycenter position and velocity of the system.

        Parameters
        ----------
        ds : xarray.Dataset
            The dataset containing the simulation data at a given time.

        Returns
        -------
        rh_barycenter : np.ndarray
            The barycenter position.
        vh_barycenter : np.ndarray
            The barycenter velocity.
        """
        rh = ds.rh.where(~np.isnan(ds.Gmass), drop=True)
        vh = ds.vh.where(~np.isnan(ds.Gmass), drop=True)
        Gm = ds.Gmass.where(~np.isnan(ds.Gmass), drop=True)
        return rh.weighted(Gm).mean(dim="name").values, vh.weighted(Gm).mean(dim="name").values

    def update_time(self):
        """
        Update the visualization to the current frame.
        """
        self.rb, self.vb = self.barycenter(self.ds.isel(time=self.iframe))
        for name in self.ds.name.data:
            if name == "":
                continue
            self.update_actors(name)
        if self.plotter_is_active:
            self.plotter.update()
        return

    def update_time_plus(self):
        """
        Advance the visualization to the next frame.
        """
        if self.iframe < self.nframes - 1:
            self.iframe += 1
        self.update_time()

    def update_time_minus(self):
        """
        Go back to the previous frame in the visualization.
        """
        if self.iframe > 0:
            self.iframe -= 1
        self.update_time()


vis = SimulationVisualizer(simdir="disruption_headon")
