import numpy as np
import pyvista as pv
import xarray as xr

import swiftest

num_movie_frames = 300


def encounter_combiner(sim):
    """
    Combines simulation data with encounter data to produce a dataset that contains the position, mass, radius, etc of both.

    It will interpolate over empty time values to fill in gaps.
    """
    names = ["Target", "Projectile"]
    # Only keep a minimal subset of necessary data from the simulation and encounter datasets
    keep_vars = ["name", "rh", "vh", "Gmass", "radius", "rot"]
    keep_names = sim.data["name"].values[1:]
    data = sim.data[keep_vars].sel(name=keep_names)
    enc = sim.encounters[keep_vars].sel(name=keep_names).load()

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
    smooth_time = np.linspace(start=ds.time.values[0], stop=ds.time.values[-1], num=int(1.2 * num_movie_frames))
    ds = ds.interp(time=smooth_time)
    ds["rotangle"] = xr.zeros_like(ds["rot"])
    ds["rot"] = ds["rot"].fillna(0.0)

    return ds


sim = swiftest.Simulation(simdir="disruption_off_axis", read_data=True)
ds = encounter_combiner(sim)
plotter = pv.Plotter()
params = {"t_index": 0}


def barycenter(ds):
    rh = ds.rh.where(~np.isnan(ds.Gmass), drop=True)
    vh = ds.vh.where(~np.isnan(ds.Gmass), drop=True)
    Gm = ds.Gmass.where(~np.isnan(ds.Gmass), drop=True)
    return rh.weighted(Gm).mean(dim="name").values, vh.weighted(Gm).mean(dim="name").values


def plot_frame(t_index):
    plotter.clear()
    plotter.enable_lightkit()
    ids = ds.isel(time=t_index)
    for n in ids.name:
        if n == "":
            continue
        rb, vb = barycenter(ids)
        center = ids.sel(name=n).rh.values - rb
        radius = ids.sel(name=n).radius.values
        vh = ids.sel(name=n).vh - vb
        vmag = ids.sel(name=n).vh.magnitude().item()
        vh = vh.values
        if np.isnan(radius):
            continue
        sphere = pv.Sphere(radius=radius, center=center)
        arrow = pv.Arrow(start=center, direction=vh, scale=vmag * 1e-6)
        plotter.add_mesh(sphere, name=n.item())
        plotter.add_mesh(arrow, name=n.item() + "_vel")


def update_time(t_index):
    plot_frame(t_index)
    plotter.update()


def update_time_plus():
    if params["t_index"] < len(sim.encounters.time) - 1:
        params["t_index"] += 1
    update_time(params["t_index"])


def update_time_minus():
    if params["t_index"] > 0:
        params["t_index"] -= 1
    update_time(params["t_index"])


plot_frame(0)
plotter.add_key_event("n", update_time_plus)
plotter.add_key_event("b", update_time_minus)

plotter.show()
