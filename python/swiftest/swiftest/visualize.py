"""
 Copyright 2023 - David Minton,
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest. 
 If not, see: https://www.gnu.org/licenses. 
"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

def _square_plot():
    figsize = (4,4)
    fig = plt.figure(figsize=figsize, dpi=300)
    plt.tight_layout(pad=0)

    ax = plt.Axes(fig, [0.1, 0.1, 0.8, 0.8])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    fig.add_axes(ax)
    return fig, ax

def select_one_collision(collisions, collision_id):
    ds = collisions.sel(collision_id=collision_id)
    bgood = ds.sel(stage='before')['name'].where(ds.sel(stage='before')['particle_type'] != 'nan',drop=True)
    agood = ds.sel(stage='after')['name'].where(ds.sel(stage='after')['particle_type'] != 'nan',drop=True)
    goodname=np.unique(np.concatenate((bgood,agood)))
    ds = ds.sel(name=goodname)
    
    return ds
    
def collisions(collisions, collision_id):
    fig, ax = _square_plot()
    
    ds = select_one_collision(collisions, collision_id)
    
    
    return