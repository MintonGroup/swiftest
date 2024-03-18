"""
Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
This file is part of Swiftest.
Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Swiftest. 
If not, see: https://www.gnu.org/licenses. 
"""

import numpy as np
import xarray as xr
from scipy.spatial.transform import Rotation as R

def magnitude(da: xr.DataArray) -> xr.DataArray:
    """
    Computes the magnitude of a vector quantity from a Dataset.
    
    Parameters
    ----------
    da : xarray DataArray
        DataArray containing the vector quantity
        
    Returns
    -------
    mag : Xarray DataArray
        DataArray containing the magnitude of the vector quantity 
    """
    dim = "space"
    ord = None
    
    if dim not in da.dims:
        raise ValueError(f"Dimension {dim} not found in DataArray")
    return xr.apply_ufunc(
        np.linalg.norm, da.where(~np.isnan(da)), input_core_dims=[[dim]], kwargs={"ord": ord, "axis": -1}, dask="allowed"
    )
    
def wrap_angle(angle):
    """
    Converts angles to be between 0 and 360 degrees.
        
    Parameters
    ----------
    angle : float
        Angle to be converted
    Returns
    -------
    angle : float
        Converted angle
    """
    while np.any(angle >= 360.0 ):
        angle[angle >= 360.0] -= 360.0
    while np.any(angle < 0.0):
        angle[angle < 0.0] += 360.0
    return angle


def follow_swift(ds, ifol=None, nskp=None):
    """
    Emulates the Swift version of tool_follow.f
    It will generate a file called follow.out containing 10 columns (angles are all in degrees)::
    
        1 2 3  4    5     6    7    8    9   10
        t,a,e,inc,capom,omega,capm,peri,apo,obar
        
    Parameters
    ----------
    ds : Xarray Dataset 
        Dataset containing orbital elements
    ifol : int, optional
        Particle number to follow. The default is None, in which case the user is prompted to enter a particle number.
    nskp : int, optional
        Print frequency. The default is None, in which case the user is prompted to enter a print frequency. 

    Returns
    -------
    fol : Xarray Dataset
        Dataset containing only the followed body with angles converted to degrees

    """
    fol = None
    if ifol is None:
        intxt = input(' Input the particle number to follow \n')
        ifol = int(intxt)
    print(f"Following particle {ifol}")
    if ifol < 0: # Negative numbers are planets
        fol = ds.where(np.invert(np.isnan(ds['Gmass'])), drop=True)
        fol = fol.where(np.invert(np.isnan(fol['a'])), drop=True) # Remove times where this body doesn't exist (but this also gets rid of the central body)
        fol = fol.isel(id = -ifol - 2)  # Take 1 off for 0-indexed arrays in Python, and take 1 more off because the central body is gone
    elif ifol > 0: # Positive numbers are test particles
        fol = ds.where(np.isnan(ds['Gmass']), drop=True).drop_vars(['Gmass', 'radius'])
        fol = fol.where(np.invert(np.isnan(fol['a'])), drop=True)
        fol = fol.isel(id = ifol - 1)  # Take 1 off for 0-indexed arrays in Python

    if nskp is None:
        intxt = input('Input the print frequency\n')
        nskp = int(intxt)
        
    fol['obar'] = fol['capom'] + fol['omega']
    fol['obar'] = fol['obar'].fillna(0)
    fol['obar'] = wrap_angle(fol['obar'])
    fol['peri'] = fol['a'] * (1.0 - fol['e'])
    fol['apo']  = fol['a'] * (1.0 + fol['e'])

    
    tslice = slice(None, None, nskp)
    try:
        with open('follow.out', 'w') as f:
            print('# 1 2 3  4    5     6    7    8    9   10', file=f)
            print('# t,a,e,inc,capom,omega,capm,peri,apo,obar', file=f)
            for t in fol.isel(time=tslice).time:
                a = fol['a'].sel(time=t).values
                e = fol['e'].sel(time=t).values
                inc = fol['inc'].sel(time=t).values
                capom = fol['capom'].sel(time=t).values
                omega = fol['omega'].sel(time=t).values
                capm = fol['capm'].sel(time=t).values
                peri = fol['peri'].sel(time=t).values
                apo = fol['apo'].sel(time=t).values
                obar = fol['obar'].sel(time=t).values
                print(f"{t.values:15.7e} {a:22.16f} {e:22.16f} {inc:22.16f} {capom:22.16f} {omega:22.16f} {capm:22.16f} {peri:22.16f} {apo:22.16f} {obar:22.16f}", file=f)
                
    except IOError:
        print(f"Error writing to follow.out")
    
    return fol


def rotate_to_vector(ds, new_pole, skip_vars=['space','Ip']):
    """
    Rotates the coordinate system such that the z-axis is aligned with an input pole.  The new pole is defined by the input vector. 
    This will change all variables in the Dataset that have the "space" dimension, except for those passed to the skip_vars parameter.
    
    Parameters
    ----------
    ds : Xarray Dataset
        Dataset containing the vector quantity
    new_pole : (3) float array
        New pole vector
    skip_vars : list of str, optional
        List of variable names to skip. The default is ['space','Ip'].
        
    Returns
    -------
    ds : Xarray Dataset
        Dataset with the new pole vector applied to all variables with the "space" dimension
    """
    
    if 'space' not in ds.dims:
        print("No space dimension in Dataset")
        return ds
    
    # Verify that the new pole is a 3-element array
    if len(new_pole) != 3:
        print("New pole must be a 3-element array")
        return ds
  
    # Normalize the new pole vector to ensure it is a unit vector
    pole_mag = np.linalg.norm(new_pole)
    unit_pole = new_pole / pole_mag
    
    # Define the original and target vectors
    target_vector = np.array([0, 0, 1])  # Rotate so that the z-axis is aligned with the new pole
    original_vector = unit_pole.reshape(1, 3)  
    
    # Use align_vectors to get the rotation that aligns the z-axis with Mars_rot
    rotation, _ = R.align_vectors(target_vector, original_vector)

    # Define a function to apply the rotation, which will be used with apply_ufunc
    def apply_rotation(vector, rotation):
        return rotation.apply(vector)
    
    # Function to apply rotation to a DataArray
    def rotate_dataarray(da, rotation):
        return xr.apply_ufunc(
            apply_rotation,
            da,
            kwargs={'rotation': rotation},
            input_core_dims=[['space']],
            output_core_dims=[['space']],
            vectorize=True,
            dask='parallelized',
            output_dtypes=[da.dtype]
        )

    # Loop through each variable in the dataset and apply the rotation if 'space' dimension is present
    for var in ds.variables:
        if 'space' in ds[var].dims and var not in skip_vars:
            ds[var] = rotate_dataarray(ds[var], rotation)
            
    return ds
        
        