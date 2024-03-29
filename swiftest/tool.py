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
        fol = ds.where(np.isnan(ds['Gmass']), drop=True).drop_vars(['Gmass', 'radius'],errors="ignore")
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
