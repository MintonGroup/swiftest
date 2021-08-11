import swiftest
import numpy as np
import os
import glob
import xarray as xr
"""
Functions that recreate the Swift/Swifter tool programs
"""

def wrap_angle(angle):
    while np.any(angle >= 2 * np.pi):
        angle[angle >= 2 * np.pi] -= 2 * np.pi
    while np.any(angle < 0.0):
        angle[angle < 0.0] += 2 * np.pi
    return angle

def follow_swift(ds, ifol=None, nskp=None):
    """
    Emulates the Swift version of tool_follow.f


    Parameters
    ----------
    ds : Dataset containing orbital elements

    Returns
    -------
    fol : Dataset containing only the followed body with angles converted to degrees
    
    Generates a file called follow.out containing 10 columns (angles are all in degrees):
        1 2 3  4    5     6    7    8    9   10
        t,a,e,inc,capom,omega,capm,peri,apo,obar

    """
    fol = None
    if ifol is None:
        intxt = input(' Input the particle number to follow \n')
        ifol = int(intxt)
    print(f"Following particle {ifol}")
    if ifol < 0: # Negative numbers are planets
        fol = ds.where(np.invert(np.isnan(ds['GMass'])), drop=True)
        fol = fol.where(np.invert(np.isnan(fol['a'])), drop=True) # Remove times where this body doesn't exist (but this also gets rid of the central body)
        fol = fol.isel(id = -ifol - 2)  # Take 1 off for 0-indexed arrays in Python, and take 1 more off because the central body is gone
    elif ifol > 0: # Positive numbers are test particles
        fol = ds.where(np.isnan(ds['GMass']), drop=True).drop_vars(['GMass', 'Radius'])
        fol = fol.where(np.invert(np.isnan(fol['a'])), drop=True)
        fol = fol.isel(id = ifol - 1)  # Take 1 off for 0-indexed arrays in Python

    if nskp is None:
        intxt = input('Input the print frequency\n')
        nskp = int(intxt)
        
    dr = 180.0 / np.pi
    fol['obar'] = fol['capom'] + fol['omega']
    fol['obar'] = fol['obar'].fillna(0)
    fol['obar'] = wrap_angle(fol['obar'])
    fol['obar'] = fol['obar'] * dr
    fol['inc'] = fol['inc'] * dr
    fol['capom'] = fol['capom'] * dr
    fol['omega'] = fol['omega'] * dr
    fol['capm'] = fol['capm'] * dr
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