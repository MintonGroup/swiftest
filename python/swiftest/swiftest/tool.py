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
"""
Functions that recreate the Swift/Swifter tool programs
"""
def magnitude(ds,x):
    dim = "space"
    ord = None
    return xr.apply_ufunc(
        np.linalg.norm, ds[x].where(~np.isnan(ds[x])), input_core_dims=[[dim]], kwargs={"ord": ord, "axis": -1}, dask="allowed"
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

def danby(M, ecc, accuracy=1e-14):
    """
    Danby's method to solve Kepler's equation.

    Parameters
    ----------
    M : float
        the mean anomaly in radians
    ecc : float
        the eccentricity
    accuracy : float
        the relative accuracy to obtain a solution. Default is 1e-12.

    Returns
    ----------
    E : float
        the eccentric anomaly in radians

    References
    __________
    Danby, J.M.A. 1988. Fundamentals of celestial mechanics. Richmond, Va., U.S.A., Willmann-Bell, 1988. 2nd ed.
    Murray, C.D., Dermott, S.F., 1999. Solar system dynamics, New York, New York. ed, Cambridge University Press.

    """

    def kepler_root(E, ecc, M, deriv):
        """
        The Kepler equation root function.

        Parameters
        ----------
        E : float
            the eccentric anomaly in radians
        ecc : float
            the eccentricity
        M : float
            the mean anomaly in radians
        deriv : int between 0 and 3
            The order of the derivative to compute

        Returns
        ----------
        deriv = 0: E - e * np.sin(E) - M
        deriv = 1: 1 - e * np.cos(E)
        deriv = 2: e * np.sin(E)
        deriv = 3: e * np.cos(E)

        Note: The function will return 0 when E is correct for a given e and M
        """

        if deriv == 0:
            return E - ecc * np.sin(E) - M
        elif deriv == 1:
            return 1.0 - ecc * np.cos(E)
        elif deriv == 2:
            return ecc * np.sin(E)
        elif deriv == 3:
            return ecc * np.cos(E)
        else:
            print(f"deriv={deriv} is undefined")
            return None

    def delta_ij(E, ecc, M, j):
        """
        Danby's intermediate delta functions. This function is recursive for j>1.

        Parameters
        ----------
        E : float
            the eccentric anomaly in radians
        ecc : float
            the eccentricity
        M : float
            the mean anomaly in radians
        j : int between 1 and 3
            The order of the delta function to compute

        Returns
        ----------
        delta_ij value used in Danby's iterative Kepler equation solver

        """
        if j == 1:
            return -kepler_root(E, ecc, M, 0) / kepler_root(E, ecc, M, 1)
        elif j == 2:
            return -kepler_root(E, ecc, M, 0) / (kepler_root(E, ecc, M, 1)
                                                 - delta_ij(E, ecc, M, 1) * kepler_root(E, ecc, M, 2) / 2.0)
        elif j == 3:
            return -kepler_root(E, ecc, M, 0) / (kepler_root(E, ecc, M, 1)
                                                 + delta_ij(E, ecc, M, 1) * kepler_root(E, ecc, M, 2) / 2.0
                                                 + delta_ij(E, ecc, M, 2) ** 2 * kepler_root(E, ecc, M,
                                                                                             3) / 6.0)
        else:
            print(f"j = {j} is not a valid input to the delta_ij function")

    # If this is a circular orbit, just return the mean anomaly as the eccentric anomaly
    if ecc < np.finfo(np.float64).tiny:  # Prevent floating point exception for 0 eccentricity orbits
        return M

    # Initial guess
    k = 0.85
    E = M + np.sign(np.sin(M)) * k * ecc
    MAXLOOPS = 50  # Maximum number of times to iterate before we give up
    for i in range(MAXLOOPS):
        Enew = E + delta_ij(E, ecc, M, 3)
        if np.abs((Enew - E) / E) < accuracy:
            return Enew
        E = Enew

    raise RuntimeError("The danby function did not converge on a solution.")



def el2xv_one(mu, a, ecc, inc, Omega, omega, M):
    """
    Compute osculating orbital elements from relative Cartesian position and velocity
    All angular measures are returned in radians
        If inclination < TINY, longitude of the ascending node is arbitrarily set to 0

        If eccentricity < sqrt(TINY), argument of pericenter is arbitrarily set to 0

          ALGORITHM:  See Fitzpatrick "Principles of Cel. Mech."

    Adapted from Martin Duncan's el2xv.f

    Parameters
    ----------
    mu : float
        Central body gravitational constant
    a : float
        semimajor axis
    ecc : float
        eccentricity
    inc : float
        inclination (degrees)
    Omega : float
        longitude of ascending node (degrees)
    omega : float
        argument of periapsis (degrees)
    M : flaot
        Mean anomaly (degrees)

    Returns
    ----------
    rvec, vvec  : tuple of float vectors

    rvec : (3) float vector
        Cartesian position vector
    vvec : (3) float vector
        Cartesian velocity vector

    """

    if ecc < 0.0:
        print("Error in el2xv! Eccentricity cannot be negative. Setting it to 0")
        e = 0.0
        iorbit_type = "ellipse"
    else:
        e = ecc
        em1 = e - 1.
        if np.abs(em1) < 2 * np.finfo(float).eps:
            iorbit_type = "parabola"
        elif e > 1.0:
            iorbit_type = "hyperbola"
        else:
            iorbit_type = "ellipse"

    def scget(angle):
        """
        Efficiently compute the sine and cosine of an input angle
        Input angle must be in radians

        Adapted from David E. Kaufmann's Swifter routine: orbel_scget.f90
        Adapted from Hal Levison's Swift routine orbel_scget.f

        Parameters
        ----------
        angle : input angle

        Returns
        -------
        sx : sin of angle
        cx : cos of angle

        """
        TWOPI = 2 * np.pi
        nper = int(angle / TWOPI)
        x = angle - nper * TWOPI
        if x < 0.0:
            x += TWOPI

        sx = np.sin(x)
        cx = np.sqrt(1.0 - sx ** 2)
        if x > np.pi / 2.0 and x < 3 * np.pi / 2.0:
            cx = -cx
        return sx, cx

    sip, cip = scget(np.deg2rad(omega))
    so, co = scget(np.deg2rad(Omega))
    si, ci = scget(np.deg2rad(inc))

    d11 = cip * co - sip * so * ci
    d12 = cip * so + sip * co * ci
    d13 = sip * si
    d21 = -sip * co - cip * so * ci
    d22 = -sip * so + cip * co * ci
    d23 = cip * si

    # Get the other quantities depending on orbit type
    if iorbit_type == "ellipse":
        E = danby(np.deg2rad(M),e)
        scap, ccap = scget(E)
        sqe = np.sqrt(1. - e**2)
        sqgma = np.sqrt(mu * a)
        xfac1 = a * (ccap - e)
        xfac2 = a * sqe * scap
        ri = 1. / (a * (1. - e * ccap))
        vfac1 = -ri * sqgma * scap
        vfac2 = ri * sqgma * sqe * ccap
    else:
        print(f"Orbit type: {iorbit_type} not yet implemented")
        xfac1 = 0.0
        xfac2 = 0.0
        vfac1 = 0.0
        vfac2 = 0.0

    rvec = np.array([d11 * xfac1 + d21 * xfac2,
                     d12 * xfac1 + d22 * xfac2,
                     d13 * xfac1 + d23 * xfac2])
    vvec = np.array([d11 * vfac1 + d21 * vfac2,
                     d12 * vfac1 + d22 * vfac2,
                     d13 * vfac1 + d23 * vfac2])

    return rvec, vvec


def el2xv_vec(mu, a, ecc, inc, Omega, omega, M):
    """

    Vectorized call to el2xv_one
    Parameters
    ----------
    mu : float
        Central body gravitational constant
    a : (n) float array
        semimajor axis
    ecc : (n) float array
        eccentricity
    inc : (n) float array
        inclination (degrees)
    Omega : (n) float array
        longitude of ascending node (degrees)
    omega : (n) float array
        argument of periapsis (degrees)
    M : (n) float array
        Mean anomaly (degrees)

    Returns
    ----------
    rvec, vvec  : tuple of float vectors

    rvec : (n,3) float rray
        Cartesian position vector
    vvec : (n,3) float array
        Cartesian velocity vector
    """
    vecfunc = np.vectorize(el2xv_one, signature='(),(),(),(),(),(),()->(3),(3)')
    return vecfunc(mu, a, ecc, inc, Omega, omega, M)

def xv2el_one(mu,rvec,vvec):
    """
    Converts from cartesian position and velocity values to orbital elements

    Parameters
    ----------
    mu : float
        Central body gravitational constant
    rvec : (3) float array
        Cartesian position vector
    vvec : (3) float array
        Cartesian velocity vector

    Returns
    ----------
    a, ecc, inc, Omega, omega, M, varpi, f, lam : tuple of floats

    a : float
        semimajor axis
    ecc : float
        eccentricity
    inc : float
        inclination (degrees)
    Omega : float
        longitude of ascending node (degrees)
    omega : float
        argument of periapsis (degrees)
    M : float
        mean anomaly (degrees)
    varpi : flaot
        longitude of periapsis (degrees)
    f : float
        true anomaly (degrees)
    lam : float
        mean longitude (degrees)

    """

    rmag = np.linalg.norm(rvec)
    vmag2 = np.vdot(vvec,vvec)
    h = np.cross(rvec,vvec)
    hmag = np.linalg.norm(h)

    rdot = np.sign(np.vdot(rvec,vvec)) * np.sqrt(vmag2 - (hmag / rmag)**2)

    a = 1.0/(2.0 / rmag - vmag2/mu)
    ecc = np.sqrt(1 - hmag**2 / (mu * a))
    inc = np.arccos(h[2]/hmag)

    goodinc = np.abs(inc) > np.finfo(np.float64).tiny
    sO = np.where(goodinc,  np.sign(h[2]) * h[0] / (hmag * np.sin(inc)),0.0)
    cO = np.where(goodinc, -np.sign(h[2]) * h[1] / (hmag * np.sin(inc)),1.0)

    Omega = np.arctan2(sO, cO)

    sof = np.where(goodinc,rvec[2] / (rmag * np.sin(inc)), rvec[1]/rmag)
    cof = np.where(goodinc,(rvec[0] / rmag + np.sin(Omega) * sof * np.cos(inc)) / np.cos(Omega), rvec[0]/rmag)

    of = np.arctan2(sof,cof)

    goodecc = ecc > np.finfo(np.float64).tiny
    sf = np.where(goodecc,a * (1.0 - ecc**2)  * rdot / (hmag * ecc), 0.0)
    cf = np.where(goodecc,(1.0 / ecc) * (a * (1.0 - ecc**2) / rmag - 1.0), 1.0 )

    f = np.arctan2(sf, cf)

    omega = of - f

    varpi = Omega + omega

    # Compute eccentric anomaly & mean anomaly in order to get mean longitude
    E = np.where(ecc > np.finfo(np.float64).tiny, np.arccos(-(rmag - a) / (a * ecc)), 0)
    if np.sign(np.vdot(rvec, vvec)) < 0.0:
        E = 2 * np.pi - E

    M = E - ecc * np.sin(E)
    lam = M + varpi

    return a, ecc, np.rad2deg(inc), np.rad2deg(Omega), np.rad2deg(omega), np.rad2deg(M), np.rad2deg(varpi), np.rad2deg(f), np.rad2deg(lam)

def xv2el_vec(mu, rvec, vvec):
    """
    Vectorized call to xv2el_one.

    Parameters
    ----------
    mu : float
        Central body gravitational constant
    rvec : (n,3) float array
        Cartesian position vector
    vvec : (n,3) float array
        Cartesian velocity vector

    Returns
    ----------
    a, ecc, inc, Omega, omega, M, varpi, f, lam : tuple of float arrays

    a : (n) float array
        semimajor axis
    ecc : (n) float array
        eccentricity
    inc : (n) float array
        inclination (degrees)
    Omega : (n) float array
        longitude of ascending node (degrees)
    omega : (n) float array
        argument of periapsis (degrees)
    M : (n) float array
        mean anomaly (degrees)
    varpi : (n) flaot array
        longitude of periapsis (degrees)
    f : (n) float array
        true anomaly (degrees)
    lam : (n) float array
        mean longitude (degrees)

    """

    vecfunc = np.vectorize(xv2el_one, signature='(),(3),(3)->(),(),(),(),(),(),(),(),()')
    return vecfunc(mu, rvec, vvec)