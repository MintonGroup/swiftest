"""
 Copyright 2023 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, Dana Singh, & Kaustub Anand
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest. 
 If not, see: https://www.gnu.org/licenses. 
"""

# python functions to read in and set up spherical harmonics coefficients for non-standard gravity terms
# using pySHTOOLS (in python) by Mark W.
# 

import swiftest
import numpy as np
from astroquery.jplhorizons import Horizons
import astropy.units as u
from astropy.coordinates import SkyCoord
import datetime
import xarray as xr
import pyshtools as pysh
from typing import (
    Literal,
    Dict,
    List,
    Any
)

def clm_from_ellipsoid(mass, density, a, b = None, c = None, lmax = 6):
    """
    Creates and returns the gravity coefficients for an ellipsoid with principal axes a, b, c upto a certain maximum degree lmax. 
    Uses pyshtools. No units necessary for a, b, & c. However, they need to be in the same units (DU).

    Parameters
    ----------
    mass : float
        mass of the central body
    density : float
        density of the central body
    a : float 
        length of the pricipal axis aligned with the x axis
    b : float, optional, default = a
        length of the pricipal axis aligned with the y axis
    c : float, optional, default = b
        length of the pricipal axis aligned with the z axis
    lmax : int, optional, default = 6
        The maximum spherical harmonic degree resolvable by the grid.

    Returns
    -------
    clm : ndarry, shape (2, lmax+1, lmax+1)
        numpy ndarray of the gravitational potential spherical harmonic coefficients. 
        This is a three-dimensional array formatted as coeffs[i, degree, order], 
        where i=0 corresponds to positive orders and i=1 to negative orders.

    """
    G = swiftest.constants.GC
    Gmass = G * mass # SHTOOLS uses an SI G value, and divides it before using the mass
            # FIND A BETTER WAY TO CHANGE UNITS

    # cap lmax to ensure fast performance without giving up accuracy
    
    lmax_limit = 6 # lmax_limit = 6 derived from Jean's Law by taking the characteristic wavelength as the radius of the CB
    if(lmax > lmax_limit):                              # FIND A BETTER WAY to judge this cut off point, i.e., relative change between coefficients
        lmax = lmax_limit
        print(f'Setting maximum spherical harmonic degree to {lmax_limit}')

    # create shape grid and convert to Spherical Harmonics (.expand())
    shape_SH = pysh.SHGrid.from_ellipsoid(lmax = lmax, a = a, b = b, c = c).expand()

    # get coefficients
    clm_class = pysh.SHGravcoeffs.from_shape(shape_SH, rho = density, gm = Gmass) # 4pi normalization
    clm = clm_class.to_array()

    return clm

def clm_from_relief(lmax, mass, density, grid):
    """
    Creates and returns the gravity coefficients for a body with a given DH grid upto a certain maximum degree lmax. 
    Uses pyshtools.

    Parameters
    ----------
    lmax : int
        The maximum spherical harmonic degree resolvable by the grid.
    mass : float
        mass of the central body
    density : float
        density of the central body
    grid : array, shape []
        DH grid of the surface relief of the body

    Returns
    -------
    clm : ndarry, shape (2, lmax+1, lmax+1)
        numpy ndarray of the gravitational potential spherical harmonic coefficients. 
        This is a three-dimensional array formatted as coeffs[i, degree, order], 
        where i=0 corresponds to positive orders and i=1 to negative orders.

    """

    G = swiftest.constants.GC
    Gmass = G * mass # SHTOOLS uses an SI G value, and divides it before using the mass
            # FIND A BETTER WAY TO CHANGE UNITS

    # cap lmax to 20 to ensure fast performance
    lmax_limit = 20
    if(lmax > lmax_limit): # FIND A BETTER WAY to judge this cut off point, i.e., relative change between coefficients
        lmax = lmax_limit
        print(f'Setting maximum spherical harmonic degree to {lmax_limit}')

    # convert to spherical harmonics
    shape_SH = pysh.SHGrid.from_array(grid).expand()
        # shape_SH = SHExpandDH(grid)

    # get coefficients
    clm_class = pysh.SHGravcoeffs.from_shape(shape_SH, rho = density, gm = Gmass) # 4pi normalization
    clm = clm_class.to_array()

    return clm
