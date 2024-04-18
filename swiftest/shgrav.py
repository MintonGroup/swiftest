"""
 Copyright 2024 - Minton Group at Purdue University
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest. 
 If not, see: https://www.gnu.org/licenses. 
"""

# python functions to read in and set up spherical harmonics coefficients for non-standard gravity terms
# using pySHTOOLS see: https://shtools.github.io/SHTOOLS/
# 

from .constants import GC
import os

if 'READTHEDOCS' in os.environ:
    # Assume pyshtools is available when building on ReadTheDocs
    PYSHTOOLS_AVAILABLE = True
else:
    try:
        import pyshtools as pysh
        PYSHTOOLS_AVAILABLE = True
    except ModuleNotFoundError:
        PYSHTOOLS_AVAILABLE = False


if PYSHTOOLS_AVAILABLE:

    def clm_from_ellipsoid(mass, density, a, b = None, c = None, lmax = 6, lref_radius = False, ref_radius = None):
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
        lref_radius : boolean, optional, default = False
            Boolean value to return the reference radius calculated by SHTOOLS
        ref_radius : float, optional, default = None
            Reference radius to scale the gravitational coefficients to

        Returns
        -------
        clm : ndarry, shape (2, lmax+1, lmax+1)
            numpy ndarray of the gravitational potential spherical harmonic coefficients. 
            This is a three-dimensional array formatted as coeffs[i, degree, order], 
            where i=0 corresponds to positive orders and i=1 to negative orders.

        """
        Gmass = GC * mass # SHTOOLS uses an SI G value, and divides it before using the mass; NO NEED TO CHANGE UNITS

        # cap lmax to ensure fast performance without giving up accuracy
        lmax_limit = 6              # lmax_limit = 6 derived from Jean's Law; characteristic wavelength = the radius of the CB
        if(lmax > lmax_limit):                           
            lmax = lmax_limit
            print(f'Setting maximum spherical harmonic degree to {lmax_limit}')

        # create shape grid 
        shape_SH = pysh.SHGrid.from_ellipsoid(lmax = lmax, a = a, b = b, c = c)

        # get gravity coefficients
        clm_class = pysh.SHGravCoeffs.from_shape(shape_SH, rho = density, gm = Gmass) # 4pi normalization
        clm = clm_class.to_array(normalization = '4pi') # export as array with 4pi normalization and not scaling by 4*pi to match normalisation

        # Return reference radius EQUALS the radius of the Central Body
        print(f'Ensure that the Central Body radius equals the reference radius.')

        if(lref_radius == True and ref_radius is None):
            ref_radius = shape_SH.expand(normalization = '4pi').coeffs[0, 0, 0]
            return clm, ref_radius
        elif(lref_radius == True and ref_radius is not None):
            clm_class = clm_class.change_ref(r0 = ref_radius)
            clm = clm_class.to_array(normalization = '4pi')
            return clm, ref_radius
        else:
            return clm

    def clm_from_relief(mass, density, grid, lmax = 6, lref_radius = False, ref_radius = None):
        """
        Creates and returns the gravity coefficients for a body with a given DH grid upto a certain maximum degree lmax. 
        Uses pyshtools.

        Parameters
        ----------
        mass : float
            mass of the central body
        density : float
            density of the central body
        grid : array, shape []
            DH grid of the surface relief of the body
        lmax : int, optional, default = 6
            The maximum spherical harmonic degree resolvable by the grid.
        lref_radius : boolean, optional, default = False
            Boolean value to return the reference radius calculated by SHTOOLS
        ref_radius : float, optional, default = None
            Reference radius to scale the gravitational coefficients to

        Returns
        -------
        clm : ndarry, shape (2, lmax+1, lmax+1)
            numpy ndarray of the gravitational potential spherical harmonic coefficients. 
            This is a three-dimensional array formatted as coeffs[i, degree, order], 
            where i=0 corresponds to positive orders and i=1 to negative orders.

        """

        Gmass = GC * mass # SHTOOLS uses an SI G value, and divides it before using the mass; NO NEED TO CHANGE UNITS

        # cap lmax to 20 to ensure fast performance
        lmax_limit = 6
        if(lmax > lmax_limit): # FIND A BETTER WAY to judge this cut off point, i.e., relative change between coefficients
            lmax = lmax_limit
            print(f'Setting maximum spherical harmonic degree to {lmax_limit}')

        # convert to spherical harmonics
        shape_SH = pysh.SHGrid.from_array(grid)

        # get coefficients
        clm_class = pysh.SHGravcoeffs.from_shape(shape_SH, rho = density, gm = Gmass) # 4pi normalization
        clm = clm_class.to_array(normalization = '4pi') # export as array with 4pi normalization

        # Return reference radius EQUALS the radius of the Central Body

        print(f'Ensure that the Central Body radius equals the reference radius.')

        if(lref_radius == True and ref_radius is None):
            ref_radius = shape_SH.expand(normalization = '4pi').coeffs[0, 0, 0]
            return clm, ref_radius
        elif(lref_radius == True and ref_radius is not None):
            clm_class = clm_class.change_ref(r0 = ref_radius)
            clm = clm_class.to_array(normalization = '4pi')
            return clm, ref_radius
        else:
            return clm

else:
    def clm_from_ellipsoid(*args, **kwargs):
        raise NotImplementedError("Sph_Harmonics is not available because pyshtools is not installed.")
    def clm_from_relief(*args, **kwargs):
        raise NotImplementedError("Sph_Harmonics is not available because pyshtools is not installed.")
    