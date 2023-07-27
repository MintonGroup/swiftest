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
import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord

# Constants in SI units
GC = const.G.value[()]
AU2M = const.au.value
GMSun = const.GM_sun.value
MSun = const.M_sun.value
RSun = const.R_sun.value
MEarth = const.M_earth.value
REarth = const.R_earth.value
GMEarth = const.GM_earth.value
JD2S = 86400
YR2S = 365.25 * JD2S
einsteinC = 299792458.0
# Solar oblatenes values: From Mecheri et al. (2004), using Corbard (b) 2002 values (Table II)
J2Sun = 2.198e-7
J4Sun = -4.805e-9
rotpoleSun = SkyCoord(ra=286.13 * u.degree, dec=63.87 * u.degree).cartesian
rotSun = (360.0 / 25.05) / JD2S  * rotpoleSun 

