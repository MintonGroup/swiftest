"""
Copyright 2025 - David Minton.

This file is part of Swiftest.
Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Swiftest.
If not, see: https://www.gnu.org/licenses.
"""

import datetime

import astropy.constants as const
import numpy as np

# Constants in SI units
GC = np.longdouble(const.G.value[()])
AU2M = np.longdouble(const.au.value)
GMSun = np.longdouble(const.GM_sun.value)
MSun = np.longdouble(const.M_sun.value)
RSun = np.longdouble(const.R_sun.value)
MEarth = np.longdouble(const.M_earth.value)
REarth = np.longdouble(const.R_earth.value)
GMEarth = np.longdouble(const.GM_earth.value)
JD2S = 86400
YR2S = np.longdouble(365.25 * JD2S)
einsteinC = np.longdouble(const.c)
CB_TYPE_NAME = "Central Body"
PL_TYPE_NAME = "Massive Body"
TP_TYPE_NAME = "Test Particle"
PL_TINY_TYPE_NAME = "Semi-Interacting Massive Body"

# The default value is Prof. Minton's Brimley/Cocoon line crossing date (aka MBCL)
_mbday = datetime.date.fromisoformat("1976-08-05")
_bcl = datetime.timedelta(days=18530)
MINTON_BCL = (_mbday + _bcl).isoformat()
