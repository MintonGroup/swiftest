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

import astropy.constants as const
import datetime

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
CB_TYPE_NAME = "Central Body"
PL_TYPE_NAME = "Massive Body"
TP_TYPE_NAME = "Test Particle"
PL_TINY_TYPE_NAME = "Semi-Interacting Massive Body"

# The default value is Prof. Minton's Brimley/Cocoon line crossing date (aka MBCL)
_mbday = datetime.date.fromisoformat('1976-08-05')
_bcl = datetime.timedelta(days=18530)
MINTON_BCL = (_mbday + _bcl).isoformat()
