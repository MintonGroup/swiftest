"""
Copyright 2024 - The Minton Group at Purdue University
This file is part of Swiftest.
Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Swiftest. 
If not, see: https://www.gnu.org/licenses. 
"""

from .constants import *
from .simulation import Simulation
from .shgrav import clm_from_ellipsoid, clm_from_relief
from .data import SwiftestDataArray, SwiftestDataset
from .init_cond import get_solar_system_body, get_solar_system_body_mass_rotation
from . import core
