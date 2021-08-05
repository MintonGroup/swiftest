import numpy as np
import astropy.constants as const

# Constants in SI units
GC = np.longdouble(const.G.value)
AU2M = np.longdouble(const.au.value)
GMSunSI = np.longdouble(const.GM_sun.value)
MSun = np.longdouble(const.M_sun.value)
RSun = np.longdouble(const.R_sun.value)
JD2S = 86400
YR2S = np.longdouble(365.25 * JD2S)
einsteinC = np.longdouble(299792458.0)
# Solar oblatenes values: From Mecheri et al. (2004), using Corbard (b) 2002 values (Table II)
J2Sun = np.longdouble(2.198e-7)
J4Sun = np.longdouble(-4.805e-9)
