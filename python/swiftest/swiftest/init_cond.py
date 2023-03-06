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
from __future__ import annotations

import swiftest
import numpy as np
import numpy.typing as npt
from astroquery.jplhorizons import Horizons
import astropy.units as u
from astropy.coordinates import SkyCoord
import datetime
import xarray as xr
from typing import (
    Literal,
    Dict,
    List,
    Any
)
def solar_system_horizons(name: str,
                          param: Dict,
                          ephemerides_start_date: str,
                          id: int | None = None):
    """
    Initializes a Swiftest dataset containing the major planets of the Solar System at a particular data from JPL/Horizons

    Parameters
    ----------
    param : dict
        Swiftest paramuration parameters. This method uses the unit conversion factors to convert from JPL's AU-day system into the system specified in the param file
    ephemerides_start_date : string
        Date to use when obtaining the ephemerides in the format YYYY-MM-DD.

    Returns
    -------
    ds : xarray dataset
    """
    # Planet ids
    planetid = {
        'Sun': '0',
        'Mercury': '1',
        'Venus': '2',
        'Earth': '3',
        'Mars': '4',
        'Jupiter': '5',
        'Saturn': '6',
        'Uranus': '7',
        'Neptune': '8',
        'Pluto': '9'
    }
    
    if name in planetid:
        ispl = True
        id = planetid[name]
    else:
        ispl = False
        print(f"\nMassive body {name} not found or not yet supported")
        print("This will be created as a massless test particle")
        if id is None:
            print("ID value required for this input type")
            return
    
    # Planet MSun/M ratio
    MSun_over_Mpl = {
        'Sun'    : np.longdouble(1.0),
        'Mercury': np.longdouble(6023600.0),
        'Venus': np.longdouble(408523.71),
        'Earth': np.longdouble(328900.56),
        'Mars': np.longdouble(3098708.),
        'Jupiter': np.longdouble(1047.3486),
        'Saturn': np.longdouble(3497.898),
        'Uranus': np.longdouble(22902.98),
        'Neptune': np.longdouble(19412.24),
        'Pluto': np.longdouble(1.35e8)
    }
    
    # Planet radii in AU
    planetradius = {
        'Sun' : swiftest.RSun / swiftest.AU2M,
        'Mercury': np.longdouble(2439.4e3) / swiftest.AU2M,
        'Venus': np.longdouble(6051.8e3) / swiftest.AU2M,
        'Earth': np.longdouble(6371.0084e3) / swiftest.AU2M,  # Earth only for radius
        'Mars': np.longdouble(3389.50e3) / swiftest.AU2M,
        'Jupiter': np.longdouble(69911e3) / swiftest.AU2M,
        'Saturn': np.longdouble(58232.0e3) / swiftest.AU2M,
        'Uranus': np.longdouble(25362.e3) / swiftest.AU2M,
        'Neptune': np.longdouble(24622.e3) / swiftest.AU2M,
        'Pluto': np.longdouble(1188.3e3 / swiftest.AU2M)
    }
    
    planetrot = {
       'Sun' : np.longdouble(360.0 / 25.05) / swiftest.JD2S, # Approximate
       'Mercury': np.longdouble(360.0 / 58.646) / swiftest.JD2S,
       'Venus': np.longdouble(360.0 / 243.0226 ) / swiftest.JD2S,
       'Earth': np.longdouble(360.0 / 0.99726968) / swiftest.JD2S,
       'Mars': np.longdouble(360.0 / 1.025957) / swiftest.JD2S,
       'Jupiter': np.longdouble(360.0 / (9.9250 / 24.0) ) / swiftest.JD2S,
       'Saturn': np.longdouble(360.0 / (10.656 / 24.0) ) / swiftest.JD2S,
       'Uranus': np.longdouble(360.0 / 0.71833) / swiftest.JD2S,
       'Neptune': np.longdouble(360.0 / 0.6713) / swiftest.JD2S,
       'Pluto': np.longdouble(360.0 / 6.387230) / swiftest.JD2S
    }
   
    planetIpz = { # Only the polar moments of inertia are used currently. Where the quantity is unkown, we just use the value of a sphere = 0.4
        'Sun' : np.longdouble(0.070),
        'Mercury' : np.longdouble(0.346),
        'Venus': np.longdouble(0.4),
        'Earth': np.longdouble(0.3307),
        'Mars': np.longdouble(0.3644),
        'Jupiter': np.longdouble(0.2756),
        'Saturn': np.longdouble(0.22),
        'Uranus': np.longdouble(0.23),
        'Neptune': np.longdouble(0.23),
        'Pluto': np.longdouble(0.4)
        }

    # Unit conversion factors
    DCONV = swiftest.AU2M / param['DU2M']
    VCONV = (swiftest.AU2M / swiftest.JD2S) / (param['DU2M'] / param['TU2S'])
    THIRDLONG = np.longdouble(1.0) / np.longdouble(3.0)
    
    # Central body value vectors
    GMcb = swiftest.GMSun * param['TU2S'] ** 2 / param['DU2M'] ** 3
    Rcb = swiftest.RSun / param['DU2M']
    J2RP2 = swiftest.J2Sun * (swiftest.RSun / param['DU2M']) ** 2
    J4RP4 = swiftest.J4Sun * (swiftest.RSun / param['DU2M']) ** 4
    
    solarpole = SkyCoord(ra=286.13 * u.degree, dec=63.87 * u.degree)
    solarrot = planetrot['Sun'] * param['TU2S']
    rotcb = solarpole.cartesian * solarrot
    rotcb = np.array([rotcb.x.value, rotcb.y.value, rotcb.z.value])
    Ipsun = np.array([0.0, 0.0, planetIpz['Sun']])

    param_tmp = param
    param_tmp['OUT_FORM'] = 'XVEL'

    rh = np.full(3,np.nan)
    vh = np.full(3,np.nan)
    a = None
    e = None
    inc = None
    capom = None
    omega = None
    capm = None
    Ip = np.full(3,np.nan)
    rot = np.full(3,np.nan)
    rhill = None
    Gmass = None
    Rpl = None
    J2 = None
    J4 = None

    if name == "Sun" : # Create central body
        print("Creating the Sun as a central body")
        Gmass = GMcb
        Rpl = Rcb
        J2 = J2RP2
        J4 = J4RP4
        if param['ROTATION']:
            Ip = Ipsun
            rot = rotcb
    else: # Fetch solar system ephemerides from Horizons
        print(f"Fetching ephemerides data for {name} from JPL/Horizons")

        # Horizons date time internal variables
        tstart = datetime.date.fromisoformat(ephemerides_start_date)
        tstep = datetime.timedelta(days=1)
        tend = tstart + tstep
        ephemerides_end_date = tend.isoformat()
        ephemerides_step = '1d'

        pldata = {}
        pldata[name] = Horizons(id=id, location='@sun',
                               epochs={'start': ephemerides_start_date, 'stop': ephemerides_end_date,
                                       'step': ephemerides_step})
        
        if param['IN_FORM'] == 'XV':
            rx = pldata[name].vectors()['x'][0] * DCONV
            ry = pldata[name].vectors()['y'][0] * DCONV
            rz = pldata[name].vectors()['z'][0] * DCONV
            vx = pldata[name].vectors()['vx'][0] * VCONV
            vy = pldata[name].vectors()['vy'][0] * VCONV
            vz = pldata[name].vectors()['vz'][0] * VCONV

            rh = np.array([rx,ry,rz])
            vh = np.array([vx,vy,vz])
        elif param['IN_FORM'] == 'EL':
            a = pldata[name].elements()['a'][0] * DCONV
            e = pldata[name].elements()['e'][0]
            inc = pldata[name].elements()['incl'][0]
            capom = pldata[name].elements()['Omega'][0]
            omega = pldata[name].elements()['w'][0]
            capm = pldata[name].elements()['M'][0]

        if ispl:
            Gmass = GMcb / MSun_over_Mpl[name]
            if param['CHK_CLOSE']:
                Rpl = planetradius[name] * DCONV

            # Generate planet value vectors
            if (param['RHILL_PRESENT']):
                rhill = pldata[name].elements()['a'][0] * DCONV * (3 * MSun_over_Mpl[name]) ** (-THIRDLONG)

            if (param['ROTATION']):
                RA = pldata[name].ephemerides()['NPole_RA'][0]
                DEC = pldata[name].ephemerides()['NPole_DEC'][0]

                rotpole = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree)
                rotrate = planetrot[name] * param['TU2S']
                rot = rotpole.cartesian * rotrate
                rot = np.array([rot.x.value, rot.y.value, rot.z.value])
                Ip = np.array([0.0, 0.0, planetIpz[name]])

        else:
            Gmass = None

    if id is None:
        id = planetid[name]

    return id,name,a,e,inc,capom,omega,capm,rh,vh,Gmass,Rpl,rhill,Ip,rot,J2,J4

def vec2xr(param: Dict, **kwargs: Any):
    """
    Converts and stores the variables of all bodies in an xarray dataset.

    Parameters
    ----------
    param : dict
        Swiftest paramuration parameters.
    name : str or array-like of str, optional
        Name or names of Bodies. If none passed, name will be "Body<id>"
    id : int or array-like of int, optional
        Unique id values. If not passed, an id will be assigned in ascending order starting from the pre-existing
        Dataset ids.
    a : float or array-like of float, optional
        semimajor axis for param['IN_FORM'] == "EL"
    e : float or array-like of float, optional
        eccentricity  for param['IN_FORM'] == "EL"
    inc : float or array-like of float, optional
        inclination for param['IN_FORM'] == "EL"
    capom : float or array-like of float, optional
        longitude of periapsis for param['IN_FORM'] == "EL"
    omega : float or array-like of float, optional
        argument of periapsis for param['IN_FORM'] == "EL"
    capm : float or array-like of float, optional
        mean anomaly for param['IN_FORM'] == "EL"
    rh : (n,3) array-like of float, optional
        Position vector array. This can be used instead of passing v1, v2, and v3 sepearately for "XV" input format
    vh : (n,3) array-like of float, optional
        Velocity vector array. This can be used instead of passing v4, v5, and v6 sepearately for "XV" input format
    Gmass : float or array-like of float, optional
        G*mass values if these are massive bodies (only one of mass or Gmass can be passed)
    radius : float or array-like of float, optional
        Radius values if these are massive bodies
    rhill : float or array-like of float, optional
        Hill's radius values if these are massive bodies
    rot:  (n,3) array-like of float, optional
        Rotation rate vectors if these are massive bodies with rotation enabled in deg/TU
    Ip: (n,3) array-like of flaot, optional
        Principal axes moments of inertia vectors if these are massive bodies with rotation enabled. This can be used
        instead of passing Ip1, Ip2, and Ip3 separately
    time : array of floats
        Time at start of simulation
    Returns
    -------
    ds : xarray dataset
    """
    scalar_dims = ['id']
    vector_dims = ['id','space']
    space_coords = np.array(["x","y","z"])

    vector_vars = ["rh","vh","Ip","rot"]
    scalar_vars = ["name","a","e","inc","capom","omega","capm","Gmass","radius","rhill","J2","J4"]
    time_vars =  ["rh","vh","Ip","rot","a","e","inc","capom","omega","capm","Gmass","radius","rhill","J2","J4"]

    # Check for valid keyword arguments
    kwargs = {k:kwargs[k] for k,v in kwargs.items() if v is not None}
    if "rot" not in kwargs and "Gmass" in kwargs:
        kwargs['rot'] = np.zeros((len(kwargs['Gmass']),3))
    if "Ip" not in kwargs and "Gmass" in kwargs:
        kwargs['Ip'] = np.full((len(kwargs['Gmass']),3), 0.4)

    if "time" not in kwargs:
        kwargs["time"] = np.array([0.0])

    valid_arguments = vector_vars + scalar_vars + ['time','id']

    kwargs = {k:v for k,v in kwargs.items() if k in valid_arguments}

    data_vars = {k:(scalar_dims,v) for k,v in kwargs.items() if k in scalar_vars}
    data_vars.update({k:(vector_dims,v) for k,v in kwargs.items() if k in vector_vars})
    ds = xr.Dataset(data_vars=data_vars,
                    coords={
                        "id":(["id"],kwargs['id']),
                        "space":(["space"],space_coords),
                    }
                    )
    time_vars = [v for v in time_vars if v in ds]
    for v in time_vars:
        ds[v] = ds[v].expand_dims({"time":1}).assign_coords({"time": kwargs['time']})

    return ds