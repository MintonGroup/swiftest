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
def solar_system_horizons(plname: str,
                          param: Dict,
                          ephemerides_start_date: str,
                          idval: int | None = None):
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
    
    if plname in planetid:
        ispl = True
        idval = planetid[plname]
    else:
        ispl = False
        print(f"\nMassive body {plname} not found or not yet supported")
        print("This will be created as a massless test particle")
        if idval is None:
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
       'Sun' : np.longdouble(2*np.pi / 25.05) / swiftest.JD2S, # Approximate
       'Mercury': np.longdouble(2*np.pi / 58.646) / swiftest.JD2S,
       'Venus': np.longdouble(2*np.pi / 243.0226 ) / swiftest.JD2S,
       'Earth': np.longdouble(2*np.pi / 0.99726968) / swiftest.JD2S,
       'Mars': np.longdouble(2*np.pi / 1.025957) / swiftest.JD2S,
       'Jupiter': np.longdouble(2*np.pi / (9.9250 / 24.0) ) / swiftest.JD2S,
       'Saturn': np.longdouble(2*np.pi / (10.656 / 24.0) ) / swiftest.JD2S,
       'Uranus': np.longdouble(2*np.pi / 0.71833) / swiftest.JD2S,
       'Neptune': np.longdouble(2*np.pi / 0.6713) / swiftest.JD2S,
       'Pluto': np.longdouble(2*np.pi / 6.387230) / swiftest.JD2S
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
    Ipsun = np.array([0.0, 0.0, planetIpz['Sun']])

    param_tmp = param
    param_tmp['OUT_FORM'] = 'XVEL'

    if plname == "Sun" : # Create central body
        print("Creating the Sun as a central body")
        v1 = None
        v2 = None
        v3 = None
        v4 = None
        v5 = None
        v6 = None
        rhill = None
        GMpl = GMcb
        Rpl = Rcb
        J2 = J2RP2
        J4 = J4RP4
        if param['ROTATION']:
            Ip1 = Ipsun[0]
            Ip2 = Ipsun[1]
            Ip3 = Ipsun[2]
            rotx = rotcb.x.value
            roty = rotcb.y.value
            rotz = rotcb.z.value
        else:
            Ip1 = None
            Ip2 = None
            Ip3 = None
            rotx = None
            roty = None
            rotz = None
    else: # Fetch solar system ephemerides from Horizons
        print(f"Fetching ephemerides data for {plname} from JPL/Horizons")

        # Horizons date time internal variables
        tstart = datetime.date.fromisoformat(ephemerides_start_date)
        tstep = datetime.timedelta(days=1)
        tend = tstart + tstep
        ephemerides_end_date = tend.isoformat()
        ephemerides_step = '1d'

        J2 = None
        J4 = None

        pldata = {}
        pldata[plname] = Horizons(id=idval, location='@sun',
                               epochs={'start': ephemerides_start_date, 'stop': ephemerides_end_date,
                                       'step': ephemerides_step})
        
        if param['IN_FORM'] == 'XV':
            v1 = pldata[plname].vectors()['x'][0] * DCONV
            v2 = pldata[plname].vectors()['y'][0] * DCONV
            v3 = pldata[plname].vectors()['z'][0] * DCONV
            v4 = pldata[plname].vectors()['vx'][0] * VCONV
            v5 = pldata[plname].vectors()['vy'][0] * VCONV
            v6 = pldata[plname].vectors()['vz'][0] * VCONV
        elif param['IN_FORM'] == 'EL':
            v1 = pldata[plname].elements()['a'][0] * DCONV
            v2 = pldata[plname].elements()['e'][0]
            v3 = pldata[plname].elements()['incl'][0]
            v4 = pldata[plname].elements()['Omega'][0]
            v5 = pldata[plname].elements()['w'][0]
            v6 = pldata[plname].elements()['M'][0]

        if ispl:
            GMpl = GMcb / MSun_over_Mpl[plname]
            if param['CHK_CLOSE']:
                Rpl = planetradius[plname] * DCONV
            else:
                Rpl = None

            # Generate planet value vectors
            if (param['RHILL_PRESENT']):
                rhill = pldata[plname].elements()['a'][0] * DCONV * (3 * MSun_over_Mpl[plname]) ** (-THIRDLONG)
            else:
                rhill = None
            if (param['ROTATION']):
                RA = pldata[plname].ephemerides()['NPole_RA'][0]
                DEC = pldata[plname].ephemerides()['NPole_DEC'][0]

                rotpole = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree)
                rotrate = planetrot[plname] * param['TU2S']
                rot = rotpole.cartesian * rotrate
                Ip = np.array([0.0, 0.0, planetIpz[plname]])
                Ip1 = Ip[0]
                Ip2 = Ip[1]
                Ip3 = Ip[2]
                rotx = rot.x.value
                roty = rot.y.value
                rotz = rot.z.value
            else:
                Ip1 = None
                Ip2 = None
                Ip3 = None
                rotx = None
                roty = None
                rotz = None
        else:
            GMpl = None

    if idval is None:
        plid = planetid[plname]
    else:
        plid = idval

    return plname,v1,v2,v3,v4,v5,v6,idval,GMpl,Rpl,rhill,Ip1,Ip2,Ip3,rotx,roty,rotz,J2,J4

def vec2xr(param: Dict,
           namevals: npt.NDArray[np.str_],
           v1: npt.NDArray[np.float_],
           v2: npt.NDArray[np.float_],
           v3: npt.NDArray[np.float_],
           v4: npt.NDArray[np.float_],
           v5: npt.NDArray[np.float_],
           v6: npt.NDArray[np.float_],
           idvals: npt.NDArray[np.int_],
           GMpl: npt.NDArray[np.float_] | None=None,
           Rpl: npt.NDArray[np.float_] | None=None,
           rhill: npt.NDArray[np.float_] | None=None,
           Ip1: npt.NDArray[np.float_] | None=None,
           Ip2: npt.NDArray[np.float_] | None=None,
           Ip3: npt.NDArray[np.float_] | None=None,
           rotx: npt.NDArray[np.float_] | None=None,
           roty: npt.NDArray[np.float_] | None=None,
           rotz: npt.NDArray[np.float_] | None=None,
           J2: npt.NDArray[np.float_] | None=None,
           J4: npt.NDArray[np.float_] | None=None,
           t: float=0.0):
    """
    Converts and stores the variables of all bodies in an xarray dataset.

    Parameters
    ----------
    param : dict
        Swiftest paramuration parameters.
    idvals : integer 
        Array of body index values.
    namevals :

    v1 : array of floats
        xh 
    v2 : array of floats
        yh
    v3 : array of floats
        zh
    v4 : array of floats
        vhxh
    v5 : array of floats
        vhyh
    v6 : array of floats
        vhzh
    GMpl : array of floats
        G*mass
    Rpl : array of floats
        radius
    rhill : array of floats
        Hill Radius
    Ip1 : array of floats
        Principal axes moments of inertia
    Ip2 : array of floats
        Principal axes moments of inertia
    Ip3 : array of floats
        Principal axes moments of inertia
    rox : array of floats
        Rotation rate vector
    roty : array of floats
        Rotation rate vector
    rotz : array of floats
        Rotation rate vector
    t : array of floats
        Time at start of simulation
    Returns
    -------
    ds : xarray dataset
    """
    if param['ROTATION']:
        if Ip1 is None:
            Ip1 = np.full_like(v1, 0.4)
        if Ip2 is None:
            Ip2 = np.full_like(v1, 0.4)
        if Ip3 is None:
            Ip3 = np.full_like(v1, 0.4)
        if rotx is None:
            rotx = np.full_like(v1, 0.0)
        if roty is None:
            roty = np.full_like(v1, 0.0)
        if rotz is None:
            rotz = np.full_like(v1, 0.0)
    
    dims = ['time', 'id', 'vec']
    infodims = ['id', 'vec']

    # The central body is always given id 0

    if GMpl is not None:
        icb = (~np.isnan(GMpl)) & (idvals == 0)
        ipl = (~np.isnan(GMpl)) & (idvals != 0)
        itp = (np.isnan(GMpl)) & (idvals != 0)
        iscb = any(icb)
        ispl = any(ipl)
        istp = any(itp)
    else:
        icb = np.full_like(idvals,False)
        ipl = np.full_like(idvals,False)
        itp = idvals != 0
        iscb = False
        ispl = False
        istp = any(itp)

    if ispl and param['CHK_CLOSE'] and Rpl is None:
        print("Massive bodies need a radius value.")
        return None
    if ispl and rhill is None and param['RHILL_PRESENT']:
        print("rhill is required.")
        return None
   
    # Be sure we use the correct input format
    old_out_form = param['OUT_FORM']
    param['OUT_FORM'] = param['IN_FORM']
    clab, plab, tlab, infolab_float, infolab_int, infolab_str = swiftest.io.make_swiftest_labels(param)
    param['OUT_FORM'] = old_out_form
    particle_type = np.empty_like(namevals)
    vec = np.vstack([v1,v2,v3,v4,v5,v6])

    if iscb:
        lab_cb = clab.copy()
        vec_cb = np.vstack([GMpl[icb],Rpl[icb],J2[icb],J4[icb]])
        if param['ROTATION']:
            vec_cb = np.vstack([vec_cb, Ip1[icb], Ip2[icb], Ip3[icb], rotx[icb], roty[icb], rotz[icb]])
        particle_type[icb] = "Central Body"
        vec_cb = np.expand_dims(vec_cb.T,axis=0) # Make way for the time dimension!
        ds_cb = xr.DataArray(vec_cb, dims=dims, coords={'time': [t], 'id': idvals[icb], 'vec': lab_cb}).to_dataset(dim='vec')
    else:
        ds_cb =  None
    if ispl:
        lab_pl = plab.copy()
        vec_pl = np.vstack([vec[:,ipl], GMpl[ipl]])
        if param['CHK_CLOSE']:
            vec_pl = np.vstack([vec_pl, Rpl[ipl]])
        if param['RHILL_PRESENT']:
            vec_pl = np.vstack([vec_pl, rhill[ipl]])
        if param['ROTATION']:
            vec_pl = np.vstack([vec_pl, Ip1[ipl], Ip2[ipl], Ip3[ipl], rotx[ipl], roty[ipl], rotz[ipl]])
        particle_type[ipl] = np.repeat("Massive Body",idvals[ipl].size)
        vec_pl = np.expand_dims(vec_pl.T,axis=0) # Make way for the time dimension!
        ds_pl = xr.DataArray(vec_pl, dims=dims, coords={'time': [t], 'id': idvals[ipl], 'vec': lab_pl}).to_dataset(dim='vec')
    else:
        ds_pl =  None
    if istp:
        lab_tp = tlab.copy()
        vec_tp = np.expand_dims(vec[:,itp].T,axis=0) # Make way for the time dimension!
        ds_tp = xr.DataArray(vec_tp, dims=dims, coords={'time': [t], 'id': idvals[itp], 'vec': lab_tp}).to_dataset(dim='vec')
        particle_type[itp] = np.repeat("Test Particle",idvals[itp].size)
    else:
        ds_tp =  None

    ds_info = xr.DataArray(np.vstack([namevals,particle_type]).T, dims=infodims, coords={'id': idvals, 'vec' : ["name", "particle_type"]}).to_dataset(dim='vec')
    ds = [d for d in [ds_cb, ds_pl, ds_tp] if d is not None]
    if len(ds) > 1:
        ds = xr.combine_by_coords(ds)
    else:
        ds = ds[0]
    ds = xr.merge([ds_info,ds])

    return ds