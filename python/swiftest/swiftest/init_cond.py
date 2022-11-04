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

import swiftest
import numpy as np
from astroquery.jplhorizons import Horizons
import astropy.units as u
from astropy.coordinates import SkyCoord
import datetime
from datetime import date
import xarray as xr

def solar_system_horizons(plname, idval, param, ephemerides_start_date, ds):
    """
    Initializes a Swiftest dataset containing the major planets of the Solar System at a particular data from JPL/Horizons

    Parameters
    ----------
    param : dict
        Swiftest paramuration parameters. This method uses the unit conversion factors to convert from JPL's AU-day system into the system specified in the param file
    ephemerides_start_date : string
        Date to use when obtaining the ephemerides in the format YYYY-MM-DD.
    ds : xarray Dataset
        Dataset to append

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
    GMcb = np.array([swiftest.GMSunSI * param['TU2S'] ** 2 / param['DU2M'] ** 3])
    Rcb = np.array([swiftest.RSun / param['DU2M']])
    J2RP2 = np.array([swiftest.J2Sun * (swiftest.RSun / param['DU2M']) ** 2])
    J4RP4 = np.array([swiftest.J4Sun * (swiftest.RSun / param['DU2M']) ** 4])
    
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
        if param['ROTATION'] == 'YES':
            Ip1 = [Ipsun[0]]
            Ip2 = [Ipsun[1]]
            Ip3 = [Ipsun[2]]
            rotx = [rotcb.x]
            roty = [rotcb.y]
            rotz = [rotcb.z]
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

        v1 = []
        v2 = []
        v3 = []
        v4 = []
        v5 = []
        v6 = []
        J2 = None
        J4 = None

        pldata = {}
        pldata[plname] = Horizons(id=idval, location='@sun',
                               epochs={'start': ephemerides_start_date, 'stop': ephemerides_end_date,
                                       'step': ephemerides_step})
        
        if param['IN_FORM'] == 'XV':
            v1.append(pldata[plname].vectors()['x'][0] * DCONV)
            v2.append(pldata[plname].vectors()['y'][0] * DCONV)
            v3.append(pldata[plname].vectors()['z'][0] * DCONV)
            v4.append(pldata[plname].vectors()['vx'][0] * VCONV)
            v5.append(pldata[plname].vectors()['vy'][0] * VCONV)
            v6.append(pldata[plname].vectors()['vz'][0] * VCONV)
        elif param['IN_FORM'] == 'EL':
            v1.append(pldata[plname].elements()['a'][0] * DCONV)
            v2.append(pldata[plname].elements()['e'][0])
            v3.append(pldata[plname].elements()['incl'][0])
            v4.append(pldata[plname].elements()['Omega'][0])
            v5.append(pldata[plname].elements()['w'][0])
            v6.append(pldata[plname].elements()['M'][0])

        if ispl:
            GMpl = []
            GMpl.append(GMcb[0] / MSun_over_Mpl[plname])
            if param['CHK_CLOSE'] == 'YES':
                Rpl = []
                Rpl.append(planetradius[plname] * DCONV)
            else:
                Rpl = None

            # Generate planet value vectors
            if (param['RHILL_PRESENT'] == 'YES'):
                rhill = []
                rhill.append(pldata[plname].elements()['a'][0] * DCONV * (3 * MSun_over_Mpl[plname]) ** (-THIRDLONG))
            else:
                rhill = None
            if (param['ROTATION'] == 'YES'):
                Ip1 = []
                Ip2 = []
                Ip3 = []
                rotx = []
                roty = []
                rotz = []
                RA = pldata[plname].ephemerides()['NPole_RA'][0]
                DEC = pldata[plname].ephemerides()['NPole_DEC'][0]

                rotpole = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree)
                rotrate = planetrot[plname] * param['TU2S']
                rot = rotpole.cartesian * rotrate
                Ip = np.array([0.0, 0.0, planetIpz[plname]])
                Ip1.append(Ip[0])
                Ip2.append(Ip[1])
                Ip3.append(Ip[2])
                rotx.append(rot.x)
                roty.append(rot.y)
                rotz.append(rot.z)
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
        plid = np.array([planetid[plname]], dtype=int)
    else:
        plid = np.array([idval], dtype=int)

    return plid,[plname],v1,v2,v3,v4,v5,v6,GMpl,Rpl,rhill,Ip1,Ip2,Ip3,rotx,roty,rotz,J2,J4

def vec2xr(param, idvals, namevals, v1, v2, v3, v4, v5, v6, GMpl=None, Rpl=None, rhill=None, Ip1=None, Ip2=None, Ip3=None, rotx=None, roty=None, rotz=None, J2=None, J4=None,t=0.0):
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
    if v1 is None: # This is the central body
        iscb = True
    else:
        iscb = False

    if param['ROTATION'] == 'YES':
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
    if not iscb and GMpl is not None:
        ispl = True
    else:
        ispl = False
   
    if ispl and param['CHK_CLOSE'] == 'YES' and Rpl is None:
        print("Massive bodies need a radius value.")
        return None
    if ispl and rhill is None and param['RHILL_PRESENT'] == 'YES':
        print("rhill is required.")
        return None
   
    # Be sure we use the correct input format
    old_out_form = param['OUT_FORM']
    param['OUT_FORM'] = param['IN_FORM']
    clab, plab, tlab, infolab_float, infolab_int, infolab_str = swiftest.io.make_swiftest_labels(param)
    param['OUT_FORM'] = old_out_form
    vec_str = np.vstack([namevals])
    label_str = ["name"]
    if iscb:
        label_float = clab.copy()
        vec_float = np.vstack([GMpl,Rpl,J2,J4])
        if param['ROTATION'] == 'YES':
            vec_float = np.vstack([vec_float, Ip1, Ip2, Ip3, rotx, roty, rotz])
        particle_type = "Central Body"
    else:
        vec_float = np.vstack([v1, v2, v3, v4, v5, v6])
        if ispl:
            label_float = plab.copy()
            vec_float = np.vstack([vec_float, GMpl])
            if param['CHK_CLOSE'] == 'YES':
                vec_float = np.vstack([vec_float, Rpl])
            if param['RHILL_PRESENT'] == 'YES':
                vec_float = np.vstack([vec_float, rhill])
            if param['ROTATION'] == 'YES':
                vec_float = np.vstack([vec_float, Ip1, Ip2, Ip3, rotx, roty, rotz])
            particle_type = np.repeat("Massive Body",idvals.size)
        else:
            label_float = tlab.copy()
            particle_type = np.repeat("Test Particle",idvals.size)
    origin_type = np.repeat("User Added Body",idvals.size)
    origin_time = np.full_like(v1,t)
    collision_id = np.full_like(idvals,0)
    origin_xhx = v1
    origin_xhy = v2
    origin_xhz = v3
    origin_vhx = v4
    origin_vhy = v5
    origin_vhz = v6
    discard_time = np.full_like(v1,-1.0)
    status = np.repeat("ACTIVE",idvals.size)
    discard_xhx = np.zeros_like(v1)
    discard_xhy = np.zeros_like(v1)
    discard_xhz = np.zeros_like(v1)
    discard_vhx = np.zeros_like(v1)
    discard_vhy = np.zeros_like(v1)
    discard_vhz = np.zeros_like(v1)
    discard_body_id = np.full_like(idvals,-1)
    info_vec_float = np.vstack([
         origin_time,
         origin_xhx,
         origin_xhy,
         origin_xhz,
         origin_vhx,
         origin_vhy,
         origin_vhz,
         discard_time,
         discard_xhx,
         discard_xhy,
         discard_xhz,
         discard_vhx,
         discard_vhy,
         discard_vhz])
    info_vec_int = np.vstack([collision_id, discard_body_id])
    info_vec_str = np.vstack([particle_type, origin_type, status])
    frame_float = info_vec_float.T
    frame_int = info_vec_int.T
    frame_str = info_vec_str.T
    if param['IN_TYPE'] == 'NETCDF_FLOAT':
        ftype=np.float32
    elif param['IN_TYPE'] == 'NETCDF_DOUBLE' or param['IN_TYPE'] == 'ASCII':
        ftype=np.float64
    da_float = xr.DataArray(frame_float, dims=infodims, coords={'id': idvals, 'vec': infolab_float}).astype(ftype)
    da_int = xr.DataArray(frame_int, dims=infodims, coords={'id': idvals, 'vec': infolab_int})
    da_str = xr.DataArray(frame_str, dims=infodims, coords={'id': idvals, 'vec': infolab_str})
    ds_float = da_float.to_dataset(dim="vec")
    ds_int = da_int.to_dataset(dim="vec")
    ds_str = da_str.to_dataset(dim="vec")
    info_ds = xr.combine_by_coords([ds_float, ds_int, ds_str])

    frame_float = np.expand_dims(vec_float.T, axis=0)
    frame_str = vec_str.T
    da_float = xr.DataArray(frame_float, dims=dims, coords={'time': [t], 'id': idvals, 'vec': label_float}).astype(ftype)
    da_str= xr.DataArray(frame_str, dims=infodims, coords={'id': idvals, 'vec': label_str})
    ds_float = da_float.to_dataset(dim="vec")
    ds_str = da_str.to_dataset(dim="vec")
    ds = xr.combine_by_coords([ds_float, ds_str,info_ds])
    return ds