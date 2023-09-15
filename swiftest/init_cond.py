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

def horizons_get_physical_properties(altid,**kwargs):
    """
    Parses the raw output from JPL Horizons in order to extract physical properties of a body if they exist
    

    Parameters
    ----------
    altid : list of str
        List of ids to use for Horizons query
    **kwargs: Any
            Additional keyword arguments to pass to the query method (see https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html)


    Returns
    -------
    MSun_over_Mpl : float
        The ratio of MSun/M of the body
    radius : float
        The radius of the body in m
    rot: (3) float vector
        The rotation rate vector oriented toward the north pole
    """

    def get_Gmass(raw_response):
        GM = [s for s in raw_response.split('\n') if 'GM' in s 
            and 'Rel' not in s 
            and 'GMT' not in s
            and 'ANGMOM' not in s]
        if len(GM) == 0:
            return None
        GM = GM[0]
        if len(GM) > 1:
            GM = GM.split('GM')[1]
            if len(GM) > 1:
                GM = GM.split('=')
                if len(GM) > 1:
                    GM = GM[1].strip().split(' ')[0].split('+')[0]
        try:
            GM = float(GM)
        except:
            GM = None
        return GM

    def get_radius(raw_response):

        radius = [s for s in raw_response.split('\n') if 'mean radius' in s.lower() or 'RAD' in s or 'Radius (km)' in s]
        if len(radius) == 0:
            return None
        radius = radius[0]
        if "Radius (km)" in radius:
            radius = radius.split("Radius (km)")[1].strip(' =').split('+')[0].split()[0].strip()
        elif "R_eff" in radius: # Some small bodies list an effective radius 
            radius = radius.split('R_eff')[1].strip(' =').split('+')[0].split()[0].strip()
        elif "RAD" in radius: # Some small bodies list the radius like this
            radius = radius.split('RAD')[1].strip(' =').split('+')[0].split()[0].strip()
            if 'x' in radius: # Triaxial ellipsoid bodies like Haumea may have multiple dimensions which need to be averaged
                radius = radius.split('x')
                try:
                    for i,v in enumerate(radius):
                        radius[i] = float(v.split()[0].strip())
                    radius = np.average(radius)
                except:
                    radius = None
        else: # Handles most major bodies
            radius = radius.split('=')[1].strip().split('+')[0].split()[0].strip() 
        try:
            radius = float(radius)
        except:
            radius = None
        return radius

    def get_rotrate(raw_response):
        raw_response=jpl.raw_response
        rotrate = [s for s in raw_response.split('\n') if 'rot. rat' in s.lower()]
        if len(rotrate) == 0:
            rotrate = [s for s in raw_response.split('\n') if 'ROTPER' in s.upper()] # Try the small body version
            if len(rotrate) > 0:
                rotrate = rotrate[0].split('ROTPER=')[1].strip()
                try:
                   rotrate = 2*np.pi / (float(rotrate) * 3600)
                except:
                   rotrate = None
            else:
                if "Synchronous" in raw_response: # Satellites have this:
                    rotrate = [s for s in raw_response.split('\n') if 'Orbital period' in s][0]
                    rotrate = rotrate.split('Orbital period')[1].replace('~',' ').replace('d',' ').replace('=',' ').strip()
                    rotrate = 2*np.pi / (float(rotrate) *  swiftest.JD2S) 
                else:   
                    rotrate = None
        else:
            rotrate = rotrate[0].lower().split('rot. rat')[1].split('=')[1].strip().split('  ')[0].strip()
            try:
                rotrate = float(rotrate)
            except:
                rotrate = None
        return rotrate

    def get_rotpole(jpl):
        RA = jpl.ephemerides()['NPole_RA'][0]
        DEC = jpl.ephemerides()['NPole_DEC'][0]

        if np.ma.is_masked(RA) or np.ma.is_masked(DEC):
            return np.array([0.0,0.0,1.0])

        rotpole = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree).cartesian
        return np.array([rotpole.x.value, rotpole.y.value, rotpole.z.value])
    
    if type(altid) != list:
        altid = [altid]

    for id in altid:
        jpl,idlist,namelist = horizons_query(id=id,ephemerides_start_date='2023-07-26',verbose=False,**kwargs)
        Rpl = get_radius(jpl.raw_response) 
        if Rpl is not None:
            Rpl *= 1e3
            break
        
    Gmass = get_Gmass(jpl.raw_response)
    if Rpl is None or Gmass is None:
        rot = np.full(3,np.nan) 
    else:
        print(f"Physical properties found for {namelist[0]}") 
        Gmass *= 1e9  # JPL returns GM in units of km**3 / s**2, so convert to SI

        rotrate = get_rotrate(jpl.raw_response)
        if rotrate is None:
            rotrate = 0.0
        else:
            rotrate = np.rad2deg(rotrate)

        rotpole = get_rotpole(jpl)
        rot = rotpole*rotrate
            
    return Gmass,Rpl,rot


def horizons_query(id, ephemerides_start_date, exclude_spacecraft=True, verbose=False,**kwargs):
    """
    Queries JPL/Horizons for a body matching the id. If one is found, a HorizonsClass object is returned for the first object that
    matches the passed id string. If more than one match is found, a list of alternate ids is also returned. If no object is found
    then None is returned.

    Parameters
    ----------
    id : string
        A string identifying which body is requested from JPL/Horizons
    ephemerides_start_date : string
        Date to use when obtaining the ephemerides in the format YYYY-MM-DD.
    exclude_spacecraft: bool (optional) - Default True
        Indicate whether spacecraft ids should be exluded from the alternate id list
    verbose: bool (optional) - Default True
        Indicate whether to print messages about the query or not

    Returns
    -------
    jpl: HorizonsClass | None
        An astroquery.jplhorizons HorizonsClass object. Or None if no match was found.
    altid: string list | None
        A list of alternate ids if more than one object matches the list
        
    """
    
    def get_altid(errstr,exclude_spacecraft=True):
        """
        Parses the error message returned from a failed Horizons query. If the query failed because of an ambiguous id, then it will
        return a list of id values that could possibly match the query.
        not found

        Parameters
        ----------
        raw_response : string
            Raw response from the JPL Horizons query

        Returns
        -------
        MSun_over_Mpl : float
        """    
        if "ID" in errstr:
            altid = errstr.split('ID')[1]
            altid = altid.split('\n')[2:-1]
            altname = []
            for n,v in enumerate(altid):
                altid[n] = v.strip().split(' ')[0]
                altname.append(' '.join(v.strip().split(' ')[1:]).strip().split('  ')[0])
            if exclude_spacecraft:
                altname = [altname[i] for i,n in enumerate(altid) if int(n) > 0] 
                altid = [n for n in altid if int(n) > 0]
                
            return altid,altname 
        else:
            return None,None
        
        
    # Horizons date time internal variables
    tstart = datetime.date.fromisoformat(ephemerides_start_date)
    tstep = datetime.timedelta(days=1)
    tend = tstart + tstep
    ephemerides_end_date = tend.isoformat()
    ephemerides_step = '1d'
    
    try:
        jpl = Horizons(id=id, location='@sun',
                            epochs={'start': ephemerides_start_date, 'stop': ephemerides_end_date,
                                    'step': ephemerides_step},**kwargs)
        eph=jpl.ephemerides()
        altid = [id]
        altname =[jpl.table['targetname'][0]]
    except Exception as e:
        altid,altname = get_altid(str(e))
        if altid is not None and len(altid) >0: # Return the first matching id
            id = altid[0]
            jpl = Horizons(id=id, location='@sun',
                        epochs={'start': ephemerides_start_date, 'stop': ephemerides_end_date,
                                'step': ephemerides_step})
            eph=jpl.ephemerides()
        else:
            print(f"Could not find {id} in the JPL/Horizons system")
            return None,None,None
    if verbose:
        print(f"Found matching body: {altname[0]} ({altid[0]})") 
        if len(altid) > 1:
            print("    Alternate matching bodies:")
            for i,v in enumerate(altid):
                if i > 0:
                    print(f"    {altname[i]} ({v})")
        
    return jpl,altid,altname
    
    
def solar_system_horizons(name: str,
                          param: Dict,
                          ephemerides_start_date: str,
                          ephemeris_id: str | None = None,
                          **kwargs: Any):
    """
    Initializes a Swiftest dataset containing the major planets of the Solar System at a particular data from JPL/Horizons

    Parameters
    ----------
    name  : string
        Name of body to add to Dataset. If `id` is not supplied, this is also what will be searche for in the JPL/Horizon's database.
        The first matching body is found (for major planets, this is the barycenter of a planet-satellite system)
    param : dict
        Swiftest paramuration parameters. This method uses the unit conversion factors to convert from JPL's AU-day system into the system specified in the param file
    ephemerides_start_date : string
        Date to use when obtaining the ephemerides in the format YYYY-MM-DD.
    ephemeris_id : string (optional)
        If passed, this is passed to Horizons instead of `name`. This can be used to find a more precise body than given by `name`. 
    **kwargs: Any
            Additional keyword arguments to pass to the query method (see https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html)

    Returns
    -------
    ds : xarray dataset
        Initial conditions of body formatted for Swiftest
    
    Notes
    --------
    When passing `name` == "Earth" or `name` == "Pluto", it a body is generated that has initial conditions matching the system
    barycenter and mass equal to the sum of Earth+Moon or Pluto+Charon. To obtain initial conditions for either Earth or Pluto alone,
    pass `ephemeris_id` == "399" for Earth or `id` == "999" for Pluto. 
    """
    
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
    
    rotcb = swiftest.rotSun * param['TU2S'] 
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

    if name == "Sun" or ephemeris_id == "0": # Create central body
        print("Creating the Sun as a central body")
        Gmass = GMcb
        Rpl = Rcb
        J2 = J2RP2
        J4 = J4RP4
        if param['ROTATION']:
            Ip = Ipsun
            rot = rotcb
    else: # Fetch solar system ephemerides from Horizons
        if ephemeris_id is None:
            ephemeris_id = name
            
        print(f"Fetching ephemerides data for {ephemeris_id} from JPL/Horizons")
        
        jpl,altid,altname = horizons_query(ephemeris_id,ephemerides_start_date,**kwargs)
        if jpl is not None:
            print(f"Found ephemerides data for {altname[0]} ({altid[0]}) from JPL/Horizons")
            if name == None:
                name = altname[0]
        else:
            return None
        
        if param['IN_FORM'] == 'XV':
            rx = jpl.vectors()['x'][0] * DCONV
            ry = jpl.vectors()['y'][0] * DCONV
            rz = jpl.vectors()['z'][0] * DCONV
            vx = jpl.vectors()['vx'][0] * VCONV
            vy = jpl.vectors()['vy'][0] * VCONV
            vz = jpl.vectors()['vz'][0] * VCONV

            rh = np.array([rx,ry,rz])
            vh = np.array([vx,vy,vz])
        elif param['IN_FORM'] == 'EL':
            a = jpl.elements()['a'][0] * DCONV
            e = jpl.elements()['e'][0]
            inc = jpl.elements()['incl'][0]
            capom = jpl.elements()['Omega'][0]
            omega = jpl.elements()['w'][0]
            capm = jpl.elements()['M'][0]

        Gmass,Rpl,rot = horizons_get_physical_properties(altid,**kwargs)
        # If the user inputs "Earth" or Pluto, then the Earth-Moon or Pluto-Charon barycenter and combined mass is used. 
        # To use the Earth or Pluto alone, simply pass "399" or "999", respectively to name
        if ephemeris_id == "Earth":
            print("Combining mass of Earth and the Moon")
            Gmass_moon,tmp,tmp = horizons_get_physical_properties(["301"],**kwargs)
            Gmass += Gmass_moon
        elif ephemeris_id == "Pluto":
            print("Combining mass of Pluto and Charon")
            Gmass_charon,tmp,tmp = horizons_get_physical_properties(["901"],**kwargs)
            Gmass += Gmass_charon 
        
        if Gmass is not None:
            # Convert from SI to system units
            Gmass *= param['TU2S'] ** 2 / param['DU2M'] ** 3
            
            if param['CHK_CLOSE']:
                Rpl /= param['DU2M']

            # Generate planet value vectors
            if (param['RHILL_PRESENT']):
                rhill = jpl.elements()['a'][0] * DCONV * (3 * Gmass / GMcb) ** (-THIRDLONG)

            if (param['ROTATION']):
                rot *= param['TU2S']
                if name in planetIpz:
                    Ip = np.array([0.0, 0.0, planetIpz[name]])
                else:
                    Ip = np.array([0.4, 0.4, 0.4])
        else:
            Gmass = None

    # Only the Sun gets assigned its own special id for now. All other ids will be sorted later
    if name == "Sun":
        id = 0
    else:
        id = -1

    return id,name,a,e,inc,capom,omega,capm,rh,vh,Gmass,Rpl,rhill,Ip,rot,J2,J4

def vec2xr(param: Dict, **kwargs: Any):
    """
    Converts and stores the variables of all bodies in an xarray dataset.

    Parameters
    ----------
    param : dict
        Swiftest simulation parameters.
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
        Position vector array.
    vh : (n,3) array-like of float, optional
        Velocity vector array. 
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
    scalar_vars = ["name","a","e","inc","capom","omega","capm","Gmass","radius","rhill","j2rp2","j4rp4"]
    time_vars =  ["rh","vh","Ip","rot","a","e","inc","capom","omega","capm","Gmass","radius","rhill","j2rp2","j4rp4"]

    # Check for valid keyword arguments
    kwargs = {k:kwargs[k] for k,v in kwargs.items() if v is not None}
  
    if "ROTATION" in param and param['ROTATION'] == True: 
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