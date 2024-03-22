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
from .data import SwiftestDataset
import numpy as np
import numpy.typing as ArrayLike
from astroquery.jplhorizons import Horizons
import astropy.units as u
from astropy.coordinates import SkyCoord
import datetime
import xarray as xr
from typing import (
    Dict,
    List,
    Any
)

def horizons_get_physical_properties(altid,jpl=None,**kwargs):
    """
    Parses the raw output from JPL Horizons in order to extract physical properties of a body if they exist
    

    Parameters
    ----------
    altid : list of str
        List of ids to use for Horizons query
    jpl : HorizonsClass
        An astroquery.jplhorizons HorizonsClass object. If None, a new query will be made.
    **kwargs: Any
        Additional keyword arguments to pass to the query method (see https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html)

    Returns
    -------
    Gmass : float
        G*Mass of the body
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
            # Try an alternative name for the Mass found in some satellite queries
            M = [s for s in raw_response.split('\n') if 'Mass' in s]
            if len(M) > 0:
                M = M[0].split('Mass')[-1].strip()
                if 'kg' in M:
                    unit_conv_str = M.split('kg')[0].strip()
                    unit_conv_str = unit_conv_str.split('^')[1].strip()
                    unit_conv = 10**int(unit_conv_str)
                    mult = M.split('=')[1].strip().split(' ')[1].strip('()')
                    mult = 10**int(mult.split('^')[1].strip())
                    M = M.split('=')[1].strip().split(' ')[0].strip()
                    M = float(M)  * mult * unit_conv
                    try:
                        return M * swiftest.GC * 1e-9 # Return units of km**3 / s**2 for consistency
                    except:
                        return None
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

        rotpole = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree,frame='icrs').transform_to('barycentricmeanecliptic').cartesian
         
        return np.array([rotpole.x.value, rotpole.y.value, rotpole.z.value])
    
    if type(altid) != list:
        altid = [altid]

    for id in altid:
        if jpl is None:
            jpl,_,namelist= horizons_query(id=id,ephemerides_start_date='2023-07-26',verbose=False,**kwargs)
        else:
            namelist = [jpl.table['targetname'][0]]
        raw_response = jpl.vectors_async().text
        Rpl = get_radius(raw_response) 
        if Rpl is not None:
            Rpl *= 1e3
            rotpole = get_rotpole(jpl)
            break
        else:
            jpl = None
    Gmass = get_Gmass(raw_response)
    if Rpl is None or Gmass is None:
        rot = np.full(3,np.nan) 
    else:
        print(f"Physical properties found for {namelist[0]}") 
        Gmass *= 1e9  # JPL returns GM in units of km**3 / s**2, so convert to SI

        rotrate = get_rotrate(raw_response)
        if rotrate is None:
            rotrate = 0.0
        else:
            rotrate = np.rad2deg(rotrate)

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
        altid: string list | None
            A list of alternate ids if more than one object matches the list
        altname: string list | None
            A list of alternate names if more than one object matches the list
            
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
        _=jpl.ephemerides()
        altid = [id]
        altname =[jpl.table['targetname'][0]]
    except Exception as e:
        altid,altname = get_altid(str(e))
        if altid is not None and len(altid) >0: # Return the first matching id
            id = altid[0]
            jpl = Horizons(id=id, location='@sun',
                        epochs={'start': ephemerides_start_date, 'stop': ephemerides_end_date,
                                'step': ephemerides_step})
            _=jpl.ephemerides()
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
                          central_body_name: str = "Sun",
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
    name : string
        Name of the body
    rh : (3,) array of np.float64
        Position vector array relative to the central body.
    vh : (3,) array of np.float64
        Velocity vector array relative to the central body.
    Gmass : np.float64
        G*mass values if these are massive bodies
    Rpl : np.float64
        Radius values if these are massive bodies
    rhill : np.float64 
        The Hill's radius values if these are massive bodies 
    Ip : (3,) array of np.float64
        Principal axes moments of inertia vectors if these are massive bodies.
    rot : (3,) array of np.float
        Rotation rate vectors if these are massive bodies in deg/TU
    j2r2 : np.float64
        J_2R^2 value for the body if known
    j4r4 : np.float64
        J_4R^4 value for the body if known
         
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
   
    # J2 and J4 for the major planets are from From Murray & Dermott (1999) Table A.4. 
    # The Sun is from Mecheri et al. (2004), using Corbard (b) 2002 values (Table II)
    planetJ2 = {
        'Sun' : np.longdouble(2.198e-7),
        'Mercury' : np.longdouble(60.0 * 1e-6),
        'Venus' : np.longdouble(4.0 * 1e-6),
        'Earth' : np.longdouble(1083.0 * 1e-6),
        'Mars' : np.longdouble(1960.0 * 1e-6),
        'Jupiter': np.longdouble(14736.0 * 1e-6),
        'Saturn': np.longdouble(16298.0 * 1e-6),
        'Uranus': np.longdouble(3343.0 * 1e-6),
        'Neptune': np.longdouble(3411.0 * 1e-6),
    }
    planetJ4 = { 
        'Sun' : np.longdouble(-4.805e-9),
        'Mercury' : np.longdouble(0.0),
        'Venus' : np.longdouble(2.0 * 1e-6),
        'Earth' : np.longdouble(-2.0 * 1e-6),
        'Mars' : np.longdouble(-19.0 * 1e-6),
        'Jupiter': np.longdouble(-587.0 * 1e-6),
        'Saturn': np.longdouble(-915.0 * 1e-6),
        'Uranus': np.longdouble(-29.0 * 1e-6),
        'Neptune': np.longdouble(-35.0 * 1e-6),
    }

    # Unit conversion factors
    DCONV = swiftest.AU2M / param['DU2M']
    VCONV = (swiftest.AU2M / swiftest.JD2S) / (param['DU2M'] / param['TU2S'])
    THIRDLONG = np.longdouble(1.0) / np.longdouble(3.0)
    
    param_tmp = param
    param_tmp['OUT_FORM'] = 'XVEL'

    rh = np.full(3,np.nan)
    vh = np.full(3,np.nan)
    Ip = np.full(3,np.nan)
    rot = np.full(3,np.nan)
    rhill = None
    Gmass = None
    Rpl = None
    j2r2 = None
    j4r4 = None
   

    if name == "Sun" or ephemeris_id == "0": # Create central body
        print("Creating the Sun as a central body")
        # Central body value vectors
        rotpoleSun = SkyCoord(ra=286.13 * u.degree, dec=63.87 * u.degree).cartesian
        rotSun = (360.0 / 25.05) / swiftest.JD2S  * rotpoleSun           
        rot = rotSun * param['TU2S'] 
        rot = np.array([rot.x.value, rot.y.value, rot.z.value])
        Gmass = swiftest.GMSun * param['TU2S'] ** 2 / param['DU2M'] ** 3
        Rpl = swiftest.RSun / param['DU2M']
        rh = np.array([0.0, 0.0, 0.0])
        vh = np.array([0.0, 0.0, 0.0])
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
        
        if central_body_name != "Sun":
            jplcb, altidcb, _ = horizons_query(central_body_name,ephemerides_start_date,**kwargs)
            GMcb,*_ = horizons_get_physical_properties(altidcb,jplcb)
            GMcb *= param['TU2S'] ** 2 / param['DU2M'] ** 3
            cbrx = jplcb.vectors()['x'][0] * DCONV
            cbry = jplcb.vectors()['y'][0] * DCONV
            cbrz = jplcb.vectors()['z'][0] * DCONV
            cbvx = jplcb.vectors()['vx'][0] * VCONV
            cbvy = jplcb.vectors()['vy'][0] * VCONV
            cbvz = jplcb.vectors()['vz'][0] * VCONV
            cbrh = np.array([cbrx,cbry,cbrz])
            cbvh = np.array([cbvx,cbvy,cbvz])
        else:
            GMcb = swiftest.GMSun * param['TU2S'] ** 2 / param['DU2M'] ** 3
            cbrh = np.zeros(3)
            cbvh = np.zeros(3)
    
        rx = jpl.vectors()['x'][0] * DCONV
        ry = jpl.vectors()['y'][0] * DCONV
        rz = jpl.vectors()['z'][0] * DCONV
        vx = jpl.vectors()['vx'][0] * VCONV
        vy = jpl.vectors()['vy'][0] * VCONV
        vz = jpl.vectors()['vz'][0] * VCONV

        rh = np.array([rx,ry,rz]) - cbrh
        vh = np.array([vx,vy,vz]) - cbvh

        Gmass,Rpl,rot = horizons_get_physical_properties(altid,jpl,**kwargs)
        # If the user inputs "Earth" or Pluto, then the Earth-Moon or Pluto-Charon barycenter and combined mass is used. 
        # To use the Earth or Pluto alone, simply pass "399" or "999", respectively to name
        if name == "Earth":
            print("Combining mass of Earth and the Moon")
            Gmass_moon,*_ = horizons_get_physical_properties(["301"],**kwargs)
            Gmass += Gmass_moon
        elif name == "Pluto":
            print("Combining mass of Pluto and Charon")
            Gmass_charon,*_ = horizons_get_physical_properties(["901"],**kwargs)
            Gmass += Gmass_charon 
        
        if Gmass is not None:
            # Convert from SI to system units
            Gmass *= param['TU2S'] ** 2 / param['DU2M'] ** 3
            
            Rpl /= param['DU2M']

            # Generate planet value vectors
            rhill = jpl.elements()['a'][0] * DCONV * (3 * Gmass / GMcb) ** (-THIRDLONG)

            rot *= param['TU2S']
            
    if name in planetIpz:
        Ip = np.array([0.0, 0.0, planetIpz[name]])
    else:
        Ip = np.array([0.4, 0.4, 0.4])
               
    if name in planetJ2:
        j2r2 = planetJ2[name] * Rpl**2 
        j4r4 = planetJ4[name] * Rpl**4

    return name,rh,vh,Gmass,Rpl,rhill,Ip,rot,j2r2,j4r4


def vec2xr(param: Dict, 
           name: str | ArrayLike[str],
           id : int | ArrayLike[int] | None = None,
           a : float | ArrayLike[float] | None = None,
           e : float | ArrayLike[float] | None = None,
           inc : float | ArrayLike[float] | None = None,
           capom : float | ArrayLike[float] | None = None,
           omega : float | ArrayLike[float] | None = None,
           capm : float | ArrayLike[float] | None = None,
           rh : ArrayLike[float] | None = None,
           vh : ArrayLike[float] | None = None,
           Gmass : float | ArrayLike[float] | None = None,
           radius : float | ArrayLike[float] | None = None,
           rhill : float | ArrayLike[float] | None = None,
           rot: ArrayLike[float] | None = None,
           rotphase: float | None = None,
           Ip: ArrayLike[float] | None = None,
           j2rp2: float | ArrayLike[float] | None = None,
           j4rp4: float | ArrayLike[float] | None = None,
           c_lm: ArrayLike[float] | None = None,
           time: ArrayLike[float] | None = None) -> SwiftestDataset:
    """
    Converts and stores the variables of all bodies in an xarray dataset.

    Parameters
    ----------
    param : dict
        Swiftest simulation parameters.
    name : str or array-like of str
        Name or names of bodies. Bodies are indexed by name, so these must be unique 
    id : int or array-like of int, optional
        Unique id values. 
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
        G*mass values if these are massive bodies 
    radius : float or array-like of float, optional
        Radius values if these are massive bodies
    rhill : float or array-like of float, optional
        Hill's radius values if these are massive bodies
    rot:  (n,3) array-like of float, optional
        Rotation rate vectors if these are massive bodies with rotation enabled in deg/TU
    rotphase : float
        rotational phase angle of the central body in degrees
    Ip: (n,3) array-like of flaot, optional
        Principal axes moments of inertia vectors if these are massive bodies with rotation enabled. This can be used
        instead of passing Ip1, Ip2, and Ip3 separately
    j2rp2 : float or array-like of float, optional
        J_2R^2 value for the body 
    j4rp4 : float or array-like of float, optional
        J_4R^4 value for the body 
    c_lm : (2, lmax + 1, lmax + 1) array of floats, optional
        Spherical Harmonics coefficients; lmax = max spherical harmonics order
    time : array of floats
        Time at start of simulation

    Returns
    -------
    ds : SwiftestDataset
        Dataset containing the variables of all bodies passed in kwargs
    """
    
    # Validate the inputs
    if name is None:
        raise ValueError("Name must be passed")
    
    if isinstance(name, str):
        nbody = 1
    else:
        nbody = len(name)
        
    def validate_scalars(var,nbody):
        if var is not None and len(var) != nbody:
            raise ValueError(f"{var} must be the same length as name")
        
    def validate_vectors(var,nbody):
        if var is not None and var.shape[-1] != 3:
            raise ValueError(f"{var} must have shape (n,3)")
    
    validate_scalars(id,nbody)
    validate_scalars(a,nbody)
    validate_scalars(e,nbody)
    validate_scalars(inc,nbody)
    validate_scalars(capom,nbody)
    validate_scalars(omega,nbody)
    validate_scalars(capm,nbody)
    validate_vectors(rh,nbody)
    validate_vectors(vh,nbody)
    validate_scalars(Gmass,nbody)
    validate_scalars(radius,nbody)
    validate_scalars(rhill,nbody)
    validate_vectors(rot,nbody)
    validate_vectors(Ip,nbody)
    validate_scalars(rotphase,nbody)
     
    scalar_dims = ['name']
    vector_dims = ['name','space']
    space_coords = np.array(["x","y","z"])

    vector_vars = ["rh","vh","Ip","rot"]
    scalar_vars = ["id","a","e","inc","capom","omega","capm","mass","Gmass","radius","rhill","j2rp2","j4rp4", "rotphase"]
    sph_vars = ["c_lm"]
    time_vars =  ["status","rh","vh","Ip","rot","a","e","inc","capom","omega","capm","mass","Gmass","radius","rhill","j2rp2","j4rp4", "rotphase"]
    
    if "ROTATION" in param and param['ROTATION'] == True: 
        if rot is None and Gmass is not None:
           rot = np.zeros((nbody,3))
        if Ip is None and Gmass is not None: 
            Ip = np.full((nbody,3), 0.4)

    if time is None:
        time = np.array([0.0])
        
    if param['CHK_CLOSE']:
        if Gmass is not None and radius is None: 
            raise ValueError("If Gmass is passed, then radius must also be passed when CHK_CLOSE is True")
        
    if Gmass is not None: 
        GU = swiftest.GC * param["TU2S"] ** 2 * param["MU2KG"] / param["DU2M"] ** 3
        mass = Gmass / GU
        
    valid_vars = vector_vars + scalar_vars + sph_vars + ['time','id']

    input_vars = {k:v for k,v in locals().items() if k in valid_vars and v is not None}

    data_vars = {k:(scalar_dims,v) for k,v in input_vars.items() if k in scalar_vars}
    data_vars.update({k:(vector_dims,v) for k,v in input_vars.items() if k in vector_vars})
    ds = xr.Dataset(data_vars=data_vars,
                    coords={
                        "name":(["name"],name),
                        "space":(["space"],space_coords),
                    }
                    )
    time_vars = [v for v in time_vars if v in ds]
    for v in time_vars:
        ds[v] = ds[v].expand_dims(dim={"time":1}, axis=0).assign_coords({"time": time})

    # create a C_lm Dataset and combine
    if c_lm is not None:
        clm_xr = xr.DataArray(data = c_lm,
                              coords = {
                                'sign':(['sign'], [1, -1]),
                                'l': (['l'], range(0, c_lm.shape[1])),
                                'm':(['m'], range(0, c_lm.shape[2]))
                              }
                             ).to_dataset(name='c_lm')

        ds = xr.combine_by_coords([ds, clm_xr])

    return SwiftestDataset(ds)
