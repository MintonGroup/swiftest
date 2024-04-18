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
from . import constants
import numpy as np
from astroquery.jplhorizons import Horizons, HorizonsClass
import astropy.units as u
from astropy.coordinates import SkyCoord
import datetime
from typing import Any, Union
import warnings

def get_solar_system_body_mass_rotation(id: str,
                                        jpl: HorizonsClass=None,
                                        ephemerides_start_date: str=constants.MINTON_BCL,
                                        verbose: bool=False,
                                        **kwargs: Any) -> dict:
    """
    Parses the raw output from JPL Horizons in order to extract physical properties of a body if they exist
    

    Parameters
    ----------
    id : string or list of strings
        A string identifying which body is requested from JPL/Horizons (or a list of strings if multiple ids are possible, such as an altid list)
    jpl : HorizonsClass
        An astroquery.jplhorizons HorizonsClass object. If None, a new query will be made.
    ephemerides_start_date : string
        Date to use when obtaining the ephemerides in the format YYYY-MM-DD. Default is constants.MINTON_BCL
    verbose : bool
        Indicates whether to print messages about the query or not. Default is False
    **kwargs: Any
        Additional keyword arguments to pass to the query method (see https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html)

    Returns
    -------
    A dictionary containing the following elements
    
    Gmass : float
        G*Mass of the body in m^3/s^2
    radius : float
        The radius of the body in m
    rot : (3) float vector
        The rotation rate vector oriented toward the north pole in deg/s
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
  
    if ephemerides_start_date is None:
        ephemerides_start_date = constants.MINTON_BCL 
        
    if jpl is None:
        if type(id) != list:
            id = [id]
        jpl, altid, namelist = horizons_query(id=id[0],ephemerides_start_date=ephemerides_start_date,verbose=False,**kwargs)
    else: 
        if type(id) != list:
            altid = [id]
        else:
            altid = id

    for i in altid:
        if jpl is None:
            jpl,_,namelist = horizons_query(id=i,ephemerides_start_date=ephemerides_start_date,verbose=False,**kwargs)
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
        mass = None
    else:
        if verbose:
            print(f"Physical properties found for {namelist[0]}") 
        Gmass *= 1e9  # JPL returns GM in units of km**3 / s**2, so convert to SI
        mass = Gmass / swiftest.GC
        rotrate = get_rotrate(raw_response)
        if rotrate is None:
            rotrate = 0.0
        else:
            rotrate = np.rad2deg(rotrate)

        rot = rotpole*rotrate
        
    return {'Gmass':Gmass,'mass':mass,'radius':Rpl,'rot':rot}


def horizons_query(id: str | int, 
                   ephemerides_start_date: str, 
                   exclude_spacecraft: bool=True, 
                   verbose: bool=False,
                   **kwargs: Any) -> Union[HorizonsClass | None, list | None, list | None]:
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
    verbose: bool (optional) - Default False
        Indicate whether to print messages about the query or not

    Returns
    -------
    jpl : HorizonsClass | None
        An astroquery.jplhorizons HorizonsClass object. Or None if no match was found.
    altid : string list | None
        A list of alternate ids if more than one object matches the list
    """
    
    def get_altid(errstr,exclude_spacecraft=True):
        """
        Parses the error message returned from a failed Horizons query. If the query failed because of an ambiguous id, then it will
        return a list of id values that could possibly match the query.
        not found

        Parameters
        ----------
        errstr : string
            The error message returned from the Horizons query
        exclude_spacecraft: bool (optional) - Default True
            Indicate whether spacecraft ids should be exluded from the alternate id list        

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
            warnings.warn(f"Could not find {id} in the JPL/Horizons system",stacklevel=2)
            return None,None,None
    if verbose:
        print(f"Found matching body: {altname[0]} ({altid[0]})") 
        if len(altid) > 1:
            print("    Alternate matching bodies:")
            for i,v in enumerate(altid):
                if i > 0:
                    print(f"    {altname[i]} ({v})")
        
    return jpl,altid,altname
    
    
def get_solar_system_body(name: str,
                          ephemeris_id: str | None = None,
                          ephemerides_start_date : str = constants.MINTON_BCL,
                          central_body_name: str = "Sun",
                          verbose: bool = True,
                          **kwargs: Any) -> dict | None:
    """
    Initializes a Swiftest dataset containing the major planets of the Solar System at a particular data from JPL/Horizons

    Parameters
    ----------
    name  : string
        Name of body to add to Dataset. If `id` is not supplied, this is also what will be searched for in the JPL/Horizon's database.
        The first matching body is found (for major planets, this is the barycenter of a planet-satellite system)
    ephemeris_id : string (optional)
        If passed, this is passed to Horizons instead of `name`. This can be used to find a more precise body than given by `name`.
    ephemerides_start_date : string
        Date to use when obtaining the ephemerides in the format YYYY-MM-DD. Default is constants.MINTON_BCL
    central_body_name : string
        Name of the central body to use when calculating the relative position and velocity vectors. Default is "Sun"
    verbose : bool
        Indicates whether to print messages about the query or not. Default is True
    **kwargs: Any
        Additional keyword arguments to pass to the query method (see https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html)

    Returns
    -------
    A dictionary containing the following elements
    
    name : string
        Name of the body
    rh : (3,) array of np.float64
        Position vector array relative to the central body in m.
    vh : (3,) array of np.float64
        Velocity vector array relative to the central body in m/s.
    Gmass : np.float64
        G*mass values if these are massive bodies in m^3/s^2
    mass : np.float64
        Mass values if these are massive bodies in kg
    radius : np.float64
        Radius values if these are massive bodies in m
    rhill : np.float64 
        The Hill's radius values if these are massive bodies in m
    Ip : (3,) array of np.float64
        Principal axes moments of inertia vectors if these are massive bodies.
    rot : (3,) array of np.float
        Rotation rate vectors if these are massive bodies in deg/s
    j2rp2 : np.float64
        J_2R^2 value for the body if known
    j4rp4 : np.float64
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
    DCONV = swiftest.AU2M 
    VCONV = (swiftest.AU2M / swiftest.JD2S) 
    
    rh = np.full(3,np.nan)
    vh = np.full(3,np.nan)
    Ip = np.full(3,np.nan)
    rot = np.full(3,np.nan)
    rhill = None
    Gmass = None
    mass = None
    Rpl = None
    j2rp2 = None
    j4rp4 = None
   
    if name == "Sun" or ephemeris_id == "0": # Create central body
        if verbose:
            print("Creating the Sun as a central body")
        # Central body value vectors
        rotpoleSun = SkyCoord(ra=286.13 * u.degree, dec=63.87 * u.degree).cartesian
        rot = (360.0 / 25.05) / constants.JD2S  * rotpoleSun           
        rot = np.array([rot.x.value, rot.y.value, rot.z.value])
        Gmass = swiftest.GMSun
        Rpl = swiftest.RSun 
        rh = np.array([0.0, 0.0, 0.0])
        vh = np.array([0.0, 0.0, 0.0])
    else: # Fetch solar system ephemerides from Horizons
        if ephemeris_id is None:
            ephemeris_id = name
           
        if verbose: 
            print(f"Fetching ephemerides data for {ephemeris_id} from JPL/Horizons")
        
        jpl,altid,altname = horizons_query(ephemeris_id,ephemerides_start_date,**kwargs)
        if jpl is not None:
            if verbose:
                print(f"Found ephemerides data for {altname[0]} ({altid[0]}) from JPL/Horizons")
            if name == None:
                name = altname[0]
        else:
            return None
        
        if central_body_name != "Sun":
            jplcb, altidcb, _ = horizons_query(central_body_name,ephemerides_start_date,**kwargs)
            GMcb = get_solar_system_body_mass_rotation(altidcb,jplcb)['Gmass']
            cbrx = jplcb.vectors()['x'][0] * DCONV
            cbry = jplcb.vectors()['y'][0] * DCONV
            cbrz = jplcb.vectors()['z'][0] * DCONV
            cbvx = jplcb.vectors()['vx'][0] * VCONV
            cbvy = jplcb.vectors()['vy'][0] * VCONV
            cbvz = jplcb.vectors()['vz'][0] * VCONV
            cbrh = np.array([cbrx,cbry,cbrz])
            cbvh = np.array([cbvx,cbvy,cbvz])
        else:
            GMcb = swiftest.GMSun 
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

        Gmass,_,Rpl,rot = get_solar_system_body_mass_rotation(altid,jpl,**kwargs).values()
        # If the user inputs "Earth" or Pluto, then the Earth-Moon or Pluto-Charon barycenter and combined mass is used. 
        # To use the Earth or Pluto alone, simply pass "399" or "999", respectively to name
        if name == "Earth":
            if verbose:
                print("Combining mass of Earth and the Moon")
            Gmass_moon = get_solar_system_body_mass_rotation(["301"],**kwargs)['Gmass']
            Gmass += Gmass_moon
        elif name == "Pluto":
            if verbose:
                print("Combining mass of Pluto and Charon")
            Gmass_charon = get_solar_system_body_mass_rotation(["901"],**kwargs)['Gmass']
            Gmass += Gmass_charon 
        
        if Gmass is not None:
            rhill = jpl.elements()['a'][0] * DCONV * (Gmass / (3*GMcb))**(1.0/3.0)
            mass = Gmass / swiftest.GC
            
    if name in planetIpz:
        Ip = np.array([0.0, 0.0, planetIpz[name]])
    else:
        Ip = np.array([0.4, 0.4, 0.4])
               
    if name in planetJ2:
        j2rp2 = planetJ2[name] * Rpl**2 
        j4rp4 = planetJ4[name] * Rpl**4
        
    return {'name':name,
            'rh':rh,
            'vh':vh,
            'Gmass':Gmass,
            'mass' : mass,
            'radius': Rpl,
            'rhill': rhill,
            'Ip': Ip,
            'rot': rot,
            'j2rp2':j2rp2,
            'j4rp4':j4rp4}