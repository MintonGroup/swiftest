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
    xarray dataset
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
    
    cbid = np.array([0])
    cvec = np.vstack([GMcb, Rcb, J2RP2, J4RP4])
    if param['ROTATION'] == 'YES':
        cvec = np.vstack([cvec, Ipsun[0], Ipsun[1], Ipsun[2], rotcb.x, rotcb.y, rotcb.z])
    
    # Horizons date time internal variables
    tstart = datetime.date.fromisoformat(ephemerides_start_date)
    tstep = datetime.timedelta(days=1)
    tend = tstart + tstep
    ephemerides_end_date = tend.isoformat()
    ephemerides_step = '1d'
   
    param_tmp = param
    param_tmp['OUT_FORM'] = 'XVEL'
    clab, plab, tlab = swiftest.io.make_swiftest_labels(param)

    dims = ['time', 'id', 'vec']
    infodims = ['id', 'vec']
    t = np.array([0.0])

    key = plname
    if ispl:
        val = planetid[key]
    else:
        val = -1
    if key == "Sun" : # Create central body
        print("Creating the Sun as a central body")
        cbframe = np.expand_dims(cvec.T, axis=0)
        cbxr = xr.DataArray(cbframe, dims=dims, coords={'time': t, 'id': cbid, 'vec': clab})
        cbxr = cbxr.assign_coords(name=('id', [key]))
        cbds = cbxr.to_dataset(dim='vec')
        ds = xr.combine_by_coords([ds, cbds])
    else: # Fetch solar system ephemerides from Horizons
        print(f"Fetching ephemerides data for {key} from JPL/Horizons")

        p1 = []
        p2 = []
        p3 = []
        p4 = []
        p5 = []
        p6 = []
        p7 = []
        p8 = []
        p9 = []
        p10 = []
        p11 = []
        p12 = []
        rhill = []
        Rpl = []
        GMpl = []
        Ip1 = []
        Ip2 = []
        Ip3 = []
        rotx = []
        roty = []
        rotz = []

        pldata = {}
        if ispl:
            pldata[key] = Horizons(id=val, id_type='majorbody', location='@sun',
                                   epochs={'start': ephemerides_start_date, 'stop': ephemerides_end_date,
                                           'step': ephemerides_step})
        else:
            pldata[key] = Horizons(id=key, id_type='smallbody', location='@sun',
                                   epochs={'start': ephemerides_start_date, 'stop': ephemerides_end_date,
                                           'step': ephemerides_step})
        
        if (param['OUT_FORM'] == 'XV' or param['OUT_FORM'] == 'XVEL'):
            p1.append(pldata[key].vectors()['x'][0] * DCONV)
            p2.append(pldata[key].vectors()['y'][0] * DCONV)
            p3.append(pldata[key].vectors()['z'][0] * DCONV)
            p4.append(pldata[key].vectors()['vx'][0] * VCONV)
            p5.append(pldata[key].vectors()['vy'][0] * VCONV)
            p6.append(pldata[key].vectors()['vz'][0] * VCONV)
            p7.append(pldata[key].elements()['a'][0] * DCONV)
            p8.append(pldata[key].elements()['e'][0])
            p9.append(pldata[key].elements()['incl'][0])
            p10.append(pldata[key].elements()['Omega'][0])
            p11.append(pldata[key].elements()['w'][0])
            p12.append(pldata[key].elements()['M'][0])
        elif param['OUT_FORM'] == 'EL':
            p1.append(pldata[key].elements()['a'][0] * DCONV)
            p2.append(pldata[key].elements()['e'][0])
            p3.append(pldata[key].elements()['incl'][0])
            p4.append(pldata[key].elements()['Omega'][0])
            p5.append(pldata[key].elements()['w'][0])
            p6.append(pldata[key].elements()['M'][0])
            p7.append(pldata[key].vectors()['x'][0] * DCONV)
            p8.append(pldata[key].vectors()['y'][0] * DCONV)
            p9.append(pldata[key].vectors()['z'][0] * DCONV)
            p10.append(pldata[key].vectors()['vx'][0] * VCONV)
            p11.append(pldata[key].vectors()['vy'][0] * VCONV)
            p12.append(pldata[key].vectors()['vz'][0] * VCONV)
        pvec = np.vstack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12])

        if ispl:
            Rpl.append(planetradius[key] * DCONV)
            GMpl.append(GMcb[0] / MSun_over_Mpl[key])
            pvec = np.vstack([pvec, GMpl, Rpl])

            # Generate planet value vectors
            if (param['RHILL_PRESENT'] == 'YES'):
                rhill.append(pldata[key].elements()['a'][0] * DCONV * (3 * MSun_over_Mpl[key]) ** (-THIRDLONG))
                pvec = np.vstack([pvec, rhill])
            if (param['ROTATION'] == 'YES'):
                RA = pldata[key].ephemerides()['NPole_RA'][0]
                DEC = pldata[key].ephemerides()['NPole_DEC'][0]

                rotpole = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree)
                rotrate = planetrot[key] * param['TU2S']
                rot = rotpole.cartesian * rotrate
                Ip = np.array([0.0, 0.0, planetIpz[key]])
                Ip1.append(Ip[0])
                Ip2.append(Ip[1])
                Ip3.append(Ip[2])
                rotx.append(rot.x)
                roty.append(rot.y)
                rotz.append(rot.z)
                pvec = np.vstack([pvec, Ip1, Ip2, Ip3, rotx, roty, rotz])
        else:
            plab = tlab.copy()

        if idval is None:
            plid = np.array([planetid[key]], dtype=int)
        else:
            plid = np.array([idval], dtype=int)

        # Prepare frames by adding an extra axis for the time coordinate
        plframe = np.expand_dims(pvec.T, axis=0)
        plxr = xr.DataArray(plframe, dims=dims, coords={'time': t, 'id': plid, 'vec': plab})
        plxr = plxr.assign_coords(name=('id', [plname]))
        plds = plxr.to_dataset(dim='vec')
        ds = xr.combine_by_coords([ds, plds])

    return ds

def vec2xr(param, idvals, namevals, v1, v2, v3, v4, v5, v6, GMpl=None, Rpl=None, rhill=None, Ip1=None, Ip2=None, Ip3=None, rotx=None, roty=None, rotz=None, t=0.0):
    
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
    if GMpl is not None:
        ispl = True
    else:
        ispl = False
   
    if ispl and Rpl is None:
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
    vec_float = np.vstack([v1, v2, v3, v4, v5, v6])
    vec_str = np.vstack([namevals])
    label_str = ["name"]
    if ispl:
        label_float = plab.copy()
        vec_float = np.vstack([vec_float, GMpl, Rpl])
        if param['RHILL_PRESENT'] == 'YES':
            vec_float = np.vstack([vec_float, rhill])
        if param['ROTATION'] == 'YES':
            vec_float = np.vstack([vec_float, Ip1, Ip2, Ip3, rotx, roty, rotz])
    else:
        label_float = tlab.copy()
    if param['OUT_STAT'] == 'APPEND': # Include particle information in initial conditions when appending an existing run
        if ispl:
            particle_type = np.repeat("Massive Body",idvals.size)
        else:
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
        da_float = xr.DataArray(frame_float, dims=infodims, coords={'id': idvals, 'vec': infolab_float})
        da_int = xr.DataArray(frame_int, dims=infodims, coords={'id': idvals, 'vec': infolab_int})
        da_str = xr.DataArray(frame_str, dims=infodims, coords={'id': idvals, 'vec': infolab_str})
        ds_float = da_float.to_dataset(dim="vec")
        ds_int = da_int.to_dataset(dim="vec")
        ds_str = da_str.to_dataset(dim="vec")
        info_ds = xr.combine_by_coords([ds_float, ds_int, ds_str])

    frame_float = np.expand_dims(vec_float.T, axis=0)
    frame_str = vec_str.T
    da_float = xr.DataArray(frame_float, dims=dims, coords={'time': [t], 'id': idvals, 'vec': label_float})
    da_str= xr.DataArray(frame_str, dims=infodims, coords={'id': idvals, 'vec': label_str})
    ds_float = da_float.to_dataset(dim="vec")
    ds_str = da_str.to_dataset(dim="vec")
    ds = xr.combine_by_coords([ds_float, ds_str])
    if param['OUT_STAT'] == 'APPEND':
        ds = xr.combine_by_coords([ds,info_ds])
    return ds