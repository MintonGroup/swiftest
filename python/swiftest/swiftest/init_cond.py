import swiftest
import numpy as np
from astroquery.jplhorizons import Horizons
import datetime
import xarray as xr

def solar_system_pl(param, ephemerides_start_date):
    """
    Initializes a Swiftest dataset containing the major planets of the Solar System at a particular data from JPL/Horizons

    Parameters
    ----------
    param : dict
        Swiftest paramuration parameters. This method uses the unit conversion factors to convert from JPL's AU-day system into the system specified in the param file
    ephemerides_start_date : string
        Date to use when obtaining the ephemerides in the format YYYY-MM-DD

    Returns
    -------
    xarray dataset
    """
    # Planet ids
    planetid = {
        'mercury': '1',
        'venus': '2',
        'earthmoon': '3',
        'mars': '4',
        'jupiter': '5',
        'saturn': '6',
        'uranus': '7',
        'neptune': '8',
        'plutocharon': '9'
    }
    
    # Planet MSun/M ratio
    MSun_over_Mpl = {
        'mercury': np.longdouble(6023600.0),
        'venus': np.longdouble(408523.71),
        'earthmoon': np.longdouble(328900.56),
        'mars': np.longdouble(3098708.),
        'jupiter': np.longdouble(1047.3486),
        'saturn': np.longdouble(3497.898),
        'uranus': np.longdouble(22902.98),
        'neptune': np.longdouble(19412.24),
        'plutocharon': np.longdouble(1.35e8)
    }
    
    # Planet radii in AU
    planetradius = {
        'mercury': np.longdouble(2439.4e3) / swiftest.AU2M,
        'venus': np.longdouble(6051.8e3) / swiftest.AU2M,
        'earthmoon': np.longdouble(6371.0084e3) / swiftest.AU2M,  # Earth only for radius
        'mars': np.longdouble(3389.50e3) / swiftest.AU2M,
        'jupiter': np.longdouble(69911e3) / swiftest.AU2M,
        'saturn': np.longdouble(58232.0e3) / swiftest.AU2M,
        'uranus': np.longdouble(25362.e3) / swiftest.AU2M,
        'neptune': np.longdouble(24622.e3) / swiftest.AU2M,
        'plutocharon': np.longdouble(1188.3e3 / swiftest.AU2M)
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
    cbid = np.array([0])
    cvec = np.vstack([GMcb, Rcb, J2RP2, J4RP4])
    
    # Horizons date time internal variables
    tstart = datetime.date.fromisoformat(ephemerides_start_date)
    tstep = datetime.timedelta(days=1)
    tend = tstart + tstep
    ephemerides_end_date = tend.isoformat()
    ephemerides_step = '1d'
    
    pldata = {}
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
    Rhill = []
    Rpl = []
    GMpl = []
    
    # Fetch solar system ephemerides from Horizons
    for key, val in planetid.items():
        pldata[key] = Horizons(id=val, id_type='majorbody', location='@sun',
                               epochs={'start': ephemerides_start_date, 'stop': ephemerides_end_date,
                                       'step': ephemerides_step})
        if param['OUT_FORM'] == 'XV':
            p1.append(pldata[key].vectors()['x'][0] * DCONV)
            p2.append(pldata[key].vectors()['y'][0] * DCONV)
            p3.append(pldata[key].vectors()['z'][0] * DCONV)
            p4.append(pldata[key].vectors()['vx'][0] * VCONV)
            p5.append(pldata[key].vectors()['vy'][0] * VCONV)
            p6.append(pldata[key].vectors()['vz'][0] * VCONV)
            p7.append(pldata[key].elements()['a'][0] * DCONV)
            p8.append(pldata[key].elements()['e'][0])
            p9.append(pldata[key].elements()['incl'][0] * np.pi / 180.0)
            p10.append(pldata[key].elements()['Omega'][0] * np.pi / 180.0)
            p11.append(pldata[key].elements()['w'][0] * np.pi / 180.0)
            p12.append(pldata[key].elements()['M'][0] * np.pi / 180.0)
        elif param['OUT_FORM'] == 'EL':
            p1.append(pldata[key].elements()['a'][0] * DCONV)
            p2.append(pldata[key].elements()['e'][0])
            p3.append(pldata[key].elements()['incl'][0] * np.pi / 180.0)
            p4.append(pldata[key].elements()['Omega'][0] * np.pi / 180.0)
            p5.append(pldata[key].elements()['w'][0] * np.pi / 180.0)
            p6.append(pldata[key].elements()['M'][0] * np.pi / 180.0)
            p7.append(pldata[key].vectors()['x'][0] * DCONV)
            p8.append(pldata[key].vectors()['y'][0] * DCONV)
            p9.append(pldata[key].vectors()['z'][0] * DCONV)
            p10.append(pldata[key].vectors()['vx'][0] * VCONV)
            p11.append(pldata[key].vectors()['vy'][0] * VCONV)
            p12.append(pldata[key].vectors()['vz'][0] * VCONV)
        Rhill.append(pldata[key].elements()['a'][0] * (3 * MSun_over_Mpl[key]) ** (-THIRDLONG))
        Rpl.append(planetradius[key] * DCONV)
        GMpl.append(GMcb[0] / MSun_over_Mpl[key])
    # Generate planet value vectors
    plid = np.fromiter(planetid.values(), dtype=int)
    pvec = np.vstack([p1, p2, p3, p4, p5, p6, GMpl, Rpl, p7, p8, p9, p10, p11, p12, Rhill])
    
    dims = ['time', 'id', 'vec']
    cb = []
    pl = []
    t = np.array([0.0])
    
    clab, plab, tlab = swiftest.io.make_swiftest_labels(param)
    if param['OUT_FORM'] == 'XV':
        plab.append('a')
        plab.append('e')
        plab.append('inc')
        plab.append('capom')
        plab.append('omega')
        plab.append('capm')
    elif param['OUT_FORM'] == 'EL':
        plab.append('px')
        plab.append('py')
        plab.append('pz')
        plab.append('vx')
        plab.append('vy')
        plab.append('vz')
    plab.append('Rhill')
    
    # Prepare frames by adding an extra axis for the time coordinate
    cbframe = np.expand_dims(cvec.T, axis=0)
    plframe = np.expand_dims(pvec.T, axis=0)
    
    # Create xarray DataArrays out of each body type
    cbxr = xr.DataArray(cbframe, dims=dims, coords={'time': t, 'id': cbid, 'vec': clab})
    plxr = xr.DataArray(plframe, dims=dims, coords={'time': t, 'id': plid, 'vec': plab})
    
    cb.append(cbxr)
    pl.append(plxr)
    
    cbda = xr.concat(cb, dim='time')
    plda = xr.concat(pl, dim='time')
    
    cbds = cbda.to_dataset(dim='vec')
    plds = plda.to_dataset(dim='vec')
    ds = xr.combine_by_coords([cbds, plds])
    return ds