import numpy as np
from scipy.io import FortranFile
import xarray as xr
from astroquery.jplhorizons import Horizons
import astropy.constants as const
import datetime

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


#I/O Routines for reading in Swifter and Swiftest parameter and binary data files
def read_swifter_param(inparfile):
    """
    Reads in a Swifter param.in file and saves it as a dictionary

    Parameters
    ----------
    inparfile : string
        File name of the input parameter file

    Returns
    -------
    param
        A dictionary containing the entries in the user parameter file
    """
    param = {
    'INPARFILE'      : inparfile,
    'NPLMAX'         : -1,
    'NTPMAX'         : -1,
    'T0'             : 0.0,
    'TSTOP'          : 0.0,
    'DT'             : 0.0,
    'PL_IN'          : "",
    'TP_IN'          : "",
    'IN_TYPE'        : "ASCII",
    'ISTEP_OUT'      : -1,
    'BIN_OUT'        : "",
    'OUT_TYPE'       : 'REAL8',
    'OUT_FORM'       : "XV",
    'OUT_STAT'       : "NEW",
    'ISTEP_DUMP'     : -1,
    'J2'             : 0.0,
    'J4'             : 0.0,
    'CHK_CLOSE'      : 'NO',
    'CHK_RMIN'       : -1.0,
    'CHK_RMAX'       : -1.0,
    'CHK_EJECT'      : -1.0,
    'CHK_QMIN'       : -1.0,
    'CHK_QMIN_COORD' : "HELIO",
    'CHK_QMIN_RANGE' : "",
    'QMIN_ALO'       : -1.0,
    'QMIN_AHI'       : -1.0,
    'ENC_OUT'        : "",
    'EXTRA_FORCE'    : 'NO',
    'BIG_DISCARD'    : 'NO',
    'RHILL_PRESENT'  : 'NO',
    'GR'             : 'NO',
    'C2'             : -1.0,
             }

    # Read param.in file
    print(f'Reading Swifter file {inparfile}')
    f = open(inparfile, 'r')
    swifterlines = f.readlines()
    f.close()
    for line in swifterlines:
        fields = line.split()
        if len(fields) > 0:
            for key in param:
                if (key == fields[0].upper()): param[key] = fields[1]
            #Special case of CHK_QMIN_RANGE requires a second input
            if (param['CHK_QMIN_RANGE'] == fields[0].upper()):
                param['QMIN_ALO'] = fields[1]
                param['QMIN_AHI'] = fields[2]

    param['NPLMAX']     = int(param['NPLMAX'])
    param['NTPMAX']     = int(param['NTPMAX'])
    param['ISTEP_OUT']  = int(param['ISTEP_OUT'])
    param['ISTEP_DUMP'] = int(param['ISTEP_DUMP'])
    param['T0']         = float(param['T0'])
    param['TSTOP']      = float(param['TSTOP'])
    param['DT']         = float(param['DT'])
    param['J2']         = float(param['J2'])
    param['J4']         = float(param['J4'])
    param['CHK_RMIN']   = float(param['CHK_RMIN'])
    param['CHK_RMAX']   = float(param['CHK_RMAX'])
    param['CHK_EJECT']  = float(param['CHK_EJECT'])
    param['CHK_QMIN']   = float(param['CHK_QMIN'])
    param['QMIN_ALO']   = float(param['QMIN_ALO'])
    param['QMIN_AHI']   = float(param['QMIN_AHI'])
    param['C2']         = float(param['C2'])
    param['INV_C2']     = param['C2']
    param['EXTRA_FORCE'] = param['EXTRA_FORCE'].upper()
    param['BIG_DISCARD'] = param['BIG_DISCARD'].upper()
    param['CHK_CLOSE']  = param['CHK_CLOSE'].upper()
    param['RHILL_PRESENT'] = param['RHILL_PRESENT'].upper()
    param['GR']         = param['GR'].upper()

    return param

def read_swiftest_param(param_file_name):
    """
    Reads in a Swiftest param.in file and saves it as a dictionary

    Parameters
    ----------
    param_file_name : string
        File name of the input parameter file

    Returns
    -------
    param : dict
        A dictionary containing the entries in the user parameter file
    """
    param = {
    'param_FILE_NAME' : param_file_name,
    'NPLMAX'         : -1,
    'NTPMAX'         : -1,
    'T0'             : 0.0,
    'TSTOP'          : 0.0,
    'DT'             : 0.0,
    'PL_IN'          : "",
    'TP_IN'          : "",
    'CB_IN'          : "",
    'IN_TYPE'        : "ASCII",
    'ISTEP_OUT'      : -1,
    'BIN_OUT'        : "",
    'OUT_TYPE'       : 'REAL8',
    'OUT_FORM'       : "XV",
    'OUT_STAT'       : "NEW",
    'ISTEP_DUMP'     : -1,
    'J2'             : 0.0,
    'J4'             : 0.0,
    'CHK_RMIN'       : -1.0,
    'CHK_RMAX'       : -1.0,
    'CHK_EJECT'      : -1.0,
    'CHK_QMIN'       : -1.0,
    'CHK_QMIN_COORD' : "HELIO",
    'CHK_QMIN_RANGE' : "",
    'QMIN_ALO'       : -1.0,
    'QMIN_AHI'       : -1.0,
    'ENC_OUT'        : "",
    'MTINY'          : -1.0,
    'MU2KG'          : -1.0,
    'TU2S'           : -1.0,
    'DU2M'           : -1.0,
    'GU'             : -1.0,
    'EXTRA_FORCE'    : 'NO',
    'BIG_DISCARD'    : 'NO',
    'CHK_CLOSE'      : 'NO',
    'FRAGMENTATION'  : 'NO',
    'MTINY_SET'      : 'NO',
    'ROTATION'       : 'NO',
    'TIDES'          : 'NO',
    'ENERGY'         : 'NO',
    'GR'             : 'NO',
    'YARKOVSKY'      : 'NO',
    'YORP'           : 'NO',
    }

    # Read param.in file
    print(f'Reading Swiftest file {param_file_name}' )
    f = open(param_file_name, 'r')
    swiftestlines = f.readlines()
    f.close()
    for line in swiftestlines:
        fields = line.split()
        if len(fields) > 0:
            for key in param:
                if (key == fields[0].upper()): param[key] = fields[1]
            #Special case of CHK_QMIN_RANGE requires a second input
            if (param['CHK_QMIN_RANGE'] == fields[0].upper()):
                param['QMIN_ALO'] = fields[1]
                param['QMIN_AHI'] = fields[2]

    param['NPLMAX']     = int(param['NPLMAX'])
    param['NTPMAX']     = int(param['NTPMAX'])
    param['ISTEP_OUT']  = int(param['ISTEP_OUT'])
    param['ISTEP_DUMP'] = int(param['ISTEP_DUMP'])
    param['T0']         = float(param['T0'])
    param['TSTOP']      = float(param['TSTOP'])
    param['DT']         = float(param['DT'])
    param['J2']         = float(param['J2'])
    param['J4']         = float(param['J4'])
    param['CHK_RMIN']   = float(param['CHK_RMIN'])
    param['CHK_RMAX']   = float(param['CHK_RMAX'])
    param['CHK_EJECT']  = float(param['CHK_EJECT'])
    param['CHK_QMIN']   = float(param['CHK_QMIN'])
    param['QMIN_ALO']   = float(param['QMIN_ALO'])
    param['QMIN_AHI']   = float(param['QMIN_AHI'])
    param['MTINY']      = float(param['MTINY'])
    param['DU2M']       = float(param['DU2M'])
    param['MU2KG']      = float(param['MU2KG'])
    param['TU2S']       = float(param['TU2S'])
    param['EXTRA_FORCE'] = param['EXTRA_FORCE'].upper()
    param['BIG_DISCARD'] = param['BIG_DISCARD'].upper()
    param['CHK_CLOSE']   = param['CHK_CLOSE'].upper()
    param['FRAGMENTATION'] = param['FRAGMENTATION'].upper()
    param['ROTATION']    = param['ROTATION'].upper()
    param['TIDES']       = param['TIDES'].upper()
    param['ENERGY']      = param['ENERGY'].upper()
    param['GR']          = param['GR'].upper()
    param['YORP']        = param['YORP'].upper()

    param['GU']         = GC / (param['DU2M']**3 / (param['MU2KG'] * param['TU2S']**2))
    param['INV_C2']     = einsteinC * param['TU2S'] / param['DU2M']
    param['INV_C2']     = param['INV_C2']**(-2)
    return param

def swifter_stream(f, param):
    """
    Reads in a Swifter bin.dat file and returns a single frame of data as a datastream

    Parameters
    ----------
    f : file object
    param : dict

    Yields
    -------
    t    : float
        Time of this frame
    npl  : int
        Number of massive bodies
    plid : int array
        IDs of massive bodies
    pvec : float array
        (npl,N) - vector of N quantities or each particle (6 of XV/EL + Mass, Radius, etc)
    plab : string list
        Labels for the pvec data
    ntp  : int
        Number of test particles
    tpid : int array
        Ids of test particles
    tvec : float array
        (ntp,N) - vector of N quantities for each particle (6 of XV/EL, etc.)
    tlab : string list
        Labels for the tvec data
    """
    while True:  # Loop until you read the end of file
        try:
            # Read single-line header
            record =  f.read_record('<f8', '<i4', '<i4', '<i4')
        except:
            break
        t = record[0]
        npl = record[1][0] - 1
        ntp = record[2][0]

        pvec = np.empty((6, npl))
        plid = np.empty(npl, dtype='int')
        tvec = np.empty((6, ntp))
        tpid = np.empty(ntp, dtype='int')
        if npl > 0:
            Mpl = np.empty(npl)
            Rpl = np.empty(npl)
            for i in range(npl):
                #Read single-line pl frame for
                record =  f.read_record('<i4', '<f8', '<f8', '(6,)<f8')
                plid[i] = record[0]
                Mpl[i] = record[1]
                Rpl[i] = record[2]
                pvec[:,i] = record[3]
        if ntp > 0:
            for i in range(ntp):
                record =  f.read_record('<i4', '(6,)<f8')
                tpid[i] = record[0]
                tvec[:,i] = record[1]

        tlab = []
        if param['OUT_FORM'] == 'XV':
            tlab.append('px')
            tlab.append('py')
            tlab.append('pz')
            tlab.append('vx')
            tlab.append('vy')
            tlab.append('vz')
        elif param['OUT_FORM'] == 'EL':
            tlab.append('a')
            tlab.append('e')
            tlab.append('inc')
            tlab.append('capom')
            tlab.append('omega')
            tlab.append('capm')
        plab = tlab.copy()
        plab.append('Mass')
        plab.append('Radius')
        pvec = np.vstack([pvec,Mpl,Rpl])

        yield t, npl, plid, pvec.T, plab, \
              ntp, tpid, tvec.T, tlab


def make_swiftest_labels(param):
    tlab = []
    if param['OUT_FORM'] == 'XV':
        tlab.append('px')
        tlab.append('py')
        tlab.append('pz')
        tlab.append('vx')
        tlab.append('vy')
        tlab.append('vz')
    elif param['OUT_FORM'] == 'EL':
        tlab.append('a')
        tlab.append('e')
        tlab.append('inc')
        tlab.append('capom')
        tlab.append('omega')
        tlab.append('capm')
    plab = tlab.copy()
    plab.append('Mass')
    plab.append('Radius')
    clab = ['Mass', 'Radius', 'J_2', 'J_4']
    if param['ROTATION'] == 'YES':
        clab.append('Ip_x')
        clab.append('Ip_y')
        clab.append('Ip_z')
        clab.append('rot_x')
        clab.append('rot_y')
        clab.append('rot_z')
        plab.append('Ip_x')
        plab.append('Ip_y')
        plab.append('Ip_z')
        plab.append('rot_x')
        plab.append('rot_y')
        plab.append('rot_z')
    if param['TIDES'] == 'YES':
        clab.append('k2')
        clab.append('Q')
        plab.append('k2')
        plab.append('Q')
    return clab, plab, tlab

def swiftest_stream(f, param):
    """
    Reads in a Swifter bin.dat file and returns a single frame of data as a datastream

    Parameters
    ----------
    f : file object
    param : dict

    Yields
    -------
    t    : float
        Time of this frame
    cbid : int array
        ID of central body (always returns 0)
    cvec : float array
        (npl,1) - vector of quantities for the massive body (Mass, Radius, J2, J4, etc)
    npl  : int
        Number of massive bodies
    plid : int array
        IDs of massive bodies
    pvec : float array
        (npl,N) - vector of N quantities or each particle (6 of XV/EL + Mass, Radius, etc)
    plab : string list
        Labels for the pvec data
    ntp  : int
        Number of test particles
    tpid : int array
        Ids of test particles
    tvec : float array
        (ntp,N) - vector of N quantities for each particle (6 of XV/EL, etc.)
    tlab : string list
        Labels for the tvec data
    """
    while True:  # Loop until you read the end of file
        try:
            # Read multi-line header
            t = f.read_reals(np.float64)  # Try first part of the header
        except:
            break
        npl = f.read_ints()
        ntp = f.read_ints()
        iout_form = f.read_reals('c')
        Mcb = f.read_reals(np.float64)
        Rcb = f.read_reals(np.float64)
        J2cb = f.read_reals(np.float64)
        J4cb = f.read_reals(np.float64)
        if param['ROTATION'] == 'YES':
            Ipcbx  = f.read_reals(np.float64)
            Ipcby  = f.read_reals(np.float64)
            Ipcbz  = f.read_reals(np.float64)
            rotcbx = f.read_reals(np.float64)
            rotcby = f.read_reals(np.float64)
            rotcbz = f.read_reals(np.float64)
        if param['TIDES'] == 'YES':
            k2cb = f.read_reals(np.float64)
            Qcb = f.read_reals(np.float64)
        if npl[0] > 0:
            plid = f.read_ints()
            p1 = f.read_reals(np.float64)
            p2 = f.read_reals(np.float64)
            p3 = f.read_reals(np.float64)
            p4 = f.read_reals(np.float64)
            p5 = f.read_reals(np.float64)
            p6 = f.read_reals(np.float64)
            Mpl = f.read_reals(np.float64)
            Rpl = f.read_reals(np.float64)
            if param['ROTATION'] == 'YES':
                Ipplx = f.read_reals(np.float64)
                Ipply = f.read_reals(np.float64)
                Ipplz = f.read_reals(np.float64)
                rotplx = f.read_reals(np.float64)
                rotply = f.read_reals(np.float64)
                rotplz = f.read_reals(np.float64)
            if param['TIDES'] == 'YES':
                k2pl = f.read_reals(np.float64)
                Qpl = f.read_reals(np.float64)
        if ntp[0] > 0:
            tpid = f.read_ints()
            t1 = f.read_reals(np.float64)
            t2 = f.read_reals(np.float64)
            t3 = f.read_reals(np.float64)
            t4 = f.read_reals(np.float64)
            t5 = f.read_reals(np.float64)
            t6 = f.read_reals(np.float64)
        cbid = np.array([0])

        clab, plab, tlab = make_swiftest_labels(param)

        if npl > 0:
            pvec = np.vstack([p1,p2,p3,p4,p5,p6,Mpl,Rpl])
        else:
            pvec = np.empty((8,0))
            plid = np.empty(0)
        if ntp > 0:
            tvec = np.vstack([t1,t2,t3,t4,t5,t6])
        else:
            tvec = np.empty((6,0))
            tpid = np.empty(0)
        cvec = np.array([Mcb,Rcb,J2cb,J4cb])
        if param['ROTATION'] == 'YES':
            cvec = np.vstack([cvec, Ipcbx, Ipcby, Ipcbz, rotcbx, rotcby, rotcbz])
            if npl > 0:
                pvec = np.vstack([pvec, Ipplx, Ipply, Ipplz, rotplx, rotply, rotplz])
        if param['TIDES'] == 'YES':
            cvec = np.vstack([cvec,k2cb,Qcb])
            if npl > 0:
                pvec = np.vstack([pvec,k2pl,Qpl])
        yield t, cbid, cvec.T, clab, \
              npl, plid, pvec.T, plab, \
              ntp, tpid, tvec.T, tlab

def swifter2xr(param):
    dims  = ['time','id', 'vec']
    pl = []
    tp = []
    with FortranFile(param['BIN_OUT'], 'r') as f:
        for t, npl, plid, pvec, plab, \
              ntp, tpid, tvec, tlab in swifter_stream(f, param):

            #Prepare frames by adding an extra axis for the time coordinate
            plframe = np.expand_dims(pvec, axis=0)
            tpframe = np.expand_dims(tvec, axis=0)

            #Create xarray DataArrays out of each body type
            plxr = xr.DataArray(plframe, dims = dims, coords = {'time' : t, 'id' : plid, 'vec' : plab})
            tpxr = xr.DataArray(tpframe, dims = dims, coords = {'time' : t, 'id' : tpid, 'vec' : tlab})

            pl.append(plxr)
            tp.append(tpxr)

        plda = xr.concat(pl, dim='time')
        tpda = xr.concat(tp, dim='time')

        plds = plda.to_dataset(dim='vec')
        tpds = tpda.to_dataset(dim='vec')
        ds = xr.combine_by_coords([plds, tpds])
    return ds

def swiftest2xr(param):
    """
    Converts a Swiftest binary data file into an xarray DataSet.

    Parameters
    ----------
    param : dict
        Swiftest paramuration parameters

    Returns
    -------
    xarray dataset
    """

    dims  = ['time','id', 'vec']
    cb = []
    pl = []
    tp = []
    with FortranFile(param['BIN_OUT'], 'r') as f:
        for t, cbid, cvec, clab, \
              npl, plid, pvec, plab, \
              ntp, tpid, tvec, tlab in swiftest_stream(f, param):

            #Prepare frames by adding an extra axis for the time coordinate
            cbframe = np.expand_dims(cvec, axis=0)
            plframe = np.expand_dims(pvec, axis=0)
            tpframe = np.expand_dims(tvec, axis=0)

            #Create xarray DataArrays out of each body type
            cbxr = xr.DataArray(cbframe, dims = dims, coords = {'time' : t, 'id' : cbid, 'vec' : clab})
            plxr = xr.DataArray(plframe, dims = dims, coords = {'time' : t, 'id' : plid, 'vec' : plab})
            tpxr = xr.DataArray(tpframe, dims = dims, coords = {'time' : t, 'id' : tpid, 'vec' : tlab})

            cb.append(cbxr)
            pl.append(plxr)
            tp.append(tpxr)

    cbda = xr.concat(cb,dim='time')
    plda = xr.concat(pl,dim='time')
    tpda = xr.concat(tp,dim='time')

    cbds = cbda.to_dataset(dim = 'vec')
    plds = plda.to_dataset(dim = 'vec')
    tpds = tpda.to_dataset(dim = 'vec')
    ds = xr.combine_by_coords([cbds, plds, tpds])
    return ds

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

    # Planet radii in meters
    planetradius = {
        'mercury': np.longdouble(2439.4e3),
        'venus': np.longdouble(6051.8e3),
        'earthmoon': np.longdouble(6371.0084e3),  # Earth only for radius
        'mars': np.longdouble(3389.50e3),
        'jupiter': np.longdouble(69911e3),
        'saturn': np.longdouble(58232.0e3),
        'uranus': np.longdouble(25362.e3),
        'neptune': np.longdouble(24622.e3),
        'plutocharon': np.longdouble(1188.3e3)
    }

    # Unit conversion factors
    DCONV = AU2M / param['DU2M']
    VCONV = (AU2M / JD2S) / (param['DU2M'] / param['TU2S'])
    THIRDLONG = np.longdouble(1.0) / np.longdouble(3.0)

    # Central body value vectors
    GMcb = np.array([GMSunSI * param['TU2S']**2 / param['DU2M']**3])
    Rcb = np.array([RSun / param['DU2M']])
    J2RP2 = np.array([J2Sun * (RSun / param['DU2M'])**2])
    J4RP4 = np.array([J4Sun * (RSun / param['DU2M'])**4])
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
        elif param['OUT_FORM'] == 'EL':
            p1.append(pldata[key].elements()['a'][0] * DCONV)
            p2.append(pldata[key].elements()['e'][0])
            p3.append(pldata[key].elements()['inc'][0] * np.pi / 180.0)
            p4.append(pldata[key].elements()['Omega'][0] * np.pi / 180.0)
            p5.append(pldata[key].elements()['w'][0] * np.pi / 180.0)
            p6.append(pldata[key].elements()['M'][0] * np.pi / 180.0)
        Rhill.append(pldata[key].elements()['a'][0] * (3 * MSun_over_Mpl[key]) ** (-THIRDLONG))
        Rpl.append(planetradius[key] * DCONV)
        GMpl.append(GMcb[0] / MSun_over_Mpl[key])
    # Generate planet value vectors
    plid = np.fromiter(planetid.values(), dtype=int)
    pvec = np.vstack([p1, p2, p3, p4, p5, p6, GMpl, Rpl])

    dims = ['time', 'id', 'vec']
    cb = []
    pl = []
    tp = []
    t = np.array([0.0])

    clab, plab, tlab = make_swiftest_labels(param)

    #Prepare frames by adding an extra axis for the time coordinate
    cbframe = np.expand_dims(cvec.T, axis=0)
    plframe = np.expand_dims(pvec.T, axis=0)

    #Create xarray DataArrays out of each body type
    cbxr = xr.DataArray(cbframe, dims = dims, coords = {'time' : t, 'id' : cbid, 'vec' : clab})
    plxr = xr.DataArray(plframe, dims = dims, coords = {'time' : t, 'id' : plid, 'vec' : plab})

    cb.append(cbxr)
    pl.append(plxr)

    cbda = xr.concat(cb,dim='time')
    plda = xr.concat(pl,dim='time')

    cbds = cbda.to_dataset(dim = 'vec')
    plds = plda.to_dataset(dim = 'vec')
    ds = xr.combine_by_coords([cbds, plds])
    return ds

def swiftest_xr2_infile(ds, param, framenum=-1):
    """
    Writes a set of Swiftest input files from a single frame of a Swiftest xarray dataset

    Parameters
    ----------
    ds : xarray dataset
        Dataset containing Swiftest n-body data in XV format
    framenum : int
        Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.
    param : dict
        Swiftest paramuration parameters. This method uses the names of the cb, pl, and tp files from the paramuration

    Returns
    -------
    A set of three input files for a Swiftest run
    """
    frame = ds.isel(time=framenum)
    cb = frame.where(frame.id == 0, drop=True)
    pl = frame.where(frame.id > 0, drop=True)
    pl = pl.where(np.invert(np.isnan(pl['Mass'])), drop=True).drop_vars(['J_2', 'J_4'])
    tp = frame.where(np.isnan(frame['Mass']), drop=True).drop_vars(['Mass', 'Radius', 'J_2', 'J_4'])

    GMSun = np.double(cb['Mass'])
    RSun = np.double(cb['Radius'])
    J2 = np.double(cb['J_2'])
    J4 = np.double(cb['J_4'])

    if param['IN_TYPE'] == 'ASCII':
        # Swiftest Central body file
        cbfile = open(param['CB_IN'], 'w')
        print(GMSun, file=cbfile)
        print(RSun, file=cbfile)
        print(J2, file=cbfile)
        print(J4, file=cbfile)
        cbfile.close()

        plfile = open(param['PL_IN'], 'w')
        print(pl.id.count().values, file=plfile)
        for i in pl.id:
            pli = pl.sel(id=i)
            print(i.values, pli['Mass'].values, file=plfile)
            print(pli['Radius'].values, file=plfile)
            print(pli['px'].values, pli['py'].values, pli['pz'].values, file=plfile)
            print(pli['vx'].values, pli['vy'].values, pli['vz'].values, file=plfile)
        plfile.close()

        # TP file
        tpfile = open(param['TP_IN'], 'w')
        print(tp.id.count().values, file=tpfile)
        for i in tp.id:
            tpi = tp.sel(id=i)
            print(i.values, file=tpfile)
            print(tpi['px'].values, tpi['py'].values, tpi['pz'].values, file=tpfile)
            print(tpi['vx'].values, tpi['vy'].values, tpi['vz'].values, file=tpfile)
        tpfile.close()
    elif param['IN_TYPE'] == 'REAL8':
        # Now make Swiftest files
        cbfile = FortranFile(swiftest_cb, 'w')
        MSun = np.double(1.0)
        cbfile.write_record(np.double(GMSun))
        cbfile.write_record(np.double(rmin))
        cbfile.write_record(np.double(J2))
        cbfile.write_record(np.double(J4))
        cbfile.close()

        plfile = FortranFile(swiftest_pl, 'w')
        plfile.write_record(npl)

        plfile.write_record(plid)
        plfile.write_record(p_pl[0])
        plfile.write_record(p_pl[1])
        plfile.write_record(p_pl[2])
        plfile.write_record(v_pl[0])
        plfile.write_record(v_pl[1])
        plfile.write_record(v_pl[2])
        plfile.write_record(mass)
        plfile.write_record(radius)
        plfile.close()
        tpfile = FortranFile(swiftest_tp, 'w')
        ntp = 1
        tpfile.write_record(ntp)
        tpfile.write_record(tpid)
        tpfile.write_record(p_tp[0])
        tpfile.write_record(p_tp[1])
        tpfile.write_record(p_tp[2])
        tpfile.write_record(v_tp[0])
        tpfile.write_record(v_tp[1])
        tpfile.write_record(v_tp[2])
    else:
        print(f"{param['IN_TYPE']} is an unknown file type")

def swifter2swiftest(inparam,outparam):
    print(f"Swifter parameter is {inparam}")
    print(f"Swiftest parameter file is {outparam}")
    swifter_param = read_swifter_param(inparam)
    swiftest_param = swifter_param.copy()
    print("Select the unit system to use:")
    print("1) MSun-AU-year")
    print("2) MSun-AU-day")
    print("3) SI: kg-m-s")
    print("4) CGS: g-cm-s")
    print("5) Set units manually")
    inval = input("> ")
    try:
        unit_type = int(inval)
    except ValueError:
        goodval = False
    else:
        goodval = (unit_type > 0 and unit_type < 6)
    if not goodval:
        print(f"{inval} is not a valid menu option")
        sys.exit(-1)
    if unit_type == 1:
        print("Unit system is MSun-AU-year")
        swiftest_param['MU2KG'] = MSun
        swiftest_param['DU2M'] = AU2M
        swiftest_param['TU2S'] = YR2S
    elif unit_type == 2: # MSun-AU-day
        print("Unit system is MSun-AU-day")
        swiftest_param['MU2KG'] = MSun
        swiftest_param['DU2M'] = AU2M
        swiftest_param['TU2S'] = JD2S
    elif unit_type == 3: # SI: kg-m-s
        print("Unit system is SI: kg-m-s")
        swiftest_param['MU2KG'] = 1.0
        swiftest_param['DU2M'] = 1.0
        swiftest_param['TU2S'] = 1.0
    elif unit_type == 4:  # CGS: g-cm-s
        print("Unit system is CGS: g-cm-s")
        swiftest_param['MU2KG'] = 1e-3
        swiftest_param['DU2M'] = 1.0e-2
        swiftest_param['TU2S'] = 1.0
    elif unit_type == 5:
        print("User-defined units.")
        print("Define each unit (mass, distance, and time) by its corresponding SI value.")
        swiftest_param['MU2KG'] = input("Mass value in kilograms: ")
        swiftest_param['DU2M'] = input("Distance value in meters: ")
        swiftest_param['TU2S'] = input("Time unit in seconds: ")

    print("Set central body radius:")
    print(f"1) Use Swifter parameter value of CHK_RMIN = {swifter_param['CHK_RMIN']}")
    print(f"2) Set value manually")
    inval = input("> ")
    try:
        cbrad_type = int(inval)
    except ValueError:
        goodval = False
    else:
        goodval = (cbrad_type > 0 and cbrad_type < 3)
    if not goodval:
        print(f"{inval} is not a valid menu option")
        sys.exit(-1)
    if cbrad_type == 1:
        cbrad = swifter_param['CHK_RMIN']
    elif cbrad_type == 2:
        cbrad = input("Enter radius of central body in simulation Distance Units: ")

    swiftest_param['PL_IN'] = 'pl.swiftest.in'
    swiftest_param['TP_IN'] = 'tp.swiftest.in'
    swiftest_param['CB_IN'] = 'cb.swiftest.in'

    plnew = open(swiftest_param['PL_IN'], 'w')

    print(f'Writing out new PL file: {swiftest_param["PL_IN"]}')
    with open(swifter_param['PL_IN'], 'r') as plold:
        line = plold.readline()
        line = line.split("!")[0]  # Ignore comments
        i_list = [i for i in line.split(" ") if i.strip()]
        npl = int(i_list[0])
        print(npl - 1, file=plnew)
        line = plold.readline()
        i_list = [i for i in line.split(" ") if i.strip()]
        GMcb = float(i_list[1]) # Store central body GM for later
        line = plold.readline()  # Ignore the two zero vector lines
        line = plold.readline()
        for n in range(1, npl):  # Loop over planets
            line = plold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            name = int(i_list[0])
            GMpl = float(i_list[1])
            print(name, GMpl,file=plnew)
            if swifter_param['CHK_CLOSE'] == 'YES':
                line = plold.readline()
                i_list = [i for i in line.split(" ") if i.strip()]
                plrad = float(i_list[0])
                print(plrad,file=plnew)
            line = plold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            xh = float(i_list[0])
            yh = float(i_list[1])
            zh = float(i_list[2])
            print(xh, yh, zh, file=plnew)
            line = plold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            vx = float(i_list[0])
            vy = float(i_list[1])
            vz = float(i_list[2])
            print(vx, vy, vz, file=plnew)

    plold.close()
    plnew.close()
    tpnew = open(swiftest_param['TP_IN'], 'w')

    print(f'Writing out new TP file: {swiftest_param["TP_IN"]}')
    with open(swifter_param['TP_IN'], 'r') as tpold:
        line = tpold.readline()
        line = line.split("!")[0]  # Ignore comments
        i_list = [i for i in line.split(" ") if i.strip()]
        ntp = int(i_list[0])
        print(ntp, file=tpnew)
        for n in range(0, ntp):  # Loop over test particles
            line = tpold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            name = int(i_list[0])
            print(name, file=tpnew)
            line = tpold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            xh = float(i_list[0])
            yh = float(i_list[1])
            zh = float(i_list[2])
            print(xh, yh, zh, file=tpnew)
            line = tpold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            vx = float(i_list[0])
            vy = float(i_list[1])
            vz = float(i_list[2])
            print(vx, vy, vz, file=tpnew)

    tpold.close()
    tpnew.close()

    print(f'Writing out new CB file: {swiftest_param["CB_IN"]}')
    # Write out new central body file
    cbnew = open(swiftest_param['CB_IN'], 'w')

    print(GMcb, file=cbnew)
    print(cbrad, file=cbnew)
    print(swifter_param['J2'], file=cbnew)
    print(swifter_param['J4'], file=cbnew)

    cbnew.close()
    return

if __name__ == '__main__':

    workingdir = '/Users/daminton/git/swiftest/examples/rmvs_swifter_comparison/9pl_18tp_encounters/'
    #inparfile = workingdir + 'param.swifter.in'
    #param = read_swifter_param(inparfile)
    #param['BIN_OUT'] = workingdir + param['BIN_OUT']

    param_file_name = workingdir + 'param.swiftest.in'
    param = read_swiftest_param(param_file_name)
    param['BIN_OUT'] = workingdir + param['BIN_OUT']
    ds = solar_system_pl(param, '2020-06-17')
    ds

    #swiftestdat = swiftest2xr(param)
    param['CB_IN'] = workingdir + 'cb_test.in'
    param['PL_IN'] = workingdir + 'pl_test.in'
    param['TP_IN'] = workingdir + 'tp_test.in'
    swiftest_xr2_infile(ds, param)
    #swifterdat = swifter2xr(param)
    #print(swiftestdat['a'])
    #print(swiftestdf.head())

    #swifterdf.plot(y='px')
