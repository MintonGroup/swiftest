import swiftest
import numpy as np
from scipy.io import FortranFile
import xarray as xr
import sys
import tempfile

newfeaturelist = ("FRAGMENTATION", "ROTATION", "TIDES", "ENERGY", "GR", "YARKOVSKY", "YORP", "IN_FORM")

def real2float(realstr):
    """
    Converts a Fortran-generated ASCII string of a real value into a numpy float type. Handles cases where double precision
    numbers in exponential notation use 'd' or 'D' instead of 'e' or 'E'
    
    Parameters
    ----------
    realstr : string
        Fortran-generated ASCII string of a real value.

    Returns
    -------
     : float
        The converted floating point value of the input string
    """
    return float(realstr.replace('d', 'E').replace('D', 'E'))


def read_swiftest_param(param_file_name, param):
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
    param['! VERSION'] = f"Swiftest parameter input from file {param_file_name}"
    
    # Read param.in file
    print(f'Reading Swiftest file {param_file_name}')
    try:
        with open(param_file_name, 'r') as f:
            for line in f.readlines():
                fields = line.split()
                if len(fields) > 0:
                    for key in param:
                        if (key == fields[0].upper()): param[key] = fields[1]
                    # Special case of CHK_QMIN_RANGE requires a second input
                    if fields[0].upper() == 'CHK_QMIN_RANGE':
                        alo = real2float(fields[1])
                        ahi = real2float(fields[2])
                        param['CHK_QMIN_RANGE'] = f"{alo} {ahi}"
        
        param['ISTEP_OUT'] = int(param['ISTEP_OUT'])
        param['ISTEP_DUMP'] = int(param['ISTEP_DUMP'])
        param['T0'] = real2float(param['T0'])
        param['TSTOP'] = real2float(param['TSTOP'])
        param['DT'] = real2float(param['DT'])
        param['CHK_RMIN'] = real2float(param['CHK_RMIN'])
        param['CHK_RMAX'] = real2float(param['CHK_RMAX'])
        param['CHK_EJECT'] = real2float(param['CHK_EJECT'])
        param['CHK_QMIN'] = real2float(param['CHK_QMIN'])
        param['DU2M'] = real2float(param['DU2M'])
        param['MU2KG'] = real2float(param['MU2KG'])
        param['TU2S'] = real2float(param['TU2S'])
        param['EXTRA_FORCE'] = param['EXTRA_FORCE'].upper()
        param['BIG_DISCARD'] = param['BIG_DISCARD'].upper()
        param['CHK_CLOSE'] = param['CHK_CLOSE'].upper()
        param['RHILL_PRESENT'] = param['RHILL_PRESENT'].upper()
        param['FRAGMENTATION'] = param['FRAGMENTATION'].upper()
        if param['FRAGMENTATION'] == 'YES' and param['PARTICLE_OUT'] == '':
            param['PARTICLE_OUT'] = 'particle.dat'
        param['ROTATION'] = param['ROTATION'].upper()
        param['TIDES'] = param['TIDES'].upper()
        param['ENERGY'] = param['ENERGY'].upper()
        param['GR'] = param['GR'].upper()
        if 'GMTINY' in param:
            param['GMTINY'] = real2float(param['GMTINY'])
    except IOError:
        print(f"{param_file_name} not found.")
    return param


def read_swifter_param(param_file_name):
    """
    Reads in a Swifter param.in file and saves it as a dictionary

    Parameters
    ----------
    param_file_name : string
        File name of the input parameter file

    Returns
    -------
    param
        A dictionary containing the entries in the user parameter file
    """
    param = {
        '! VERSION': f"Swifter parameter input from file {param_file_name}",
        'T0': "0.0",
        'TSTOP': "0.0",
        'DT': "0.0",
        'PL_IN': "",
        'TP_IN': "",
        'IN_TYPE': "ASCII",
        'ISTEP_OUT': "-1",
        'ISTEP_DUMP': "-1",
        'BIN_OUT': "bin.dat",
        'OUT_TYPE': "REAL8",
        'OUT_FORM': "XV",
        'OUT_STAT': "NEW",
        'J2': "0.0",
        'J4': "0.0",
        'CHK_CLOSE': 'NO',
        'CHK_RMIN': "-1.0",
        'CHK_RMAX': "-1.0",
        'CHK_EJECT': "-1.0",
        'CHK_QMIN': "-1.0",
        'CHK_QMIN_COORD': "HELIO",
        'CHK_QMIN_RANGE': "",
        'ENC_OUT': "",
        'EXTRA_FORCE': 'NO',
        'BIG_DISCARD': 'NO',
        'RHILL_PRESENT': 'NO',
        'C': "-1.0",
    }
    
    # Read param.in file
    print(f'Reading Swifter file {param_file_name}')
    try:
        with open(param_file_name, 'r') as f:
            for line in f.readlines():
                fields = line.split()
                if len(fields) > 0:
                    for key in param:
                        if (key == fields[0].upper()): param[key] = fields[1]
                    # Special case of CHK_QMIN_RANGE requires a second input
                    if fields[0].upper() == 'CHK_QMIN_RANGE':
                        alo = real2float(fields[1])
                        ahi = real2float(fields[2])
                        param['CHK_QMIN_RANGE'] = f"{alo} {ahi}"
    
        param['ISTEP_OUT'] = int(param['ISTEP_OUT'])
        param['ISTEP_DUMP'] = int(param['ISTEP_DUMP'])
        param['T0'] = real2float(param['T0'])
        param['TSTOP'] = real2float(param['TSTOP'])
        param['DT'] = real2float(param['DT'])
        param['J2'] = real2float(param['J2'])
        param['J4'] = real2float(param['J4'])
        param['CHK_RMIN'] = real2float(param['CHK_RMIN'])
        param['CHK_RMAX'] = real2float(param['CHK_RMAX'])
        param['CHK_EJECT'] = real2float(param['CHK_EJECT'])
        param['CHK_QMIN'] = real2float(param['CHK_QMIN'])
        param['EXTRA_FORCE'] = param['EXTRA_FORCE'].upper()
        param['BIG_DISCARD'] = param['BIG_DISCARD'].upper()
        param['CHK_CLOSE'] = param['CHK_CLOSE'].upper()
        param['RHILL_PRESENT'] = param['RHILL_PRESENT'].upper()
        if param['C'] != '-1.0':
            param['C'] = real2float(param['C'])
        else:
            param.pop('C', None)
    except IOError:
        print(f"{param_file_name} not found.")

    return param


def read_swift_param(param_file_name, startfile="swift.in"):
    """
    Reads in a Swift param.in file and saves it as a dictionary

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
        '! VERSION': f"Swift parameter input from file {param_file_name}",
        'T0': 0.0,
        'TSTOP': 0.0,
        'DT': 0.0,
        'DTOUT': 0.0,
        'DTDUMP': 0.0,
        'L1': "F",
        'L1': "F",
        'L2': "F",
        'L3': "F",
        'L4': "F",
        'L5': "F",
        'L6': "F",
        'RMIN': -1,
        'RMAX': -1,
        'RMAXU': -1,
        'QMIN': -1,
        'LCLOSE': "F",
        'BINARY_OUTPUTFILE': "bin.dat",
        'STATUS_FLAG_FOR_OPEN_STATEMENTS': "NEW",
    }
    
    try:
        with open(startfile, 'r') as f:
            line = f.readline()
            plname = f.readline().split()[0]
            tpname = f.readline().split()[0]
    except:
        plname = "pl.in"
        tpname = "tp.in"
    param['PL_IN'] = plname
    param['TP_IN'] = tpname

    # Read param.in file
    print(f'Reading Swift file {param_file_name}')
    try:
        with open(param_file_name, 'r') as f:
            line = f.readline().split()
            for i, l in enumerate(line):
                line[i] = l
            param['T0'] = real2float(line[0])
            param['TSTOP'] = real2float(line[1])
            param['DT'] = real2float(line[2])
            line = f.readline().split()
            for i, l in enumerate(line):
                line[i] = l
            param['DTOUT'] = real2float(line[0])
            param['DTDUMP'] = real2float(line[1])
            line = f.readline().split()
            param['L1'] = line[0].upper()
            param['L2'] = line[1].upper()
            param['L3'] = line[2].upper()
            param['L4'] = line[3].upper()
            param['L5'] = line[4].upper()
            param['L6'] = line[5].upper()
            if param['L2'] == "T":
                line = f.readline().split()
                for i, l in enumerate(line):
                    line[i] = l
                param['RMIN'] = real2float(line[0])
                param['RMAX'] = real2float(line[1])
                param['RMAXU'] = real2float(line[2])
                param['QMIN'] = real2float(line[3])
                param['LCLOSE'] = line[4].upper()
            param['BINARY_OUTPUTFILE'] = f.readline().strip()
            param['STATUS_FLAG_FOR_OPEN_STATEMENTS'] = f.readline().strip().upper()
    except IOError:
       print(f"{param_file_name} not found.")
    
    return param


def write_swift_param(param, param_file_name):
    outfile = open(param_file_name, 'w')
    print(param['T0'], param['TSTOP'], param['DT'], file=outfile)
    print(param['DTOUT'], param['DTDUMP'], file=outfile)
    print(param['L1'], param['L2'], param['L3'], param['L4'], param['L5'], param['L6'], file=outfile)
    print(param['RMIN'], param['RMAX'], param['RMAXU'], param['QMIN'], param['LCLOSE'], file=outfile)
    print(param['BINARY_OUTPUTFILE'], file=outfile)
    print(param['STATUS_FLAG_FOR_OPEN_STATEMENTS'], file=outfile)
    outfile.close()
    return


def write_labeled_param(param, param_file_name):
    outfile = open(param_file_name, 'w')
    keylist = ['! VERSION',
               'T0',
               'TSTOP',
               'DT',
               'ISTEP_OUT',
               'ISTEP_DUMP',
               'OUT_FORM',
               'OUT_TYPE',
               'OUT_STAT',
               'IN_TYPE',
               'PL_IN',
               'TP_IN',
               'CB_IN',
               'BIN_OUT',
               'CHK_QMIN',
               'CHK_RMIN',
               'CHK_RMAX',
               'CHK_EJECT',
               'CHK_QMIN_COORD',
               'CHK_QMIN_RANGE',
               'MU2KG',
               'TU2S',
               'DU2M' ]
    ptmp = param.copy()
    # Print the list of key/value pairs in the preferred order
    for key in keylist:
        val = ptmp.pop(key, None)
        if val is not None: print(f"{key:<16} {val}", file=outfile)
    # Print the remaining key/value pairs in whatever order
    for key, val in ptmp.items():
        if val != "": print(f"{key:<16} {val}", file=outfile)
    outfile.close()
    return


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
        (npl,N) - vector of N quantities or each particle (6 of XV/EL + Gmass, radius, etc)
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
            record = f.read_record('<f8', '<i4', '<i4', '<i4')
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
            GMpl = np.empty(npl)
            Rpl = np.empty(npl)
            for i in range(npl):
                # Read single-line pl frame for
                record = f.read_record('<i4', '<f8', '<f8', '(6,)<f8')
                plid[i] = record[0]
                GMpl[i] = record[1]
                Rpl[i] = record[2]
                pvec[:, i] = record[3]
        if ntp > 0:
            for i in range(ntp):
                record = f.read_record('<i4', '(6,)<f8')
                tpid[i] = record[0]
                tvec[:, i] = record[1]
        
        tlab = []
        if param['OUT_FORM'] == 'XV' or param['OUT_FORM'] == 'XVEL':
            tlab.append('xhx')
            tlab.append('xhy')
            tlab.append('xhz')
            tlab.append('vhx')
            tlab.append('vhy')
            tlab.append('vhz')
        if param['OUT_FORM'] == 'EL' or param['OUT_FORM'] == 'XVEL':
            tlab.append('a')
            tlab.append('e')
            tlab.append('inc')
            tlab.append('capom')
            tlab.append('omega')
            tlab.append('capm')
        plab = tlab.copy()
        plab.append('Gmass')
        plab.append('radius')
        pvec = np.vstack([pvec, GMpl, Rpl])
        
        yield t, npl, plid, pvec.T, plab, \
              ntp, tpid, tvec.T, tlab


def make_swiftest_labels(param):
    tlab = []
    if param['OUT_FORM'] == 'XV' or param['OUT_FORM'] == 'XVEL':
        tlab.append('xhx')
        tlab.append('xhy')
        tlab.append('xhz')
        tlab.append('vhx')
        tlab.append('vhy')
        tlab.append('vhz')
   
    if param['OUT_FORM'] == 'EL' or param['OUT_FORM'] == 'XVEL':
        tlab.append('a')
        tlab.append('e')
        tlab.append('inc')
        tlab.append('capom')
        tlab.append('omega')
        tlab.append('capm')
    plab = tlab.copy()
    plab.append('Gmass')
    plab.append('radius')
    if param['RHILL_PRESENT'] == 'YES':
        plab.append('rhill')
    clab = ['Gmass', 'radius', 'J_2', 'J_4']
    if param['ROTATION'] == 'YES':
        clab.append('Ip1')
        clab.append('Ip2')
        clab.append('Ip3')
        clab.append('rotx')
        clab.append('roty')
        clab.append('rotz')
        plab.append('Ip1')
        plab.append('Ip2')
        plab.append('Ip3')
        plab.append('rotx')
        plab.append('roty')
        plab.append('rotz')
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
        (npl,1) - vector of quantities for the massive body (Gmass, radius, J2, J4, etc)
    npl  : int
        Number of massive bodies
    plid : int array
        IDs of massive bodies
    pvec : float array
        (npl,N) - vector of N quantities or each particle (6 of XV/EL + Gmass, radius, etc)
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
    NAMELEN = 32
    while True:  # Loop until you read the end of file
        try:
            # Read multi-line header
            t = f.read_reals(np.float64)  # Try first part of the header
        except:
            break
        npl = f.read_ints()
        ntp = f.read_ints()
        iout_form = f.read_reals('c')
        cbid = f.read_ints()
        dtstr = f'a{NAMELEN}'
        cbnames = f.read_record(dtstr)
        cbnames = [np.char.strip(str(cbnames[0], encoding='utf-8'))]
        Mcb = f.read_reals(np.float64)
        Rcb = f.read_reals(np.float64)
        J2cb = f.read_reals(np.float64)
        J4cb = f.read_reals(np.float64)
        if param['ROTATION'] == 'YES':
            Ipcbx = f.read_reals(np.float64)
            Ipcby = f.read_reals(np.float64)
            Ipcbz = f.read_reals(np.float64)
            rotcbx = f.read_reals(np.float64)
            rotcby = f.read_reals(np.float64)
            rotcbz = f.read_reals(np.float64)
        if param['TIDES'] == 'YES':
            k2cb = f.read_reals(np.float64)
            Qcb = f.read_reals(np.float64)
        if npl[0] > 0:
            plid = f.read_ints()
            dtstr = f'({npl[0]},)a{NAMELEN}'
            names = f.read_record(np.dtype(dtstr))
            plnames = []
            for i in range(npl[0]):
               plnames.append(np.char.strip(str(names[i], encoding='utf-8')))
            p1 = f.read_reals(np.float64)
            p2 = f.read_reals(np.float64)
            p3 = f.read_reals(np.float64)
            p4 = f.read_reals(np.float64)
            p5 = f.read_reals(np.float64)
            p6 = f.read_reals(np.float64)
            if param['OUT_FORM'] == 'XVEL':
               p7 = f.read_reals(np.float64)
               p8 = f.read_reals(np.float64)
               p9 = f.read_reals(np.float64)
               p10 = f.read_reals(np.float64)
               p11 = f.read_reals(np.float64)
               p12 = f.read_reals(np.float64) 
            GMpl = f.read_reals(np.float64)
            if param['RHILL_PRESENT'] == 'YES':
                rhill = f.read_reals(np.float64)
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
            dtstr = f'({ntp[0]},)a{NAMELEN}'
            names = f.read_record(np.dtype(dtstr))
            tpnames = []
            for i in range(npl[0]):
               tpnames.append(np.char.strip(str(names[i], encoding='utf-8')))
            t1 = f.read_reals(np.float64)
            t2 = f.read_reals(np.float64)
            t3 = f.read_reals(np.float64)
            t4 = f.read_reals(np.float64)
            t5 = f.read_reals(np.float64)
            t6 = f.read_reals(np.float64)
            if param['OUT_FORM'] == 'XVEL':
               t7 = f.read_reals(np.float64)
               t8 = f.read_reals(np.float64)
               t9 = f.read_reals(np.float64)
               t10 = f.read_reals(np.float64)
               t11 = f.read_reals(np.float64)
               t12 = f.read_reals(np.float64) 
        
        clab, plab, tlab = make_swiftest_labels(param)

        if npl > 0:
            pvec = np.vstack([p1, p2, p3, p4, p5, p6])
            if param['OUT_FORM'] == 'XVEL':
               pvec = np.vstack([pvec, p7, p8, p9, p10, p11, p12])
            pvec = np.vstack([pvec, GMpl, Rpl])
        else:
            if param['OUT_FORM'] == 'XVEL':
               pvec = np.empty((14, 0))
            else:
               pvec = np.empty((8, 0))
            plid = np.empty(0)
            plnames = np.empty(0)
        if ntp > 0:
            tvec = np.vstack([t1, t2, t3, t4, t5, t6])
            if param['OUT_FORM'] == 'XVEL':
               tvec = np.vstack([tvec, t7, t8, t9, t10, t11, t12])
        else:
            if param['OUT_FORM'] == 'XVEL':
               tvec = np.empty((12, 0))
            else:
               tvec = np.empty((6, 0))
            tpid = np.empty(0)
            tpnames = np.empty(0)
        cvec = np.array([Mcb, Rcb, J2cb, J4cb])
        if param['RHILL_PRESENT'] == 'YES':
           if npl > 0:
              pvec = np.vstack([pvec, rhill])
        if param['ROTATION'] == 'YES':
            cvec = np.vstack([cvec, Ipcbx, Ipcby, Ipcbz, rotcbx, rotcby, rotcbz])
            if npl > 0:
                pvec = np.vstack([pvec, Ipplx, Ipply, Ipplz, rotplx, rotply, rotplz])
        if param['TIDES'] == 'YES':
            cvec = np.vstack([cvec, k2cb, Qcb])
            if npl > 0:
                pvec = np.vstack([pvec, k2pl, Qpl])
        yield t, cbid, cbnames, cvec.T, clab, \
              npl, plid, plnames, pvec.T, plab, \
              ntp, tpid, tpnames, tvec.T, tlab


def swifter2xr(param):
    """
    Converts a Swifter binary data file into an xarray DataSet.

    Parameters
    ----------
    param : dict
        Swifter parameters

    Returns
    -------
    xarray dataset
    """
    dims = ['time', 'id', 'vec']
    pl = []
    tp = []
    with FortranFile(param['BIN_OUT'], 'r') as f:
        for t, npl, plid, pvec, plab, \
            ntp, tpid, tvec, tlab in swifter_stream(f, param):
            # Prepare frames by adding an extra axis for the time coordinate
            plframe = np.expand_dims(pvec, axis=0)
            tpframe = np.expand_dims(tvec, axis=0)

            # Create xarray DataArrays out of each body type
            plxr = xr.DataArray(plframe, dims=dims, coords={'time': t, 'id': plid, 'vec': plab})
            tpxr = xr.DataArray(tpframe, dims=dims, coords={'time': t, 'id': tpid, 'vec': tlab})
            
            pl.append(plxr)
            tp.append(tpxr)
            sys.stdout.write('\r' + f"Reading in time {t[0]:.3e}")
            sys.stdout.flush()
        
        plda = xr.concat(pl, dim='time')
        tpda = xr.concat(tp, dim='time')
        
        plds = plda.to_dataset(dim='vec')
        tpds = tpda.to_dataset(dim='vec')
        print('\nCreating Dataset')
        ds = xr.combine_by_coords([plds, tpds])
        print(f"Successfully converted {ds.sizes['time']} output frames.")
    return ds


def swiftest2xr(param):
    """
    Converts a Swiftest binary data file into an xarray DataSet.

    Parameters
    ----------
    param : dict
        Swiftest parameters

    Returns
    -------
    xarray dataset
    """
    if ((param['OUT_TYPE'] == 'REAL8') or (param['OUT_TYPE'] == 'REAL4')): 
        dims = ['time', 'id', 'vec']
        cb = []
        pl = []
        tp = []
        cbn = None
        try:
            with FortranFile(param['BIN_OUT'], 'r') as f:
                for t, cbid, cbnames, cvec, clab, \
                    npl, plid, plnames, pvec, plab, \
                    ntp, tpid, tpnames, tvec, tlab in swiftest_stream(f, param):
                    # Prepare frames by adding an extra axis for the time coordinate
                    cbframe = np.expand_dims(cvec, axis=0)
                    plframe = np.expand_dims(pvec, axis=0)
                    tpframe = np.expand_dims(tvec, axis=0)


                    # Create xarray DataArrays out of each body type
                    cbxr = xr.DataArray(cbframe, dims=dims, coords={'time': t, 'id': cbid, 'vec': clab})
                    cbxr = cbxr.assign_coords(name=("id", cbnames))
                    plxr = xr.DataArray(plframe, dims=dims, coords={'time': t, 'id': plid, 'vec': plab})
                    plxr = plxr.assign_coords(name=("id", plnames))
                    tpxr = xr.DataArray(tpframe, dims=dims, coords={'time': t, 'id': tpid, 'vec': tlab})
                    tpxr = tpxr.assign_coords(name=("id", tpnames))

                    cb.append(cbxr)
                    pl.append(plxr)
                    tp.append(tpxr)

                    sys.stdout.write('\r' + f"Reading in time {t[0]:.3e}")
                    sys.stdout.flush()
        except IOError:
            print(f"Error encountered reading in {param['BIN_OUT']}")

        cbda = xr.concat(cb, dim='time')
        plda = xr.concat(pl, dim='time')
        tpda = xr.concat(tp, dim='time')
    
        cbds = cbda.to_dataset(dim='vec')
        plds = plda.to_dataset(dim='vec')
        tpds = tpda.to_dataset(dim='vec')
        print('\nCreating Dataset')
        ds = xr.combine_by_coords([cbds, plds, tpds])

    elif ((param['OUT_TYPE'] == 'NETCDF_DOUBLE') or (param['OUT_TYPE'] == 'NETCDF_FLOAT')):
        print('\nCreating Dataset')
        ds = xr.open_dataset(param['BIN_OUT'])
        ds = clean_string_values(param, ds)
    else:
        print(f"Error encountered. OUT_TYPE {param['OUT_TYPE']} not recognized.")
        return None
    print(f"Successfully converted {ds.sizes['time']} output frames.")

    if param['PARTICLE_OUT'] != "":
        infoxr = swiftest_particle_2xr(param)
        print('\nAdding particle info to Dataset')
        ds = xr.merge([ds, infoxr])

    return ds


def clean_string_values(param, ds):
    """
    Cleans up the string values in the DataSet that have artifacts as a result of coming from NetCDF Fortran

    Parameters
    ----------
    param : dict
    ds    : xarray dataset
 
    Returns
    -------
    ds : xarray dataset with the strings cleaned up
    """  
    ds['name'] =  ds['name'].str.decode(encoding='utf-8')
    return ds


def swiftest_particle_stream(f):
   """
   Reads in a Swiftest particle.dat file and returns a single frame of particle data as a datastream

   Parameters
   ----------
   f : file object
   param : dict

   Yields
   -------
   plid : int
      ID of massive bodie
   origin_type : string
      The origin type for the body (Initial conditions, disruption, supercatastrophic, hit and run, etc)
   origin_xh : float array
      The origin heliocentric position vector
   origin_vh : float array
      The origin heliocentric velocity vector
   """
   while True:  # Loop until you read the end of file
      try:
         # Read multi-line header
         plid = f.read_ints()  # Try first part of the header
      except:
         break
      origin_rec = f.read_record(np.dtype('a32'), np.dtype(('<f8', (7))))
      origin_type = np.char.strip(str(origin_rec[0], encoding='utf-8'))
      origin_vec = origin_rec[1]
      yield plid, origin_type, origin_vec


def swiftest_particle_2xr(param):
   """Reads in the Swiftest SyMBA-generated PARTICLE_OUT  and converts it to an xarray Dataset"""
   veclab = ['time_origin', 'xhx_origin', 'py_origin', 'pz_origin', 'vhx_origin', 'vhy_origin', 'vhz_origin']
   id_list = []
   origin_type_list = []
   origin_vec_list = []

   try:
      with FortranFile(param['PARTICLE_OUT'], 'r') as f:
          for id, origin_type, origin_vec in swiftest_particle_stream(f):
             id_list.append(id)
             origin_type_list.append(origin_type)
             origin_vec_list.append(origin_vec)
   except IOError:
       print(f"Error reading in {param['PARTICLE_OUT']} ")

   id_list =  np.asarray(id_list)[:,0]
   origin_type_list = np.asarray(origin_type_list)
   origin_vec_list = np.vstack(origin_vec_list)

   typeda = xr.DataArray(origin_type_list, dims=['id'], coords={'id' : id_list})
   vecda = xr.DataArray(origin_vec_list, dims=['id', 'vec'], coords={'id' : id_list, 'vec' : veclab})

   infoxr = vecda.to_dataset(dim='vec')
   infoxr['origin_type'] = typeda

   return infoxr


def swiftest_xr2infile(ds, param, framenum=-1):
    """
    Writes a set of Swiftest input files from a single frame of a Swiftest xarray dataset

    Parameters
    ----------
    ds : xarray dataset
        Dataset containing Swiftest n-body data in XV format
    framenum : int
        Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.
    param : dict
        Swiftest input parameters. This method uses the names of the cb, pl, and tp files from the input

    Returns
    -------
    A set of three input files for a Swiftest run
    """
    frame = ds.isel(time=framenum)
    cb = frame.where(frame.id == 1, drop=True)
    pl = frame.where(frame.id > 1, drop=True)
    pl = pl.where(np.invert(np.isnan(pl['Gmass'])), drop=True).drop_vars(['J_2', 'J_4'])
    tp = frame.where(np.isnan(frame['Gmass']), drop=True).drop_vars(['Gmass', 'radius', 'J_2', 'J_4'])
    
    GMSun = np.double(cb['Gmass'])
    RSun = np.double(cb['radius'])
    J2 = np.double(cb['J_2'])
    J4 = np.double(cb['J_4'])
    cbname = cb['name'].values[0]
    if param['ROTATION'] == 'YES':
        Ip1cb = np.double(cb['Ip1'])
        Ip2cb = np.double(cb['Ip2'])
        Ip3cb = np.double(cb['Ip3'])
        rotxcb = np.double(cb['rotx'])
        rotycb = np.double(cb['roty'])
        rotzcb = np.double(cb['rotz'])
    cbid = int(0)
    
    if param['IN_TYPE'] == 'ASCII':
        # Swiftest Central body file
        cbfile = open(param['CB_IN'], 'w')
        print(cbname, file=cbfile)
        print(GMSun, file=cbfile)
        print(RSun, file=cbfile)
        print(J2, file=cbfile)
        print(J4, file=cbfile)
        if param['ROTATION'] == 'YES':
            print(Ip1cb, Ip2cb, Ip3cb, file=cbfile)
            print(rotxcb, rotycb, rotzcb, file=cbfile)
        cbfile.close()
        
        plfile = open(param['PL_IN'], 'w')
        print(pl.id.count().values, file=plfile)
        for i in pl.id:
            pli = pl.sel(id=i)
            if param['RHILL_PRESENT'] == 'YES':
               print(pli['name'].values, pli['Gmass'].values, pli['rhill'].values, file=plfile)
            else:
               print(pli['name'].values, pli['Gmass'].values, file=plfile)
            print(pli['radius'].values, file=plfile)
            if param['IN_FORM'] == 'XV':
                print(pli['xhx'].values, pli['xhy'].values, pli['xhz'].values, file=plfile)
                print(pli['vhx'].values, pli['vhy'].values, pli['vhz'].values, file=plfile)
            elif param['IN_FORM'] == 'EL':
                print(pli['a'].values, pli['e'].values, pli['inc'].values, file=plfile)
                print(pli['capom'].values, pli['omega'].values, pli['capm'].values, file=plfile)
            else:
                print(f"{param['IN_FORM']} is not a valid input format type.")
            if param['ROTATION'] == 'YES':
                print(pli['Ip1'].values, pli['Ip2'].values, pli['Ip3'].values, file=plfile)
                print(pli['rotx'].values, pli['roty'].values, pli['rotz'].values, file=plfile)
        plfile.close()
        
        # TP file
        tpfile = open(param['TP_IN'], 'w')
        print(tp.id.count().values, file=tpfile)
        for i in tp.id:
            tpi = tp.sel(id=i)
            print(tpi['name'].values, file=tpfile)
            if param['IN_FORM'] == 'XV':
                print(tpi['xhx'].values, tpi['xhy'].values, tpi['xhz'].values, file=tpfile)
                print(tpi['vhx'].values, tpi['vhy'].values, tpi['vhz'].values, file=tpfile)
            elif param['IN_FORM'] == 'EL':
                print(tpi['a'].values, tpi['e'].values, tpi['inc'].values, file=tpfile)
                print(tpi['capom'].values, tpi['omega'].values, tpi['capm'].values, file=tpfile)
            else:
                print(f"{param['IN_FORM']} is not a valid input format type.")
        tpfile.close()
    elif param['IN_TYPE'] == 'REAL8':
        # Now make Swiftest files
        cbfile = FortranFile(param['CB_IN'], 'w')
        cbfile.write_record(cbid)
        cbfile.write_record(np.double(GMSun))
        cbfile.write_record(np.double(RSun))
        cbfile.write_record(np.double(J2))
        cbfile.write_record(np.double(J4))
        if param['ROTATION'] == 'YES':
            cbfile.write_record(np.double(Ip1cb))
            cbfile.write_record(np.double(Ip2cb))
            cbfile.write_record(np.double(Ip3cb))
            cbfile.write_record(np.double(rotxcb))
            cbfile.write_record(np.double(rotycb))
            cbfile.write_record(np.double(rotzcb))

        cbfile.close()
        
        plfile = FortranFile(param['PL_IN'], 'w')
        npl = pl.id.count().values
        plid = pl.id.values
        if param['IN_FORM'] == 'XV':
            v1 = pl['xhx'].values
            v2 = pl['xhy'].values
            v3 = pl['xhz'].values
            v4 = pl['vhx'].values
            v5 = pl['vhy'].values
            v6 = pl['vhz'].values
        elif param['IN_FORM'] == 'EL':
            v1 = pl['a'].values
            v2 = pl['e'].values
            v3 = pl['inc'].values
            v4 = pl['capom'].values
            v5 = pl['omega'].values
            v6 = pl['capm'].values
        else:
            print(f"{param['IN_FORM']} is not a valid input format type.")
        Gmass = pl['Gmass'].values
        radius = pl['radius'].values  
        
        plfile.write_record(npl)
        plfile.write_record(plid)
        plfile.write_record(v1)
        plfile.write_record(v2)
        plfile.write_record(v3)
        plfile.write_record(v4)
        plfile.write_record(v5)
        plfile.write_record(v6)
        plfile.write_record(Gmass)
        if param['RHILL_PRESENT'] == 'YES':
            plfile.write_record(pl['rhill'].values)
        plfile.write_record(radius)
        if param['ROTATION'] == 'YES':
            plfile.write_record(pl['Ip1'].values)
            plfile.write_record(pl['Ip2'].values)
            plfile.write_record(pl['Ip3'].values)
            plfile.write_record(pl['rotx'].values)
            plfile.write_record(pl['roty'].values)
            plfile.write_record(pl['rotz'].values)
        plfile.close()
        tpfile = FortranFile(param['TP_IN'], 'w')
        ntp = tp.id.count().values
        tpid = tp.id.values
        if param['IN_FORM'] == 'XV':
            v1 = tp['xhx'].values
            v2 = tp['xhy'].values
            v3 = tp['xhz'].values
            v4 = tp['vhx'].values
            v5 = tp['vhy'].values
            v6 = tp['vhz'].values
        elif param['IN_FORM'] == 'EL':
            v1 = tp['a'].values
            v2 = tp['e'].values
            v3 = tp['inc'].values
            v4 = tp['capom'].values
            v5 = tp['omega'].values
            v6 = tp['capm'].values
        else:
            print(f"{param['IN_FORM']} is not a valid input format type.")
        tpfile.write_record(ntp)
        tpfile.write_record(tpid)
        tpfile.write_record(v1)
        tpfile.write_record(v2)
        tpfile.write_record(v3)
        tpfile.write_record(v4)
        tpfile.write_record(v5)
        tpfile.write_record(v6)
    else:
        print(f"{param['IN_TYPE']} is an unknown file type")


def swifter_xr2infile(ds, param, framenum=-1):
    """
    Writes a set of Swifter input files from a single frame of a Swiftest xarray dataset

    Parameters
    ----------
    ds : xarray dataset
        Dataset containing Swifter n-body data in XV format
    framenum : int
        Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.
    param : dict
        Swifter input parameters. This method uses the names of the pl and tp files from the input

    Returns
    -------
    A set of input files for a Swifter run
    """
    frame = ds.isel(time=framenum)
    cb = frame.where(frame.id == 1, drop=True)
    pl = frame.where(frame.id > 1, drop=True)
    pl = pl.where(np.invert(np.isnan(pl['Gmass'])), drop=True).drop_vars(['J_2', 'J_4'])
    tp = frame.where(np.isnan(frame['Gmass']), drop=True).drop_vars(['Gmass', 'radius', 'J_2', 'J_4'])
    
    GMSun = np.double(cb['Gmass'])
    RSun = np.double(cb['radius'])
    param['J2'] = np.double(cb['J_2'])
    param['J4'] = np.double(cb['J_4'])
    
    if param['IN_TYPE'] == 'ASCII':
        # Swiftest Central body file
        plfile = open(param['PL_IN'], 'w')
        print(pl.id.count().values + 1, file=plfile)
        print(cb.id.values[0], GMSun, file=plfile)
        print('0.0 0.0 0.0', file=plfile)
        print('0.0 0.0 0.0', file=plfile)
        for i in pl.id:
            pli = pl.sel(id=i)
            if param['RHILL_PRESENT'] == "YES":
                print(i.values, pli['Gmass'].values, pli['rhill'].values, file=plfile)
            else:
                print(i.values, pli['Gmass'].values, file=plfile)
            if param['CHK_CLOSE'] == "YES":
                print(pli['radius'].values, file=plfile)
            print(pli['xhx'].values, pli['xhy'].values, pli['xhz'].values, file=plfile)
            print(pli['vhx'].values, pli['vhy'].values, pli['vhz'].values, file=plfile)
        plfile.close()
        
        # TP file
        tpfile = open(param['TP_IN'], 'w')
        print(tp.id.count().values, file=tpfile)
        for i in tp.id:
            tpi = tp.sel(id=i)
            print(i.values, file=tpfile)
            print(tpi['xhx'].values, tpi['xhy'].values, tpi['xhz'].values, file=tpfile)
            print(tpi['vhx'].values, tpi['vhy'].values, tpi['vhz'].values, file=tpfile)
        tpfile.close()
    else:
        # Now make Swiftest files
        print(f"{param['IN_TYPE']} is an unknown input file type")


def swift2swifter(swift_param, plname="", tpname="", conversion_questions={}):
    swifter_param = {}
    intxt = conversion_questions.get('RHILL', None)
    if not intxt:
        intxt = input("Is this a SyMBA input file with RHILL values in pl.in? (y/N)> ")
    if intxt.upper() == 'Y':
        isSyMBA = True
        swifter_param['RHILL_PRESENT'] = 'YES'
    else:
        isSyMBA = False
        swifter_param['RHILL_PRESENT'] = 'NO'
        
    isDouble = conversion_questions.get('DOUBLE', None)
    if not isDouble:
        print("Use single precision or double precision for real outputs?")
        print(" 1) Single (real*4)")
        print("*2) Double (real*8)")
        intxt = input("> ")
        if intxt == '1':
            isDouble = False
        else:
            isDouble = True
        
    # Convert the parameter file values
    swifter_param['T0'] = swift_param['T0']
    swifter_param['TSTOP'] = swift_param['TSTOP']
    swifter_param['DT'] = swift_param['DT']
    swifter_param['ISTEP_OUT'] = int(swift_param['DTOUT'] / swift_param['DT'])
    swifter_param['ISTEP_DUMP'] = int(swift_param['DTDUMP'] / swift_param['DT'])
    swifter_param['BIN_OUT'] = swift_param['BINARY_OUTPUTFILE']
    swifter_param['OUT_STAT'] = swift_param['STATUS_FLAG_FOR_OPEN_STATEMENTS']
    
    if swift_param['L5'] == "T":
        swifter_param['OUT_FORM'] = 'EL'
    else:
        swifter_param['OUT_FORM'] = 'XV'
    
    if swift_param['LCLOSE'] == "T":
        swifter_param['CHK_CLOSE'] = "YES"
    else:
        swifter_param['CHK_CLOSE'] = "NO"
        
    swifter_param['CHK_RMIN'] = swift_param['RMIN']
    swifter_param['CHK_RMAX'] = swift_param['RMAX']
    swifter_param['CHK_QMIN'] = swift_param['QMIN']
    if swift_param['QMIN'] != '-1':
        
        swifter_param['CHK_QMIN_COORD'] = conversion_questions.get('CHK_QMIN_COORD', None)
        if not swifter_param['CHK_QMIN_COORD']:
            print("CHK_QMIN_COORD value:")
            print("*1) HELIO")
            print(" 2) BARY")
            intxt = input("> ")
            if intxt == '2':
                swifter_param['CHK_QMIN_COORD'] = "BARY"
            else:
                swifter_param['CHK_QMIN_COORD'] = "HELIO"

        swifter_param['CHK_QMIN_RANGE'] = conversion_questions.get('CHK_QMIN_RANGE', None)
        if not swifter_param['CHK_QMIN_RANGE']:
            alo = input(f"Lower bound on CHK_QMIN_RANGE [{swift_param['RMIN']}]: ")
            if alo == '':
                alo = swift_param['RMIN']
            ahi = input(f"Upper bound on CHK_QMIN_RANGE: [{swift_param['RMAXU']}]: ")
            if ahi == '':
                ahi = swift_param['RMAXU']
            swifter_param['CHK_QMIN_RANGE'] = f"{alo} {ahi}"

    
    swifter_param['ENC_OUT'] = conversion_questions.get('ENC_OUT', None)
    if not swifter_param['ENC_OUT']:
        swifter_param['ENC_OUT'] = input("ENC_OUT: Encounter file name: [enc.dat]> ")
        if swifter_param['ENC_OUT'] == '':
            swifter_param['ENC_OUT'] = "enc.dat"

    intxt = conversion_questions.get('EXTRA_FORCE', None)
    if not intxt:
        intxt = input("EXTRA_FORCE: Use additional user-specified force routines? (y/N)> ")
    if intxt.upper() == 'Y':
        swifter_param['EXTRA_FORCE'] = 'YES'
    else:
        swifter_param['EXTRA_FORCE'] = 'NO'

    intxt = conversion_questions.get('BIG_DISCARD', None)
    if not intxt:
        intxt = input("BIG_DISCARD: include data for all bodies > GMTINY for each discard record? (y/N)> ")
    if intxt.upper() == 'Y':
        swifter_param['BIG_DISCARD'] = 'YES'
    else:
        swifter_param['BIG_DISCARD'] = 'NO'
    
    # Convert the PL file
    if plname == '':
        plname = input("PL_IN: Name of new planet input file: [pl.swifter.in]> ")
        if plname == '':
            plname = "pl.swifter.in"
    swifter_param['PL_IN'] = plname
    try:
        plnew = open(swifter_param['PL_IN'], 'w')
    except IOError:
        print(f"Cannot write to file {swifter_param['PL_IN']}")
        return swift_param
    
    print(f"Converting PL file: {swift_param['PL_IN']} -> {swifter_param['PL_IN']}")
    try:
        with open(swift_param['PL_IN'], 'r') as plold:
            line = plold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            npl = int(i_list[0])
            print(npl, file=plnew)
            line = plold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            GMcb = real2float(i_list[0])
            if swift_param['L1'] == "T":
                swifter_param['J2'] = real2float(i_list[1])
                swifter_param['J4'] = real2float(i_list[2])
            else:
                swifter_param['J2'] = 0.0
                swifter_param['J4'] = 0.0
            print(1, GMcb, file=plnew)
            print(0.0, 0.0, 0.0, file=plnew)
            print(0.0, 0.0, 0.0, file=plnew)
            line = plold.readline()  # Ignore the two zero vector lines
            line = plold.readline()
            for n in range(1, npl):  # Loop over planets
                line = plold.readline()
                i_list = [i for i in line.split(" ") if i.strip()]
                GMpl = real2float(i_list[0])
                if isSyMBA:
                    rhill = real2float(i_list[1])
                    if swift_param['LCLOSE'] == "T":
                        plrad = real2float(i_list[2])
                else:
                    if swift_param['LCLOSE'] == "T":
                        plrad = real2float(i_list[1])
                if swifter_param['RHILL_PRESENT'] == 'YES':
                    print(n + 1, GMpl, rhill, file=plnew)
                else:
                    print(n + 1, GMpl, file=plnew)
                if swifter_param['CHK_CLOSE'] == 'YES':
                    print(plrad, file=plnew)
                line = plold.readline()
                i_list = [i for i in line.split(" ") if i.strip()]
                xh = real2float(i_list[0])
                yh = real2float(i_list[1])
                zh = real2float(i_list[2])
                print(xh, yh, zh, file=plnew)
                line = plold.readline()
                i_list = [i for i in line.split(" ") if i.strip()]
                vhx = real2float(i_list[0])
                vhy = real2float(i_list[1])
                vhz = real2float(i_list[2])
                print(vhx, vhy, vhz, file=plnew)
        plnew.close()
        plold.close()
    except IOError:
        print(f"Error converting PL file")
        
    # Convert the TP file
    if tpname == '':
        tpname = input("TP_IN: Name of new test particle input file: [tp.swifter.in]> ")
        if tpname == '':
            tpname = "tp.swifter.in"
    swifter_param['TP_IN'] = tpname

    try:
        tpnew = open(swifter_param['TP_IN'], 'w')
    except IOError:
        print(f"Cannot write to file {swifter_param['TP_IN']}")

    print(f"Converting TP file: {swift_param['TP_IN']} -> {swifter_param['TP_IN']}")
    try:
        print(f'Writing out new TP file: {swifter_param["TP_IN"]}')
        with open(swift_param['TP_IN'], 'r') as tpold:
            line = tpold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            ntp = int(i_list[0])
            print(ntp, file=tpnew)
            for n in range(0, ntp):  # Loop over test particles
                print(npl + n + 1, file=tpnew)
                line = tpold.readline()
                i_list = [i for i in line.split(" ") if i.strip()]
                xh = real2float(i_list[0])
                yh = real2float(i_list[1])
                zh = real2float(i_list[2])
                print(xh, yh, zh, file=tpnew)
                line = tpold.readline()
                i_list = [i for i in line.split(" ") if i.strip()]
                vhx = real2float(i_list[0])
                vhy = real2float(i_list[1])
                vhz = real2float(i_list[2])
                print(vhx, vhy, vhz, file=tpnew)
                # Ignore STAT lines
                line = tpold.readline()
                line = tpold.readline()
    except IOError:
        print(f"Error converting TP file")
    swifter_param['! VERSION'] = "Swifter parameter file converted from Swift"
    
    return swifter_param

def swifter2swiftest(swifter_param, plname="", tpname="", cbname="", conversion_questions={}):
    swiftest_param = swifter_param.copy()
    # Pull additional feature status from the conversion_questions dictionary
    
    for key in newfeaturelist:
        swiftest_param[key] = conversion_questions.get(key, "NO")
    # Convert the PL file
    if plname == '':
        plname = input("PL_IN: Name of new planet input file: [pl.swiftest.in]> ")
        if plname == '':
            plname = "pl.swiftest.in"
    swiftest_param['PL_IN'] = plname
   
    try:
        plnew = open(swiftest_param['PL_IN'], 'w')
    except IOError:
        print(f"Cannot write to file {swiftest_param['PL_IN']}")
        return swifter_param

    print(f"Converting PL file: {swifter_param['PL_IN']} -> {swiftest_param['PL_IN']}")
    try:
        with open(swifter_param['PL_IN'], 'r') as plold:
            line = plold.readline()
            line = line.split("!")[0]  # Ignore comments
            i_list = [i for i in line.split(" ") if i.strip()]
            npl = int(i_list[0])
            print(npl - 1, file=plnew)
            line = plold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            GMcb = real2float(i_list[1])  # Store central body GM for later
            line = plold.readline()  # Ignore the two zero vector lines
            line = plold.readline()
            for n in range(1, npl):  # Loop over planets
                line = plold.readline()
                i_list = [i for i in line.split(" ") if i.strip()]
                idnum = int(i_list[0])
                GMpl = real2float(i_list[1])
                if swifter_param['RHILL_PRESENT'] == 'YES':
                   rhill = real2float(i_list[2])
                   print(idnum, GMpl, rhill, file=plnew)
                else:
                   print(idnum, GMpl, file=plnew)
                if swifter_param['CHK_CLOSE'] == 'YES':
                    line = plold.readline()
                    i_list = [i for i in line.split(" ") if i.strip()]
                    plrad = real2float(i_list[0])
                    print(plrad, file=plnew)
                line = plold.readline()
                i_list = [i for i in line.split(" ") if i.strip()]
                xh = real2float(i_list[0])
                yh = real2float(i_list[1])
                zh = real2float(i_list[2])
                print(xh, yh, zh, file=plnew)
                line = plold.readline()
                i_list = [i for i in line.split(" ") if i.strip()]
                vhx = real2float(i_list[0])
                vhy = real2float(i_list[1])
                vhz = real2float(i_list[2])
                print(vhx, vhy, vhz, file=plnew)
        plnew.close()
        plold.close()
    except IOError:
        print(f"Error converting PL file")
        
    # Convert the TP file
    if tpname == '':
        tpname = input("TP_IN: Name of new planet input file: [tp.swiftest.in]> ")
        if tpname == '':
            tpname = "tp.swiftest.in"
    swiftest_param['TP_IN'] = tpname

    try:
        tpnew = open(swiftest_param['TP_IN'], 'w')
    except IOError:
        print(f"Cannot write to file {swiftest_param['TP_IN']}")
        return swifter_param
  
    print(f"Converting TP file: {swifter_param['TP_IN']} -> {swiftest_param['TP_IN']}")
    try:
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
                xh = real2float(i_list[0])
                yh = real2float(i_list[1])
                zh = real2float(i_list[2])
                print(xh, yh, zh, file=tpnew)
                line = tpold.readline()
                i_list = [i for i in line.split(" ") if i.strip()]
                vhx = real2float(i_list[0])
                vhy = real2float(i_list[1])
                vhz = real2float(i_list[2])
                print(vhx, vhy, vhz, file=tpnew)
        
        tpold.close()
        tpnew.close()
    except IOError:
        print(f"Error converting TP file")
        
    # Create the CB file
    if cbname == '':
        cbname = input("CB_IN: Name of new planet input file: [cb.swiftest.in]> ")
        if cbname == '':
            cbname = "cb.swiftest.in"
    swiftest_param['CB_IN'] = cbname
   
    unit_system = conversion_questions.get('UNITS', '')
    unit_type = 5
    if not unit_system:
        print(f"\nCentral body G*M = {GMcb}\n")
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
    if unit_type == 1 or unit_system.upper() == 'MSUN-AU-YEAR':
        print("Unit system is MSun-AU-year")
        swiftest_param['MU2KG'] = swiftest.MSun
        swiftest_param['DU2M'] = swiftest.AU2M
        swiftest_param['TU2S'] = swiftest.YR2S
    elif unit_type == 2 or unit_system.upper() == 'MSUN-AU-DAY':
        print("Unit system is MSun-AU-day")
        swiftest_param['MU2KG'] = swiftest.MSun
        swiftest_param['DU2M'] = swiftest.AU2M
        swiftest_param['TU2S'] = swiftest.JD2S
    elif unit_type == 3 or unit_system.upper() == 'SI':
        print("Unit system is SI: kg-m-s")
        swiftest_param['MU2KG'] = 1.0
        swiftest_param['DU2M'] = 1.0
        swiftest_param['TU2S'] = 1.0
    elif unit_type == 4 or unit_system.upper() == 'CGS':
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
    GU = swiftest.GC / (swiftest_param['DU2M'] ** 3 / (swiftest_param['MU2KG'] * swiftest_param['TU2S'] ** 2))
    print(f"Central body mass: {GMcb / GU} MU ({(GMcb / GU) * swiftest_param['MU2KG']} kg)")
   
    cbrad = conversion_questions.get('CBRAD', None)
    if not cbrad:
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
            cbrad = real2float(cbrad.strip())
    
    print(f'Writing out new CB file: {swiftest_param["CB_IN"]}')
    # Write out new central body file
    try:
        cbnew = open(swiftest_param['CB_IN'], 'w')
        print(GMcb, file=cbnew)
        print(cbrad, file=cbnew)
        print(swifter_param['J2'], file=cbnew)
        print(swifter_param['J4'], file=cbnew)
        cbnew.close()
    except IOError:
        print(f"Cannot write to file {swiftest_param['CB_IN']}")
        return swifter_param
   
    GMTINY = conversion_questions.get('GMTINY', None)
    if not GMTINY:
        GMTINY = input(f"Value of GMTINY if this is a SyMBA simulation (enter nothing if this is not a SyMBA parameter file)> ")
    if GMTINY != '' and real2float(GMTINY.strip()) > 0:
        swiftest_param['GMTINY'] = real2float(GMTINY.strip())
        
    # Remove the unneeded parameters
    if 'C' in swiftest_param:
        swiftest_param['GR'] = 'YES'
        swiftest_param.pop('C', None)
    swiftest_param.pop('J2', None)
    swiftest_param.pop('J4', None)

    swiftest_param['DISCARD_OUT'] = conversion_questions.get('DISCARD_OUT', '')
    if not swiftest_param['DISCARD_OUT']:
        swiftest_param['DISCARD_OUT'] = input("DISCARD_OUT: Discard file name: [discard.out]> ")
        if swiftest_param['DISCARD_OUT'] == '':
            swiftest_param['DISCARD_OUT'] = "discard.out"

    swiftest_param['PARTICLE_OUT'] = conversion_questions.get('PARTICLE_OUT', '')
    if not swiftest_param['PARTICLE_OUT']:
        swiftest_param['PARTICLE_OUT'] = input("PARTICLE_OUT: Particle info file name (Only used by SyMBA): []> ")

    swiftest_param['! VERSION'] = "Swiftest parameter file converted from Swifter"
    return swiftest_param

def swift2swiftest(swift_param, plname="", tpname="", cbname="", conversion_questions={}):
    if plname == '':
        plname = input("PL_IN: Name of new planet input file: [pl.swiftest.in]> ")
        if plname == '':
            plname = "pl.swiftest.in"
    pltmp = tempfile.NamedTemporaryFile()
    pltmpname = pltmp.name

    if tpname == '':
        tpname = input("TP_IN: Name of new planet input file: [tp.swiftest.in]> ")
        if tpname == '':
            tpname = "tp.swiftest.in"
    tptmp = tempfile.NamedTemporaryFile()
    tptmpname = tptmp.name

    # Create the CB file
    if cbname == '':
        cbname = input("CB_IN: Name of new planet input file: [cb.swiftest.in]> ")
        if cbname == '':
            cbname = "cb.swiftest.in"
    swifter_param = swift2swifter(swift_param, pltmpname, tptmpname, conversion_questions)
    swiftest_param = swifter2swiftest(swifter_param, plname, tpname, cbname, conversion_questions)
    swiftest_param['! VERSION'] = "Swiftest parameter file converted from Swift"
    return swiftest_param

def swiftest2swifter_param(swiftest_param, J2=0.0, J4=0.0):
    swifter_param = swiftest_param
    CBIN = swifter_param.pop("CB_IN", None)
    GMTINY = swifter_param.pop("GMTINY", None)
    DISCARD_OUT = swifter_param.pop("DISCARD_OUT", None)
    MU2KG = swifter_param.pop("MU2KG", 1.0)
    DU2M = swifter_param.pop("DU2M", 1.0)
    TU2S = swifter_param.pop("TU2S", 1.0)
    GR = swifter_param.pop("GR", None)
    if GR is not None:
        if GR == 'YES':
           swifter_param['C'] =  swiftest.einsteinC * np.longdouble(TU2S) / np.longdouble(DU2M)
    for key in newfeaturelist:
       tmp = swifter_param.pop(key, None)
    swifter_param['J2'] = J2
    swifter_param['J4'] = J4
    if swifter_param['OUT_STAT'] == "REPLACE":
        swifter_param['OUT_STAT'] = "UNKNOWN"
    swifter_param['! VERSION'] = "Swifter parameter file converted from Swiftest"

    return swifter_param