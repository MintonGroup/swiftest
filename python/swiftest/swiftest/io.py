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
from scipy.io import FortranFile
import xarray as xr
import sys
import tempfile
import re

# This defines features that are new in Swiftest and not in Swifter (for conversion between param.in files)
newfeaturelist = ("RESTART",
                  "FRAGMENTATION",
                  "ROTATION",
                  "TIDES",
                  "ENERGY",
                  "GR",
                  "YARKOVSKY",
                  "YORP",
                  "IN_FORM",
                  "SEED",
                  "INTERACTION_LOOPS",
                  "ENCOUNTER_CHECK",
                  "TSTART",
                  "DUMP_CADENCE")

# This list defines features that are booleans, so must be converted to/from string when writing/reading from file
bool_param = ["RESTART",
              "CHK_CLOSE",
              "EXTRA_FORCE",
              "RHILL_PRESENT",
              "BIG_DISCARD",
              "FRAGMENTATION",
              "ROTATION",
              "TIDES",
              "ENERGY",
              "GR",
              "YARKOVSKY",
              "YORP"]

int_param = ["ISTEP_OUT", "DUMP_CADENCE"]
float_param = ["T0", "TSTART", "TSTOP", "DT", "CHK_RMIN", "CHK_RMAX", "CHK_EJECT", "CHK_QMIN", "DU2M", "MU2KG",
               "TU2S", "MIN_GMFRAG", "GMTINY"]

upper_str_param = ["OUT_TYPE","OUT_FORM","OUT_STAT","IN_TYPE","IN_FORM"]

# This defines Xarray Dataset variables that are strings, which must be processed due to quirks in how NetCDF-Fortran
# handles strings differently than Python's Xarray.
string_varnames = ["name", "particle_type", "status", "origin_type"]
char_varnames = ["space"]
int_varnames = ["id", "ntp", "npl", "nplm", "discard_body_id", "collision_id"]

def bool2yesno(boolval):
    """
    Converts a boolean into a string of either "YES" or "NO".

    Parameters
    ----------
    boolval : bool
        Input value

    Returns
    -------
    {"YES","NO"}

    """
    if boolval:
        return "YES"
    else:
        return "NO"

def bool2tf(boolval):
    """
    Converts a boolean into a string of either "T" or "F".

    Parameters
    ----------
    boolval : bool
        Input value

    Returns
    -------
    {"T","F"}

    """
    if boolval:
        return "T"
    else:
        return "F"

def str2bool(input_str):
    """
    Converts a string into an equivalent boolean.

    Parameters
    ----------
    input_str : {"YES", "Y", "T", "TRUE", ".TRUE.", "NO", "N", "F", "FALSE", ".FALSE."}
        Input string. Input is case-insensitive.

    Returns
    -------
    {True, False}

    """
    if type(input_str) is bool:
       return input_str
    valid_true = ["YES", "Y", "T", "TRUE", ".TRUE."]
    valid_false = ["NO", "N", "F", "FALSE", ".FALSE."]
    if input_str.upper() in valid_true:
        return True
    elif input_str.upper() in valid_false:
        return False
    else:
        raise ValueError(f"{input_str} is not recognized as boolean")



def real2float(realstr):
    """
    Converts a Fortran-generated ASCII string of a real value into a numpy float type. Handles cases where double precision
    numbers in exponential notation use 'd' or 'D' instead of 'e' or 'E'
    
    Parameters
    ----------
    realstr : str
        Fortran-generated ASCII string of a real value.

    Returns
    -------
     : float
        The converted floating point value of the input string
    """
    return float(realstr.replace('d', 'E').replace('D', 'E'))


def read_swiftest_param(param_file_name, param, verbose=True):
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
    param['! VERSION'] = f"Swiftest parameter input file"

    # Read param.in file
    if verbose: print(f'Reading Swiftest file {param_file_name}')
    try:
        with open(param_file_name, 'r') as f:
            for line in f.readlines():
                fields = line.split()
                if len(fields) > 0:
                    if fields[0][0] != '!':
                        key = fields[0].upper()
                        param[key] = fields[1]
                    # Special case of CHK_QMIN_RANGE requires a second input
                    if fields[0].upper() == 'CHK_QMIN_RANGE':
                        alo = real2float(fields[1])
                        ahi = real2float(fields[2])
                        param['CHK_QMIN_RANGE'] = f"{alo} {ahi}"

        for uc in upper_str_param:
            if uc in param:
                param[uc] = param[uc].upper()

        for i in int_param:
            if i in param and type(param[i]) != int:
                param[i] = int(float(param[i]))

        for f in float_param:
            if f in param and type(param[f]) is str:
                param[f] = real2float(param[f])

        for b in bool_param:
            if b in param:
                param[b] = str2bool(param[b])
    except IOError:
        print(f"{param_file_name} not found.")
    return param

def reorder_dims(ds):

    # Re-order dimension coordinates so that they are in the same order as the Fortran side
    idx = ds.indexes
    if "id" in idx:
        dim_order = ["time", "space", "id"]
    elif "name" in idx:
        dim_order = ["time", "space", "name"]
    else:
        dim_order = idx
    idx = {index_name: idx[index_name] for index_name in dim_order}
    ds = ds.reindex(idx)
    return ds
def read_swifter_param(param_file_name, verbose=True):
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
        'CHK_CLOSE': False,
        'CHK_RMIN': "-1.0",
        'CHK_RMAX': "-1.0",
        'CHK_EJECT': "-1.0",
        'CHK_QMIN': "-1.0",
        'CHK_QMIN_COORD': "HELIO",
        'CHK_QMIN_RANGE': "",
        'ENC_OUT': "",
        'EXTRA_FORCE': False,
        'BIG_DISCARD': False,
        'RHILL_PRESENT': False,
        'C': "-1.0",
    }
    
    # Read param.in file
    if verbose: print(f'Reading Swifter file {param_file_name}')
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
        for b in bool_param:
            if b in param:
                param[b] = str2bool(param[b])
        if param['C'] != '-1.0':
            param['C'] = real2float(param['C'])
        else:
            param.pop('C', None)
    except IOError:
        print(f"{param_file_name} not found.")

    return param


def read_swift_param(param_file_name, startfile="swift.in", verbose=True):
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
    if verbose: print(f'Reading Swift file {param_file_name}')
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
    """
    Writes a Swift param.in file.

    Parameters
    ----------
    param : dictionary
        The entries in the user parameter file
    param_file_name : string
        File name of the input parameter file

    Returns
    -------
    Prints a text file containing the parameter information.
    """
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
    """
    Writes a Swifter/Swiftest param.in file.

    Parameters
    ----------
    param : dictionary
        The entries in the user parameter file
    param_file_name : string
        File name of the input parameter file

    Returns
    -------
    Prints a text file containing the parameter information.
    """
    outfile = open(param_file_name, 'w')
    keylist = ['! VERSION',
               'T0',
               'TSTART',
               'TSTOP',
               'DT',
               'ISTEP_OUT',
               'DUMP_CADENCE',
               'NC_IN',
               'PL_IN',
               'TP_IN',
               'CB_IN',
               'IN_TYPE',
               'IN_FORM',
               'BIN_OUT',
               'OUT_FORM',
               'OUT_TYPE',
               'OUT_STAT',
               'CHK_QMIN',
               'CHK_RMIN',
               'CHK_RMAX',
               'CHK_EJECT',
               'CHK_QMIN_COORD',
               'CHK_QMIN_RANGE',
               'MU2KG',
               'TU2S',
               'DU2M',
               'GMTINY',
               'FRAGMENTATION',
               'MIN_GMFRAG',
               'RESTART']
    ptmp = param.copy()
    # Print the list of key/value pairs in the preferred order
    for key in keylist:
        val = ptmp.pop(key, None)
        if val is not None:
            if type(val) is bool:
                print(f"{key:<16} {bool2yesno(val)}", file=outfile)
            else:
                print(f"{key:<16} {val}", file=outfile)
    # Print the remaining key/value pairs in whatever order
    for key, val in ptmp.items():
        if val != "":
            if type(val) is bool:
                print(f"{key:<16} {bool2yesno(val)}", file=outfile)
            else:
                print(f"{key:<16} {val}", file=outfile)
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
    """
    Creates the lables for the variables to be included in the output file.

    Parameters
    ----------
    param : dictionary
        The entries in the user parameter file

    Returns
    -------
    clab : string list
        Labels for the cvec data
    plab : string list
        Labels for the pvec data
    tlab : string list
        Labels for the tvec data
    infolab_float : string list
        Labels for floating point data
    infolab_int : 
        Labels for integer data
    infolab_str :   
        Labels for string data
    """
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
    if param['CHK_CLOSE']:
       plab.append('radius')
    if param['RHILL_PRESENT']:
        plab.append('rhill')
    clab = ['Gmass', 'radius', 'j2rp2', 'j4rp4']
    if param['ROTATION']:
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
    if param['TIDES']:
        clab.append('k2')
        clab.append('Q')
        plab.append('k2')
        plab.append('Q')

    infolab_float = [
        "origin_time",
        "origin_xhx",
        "origin_xhy",
        "origin_xhz",
        "origin_vhx",
        "origin_vhy",
        "origin_vhz",
        "discard_time",
        "discard_xhx",
        "discard_xhy",
        "discard_xhz",
        "discard_vhx",
        "discard_vhy",
        "discard_vhz",
        ]
    infolab_int = [
        "collision_id",
        "discard_body_id"
        ]
    infolab_str = [
        "particle_type",
        "origin_type",
        "status",
        ]

    return clab, plab, tlab, infolab_float, infolab_int, infolab_str


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
        if param['ROTATION']:
            Ipcbx = f.read_reals(np.float64)
            Ipcby = f.read_reals(np.float64)
            Ipcbz = f.read_reals(np.float64)
            rotcbx = f.read_reals(np.float64)
            rotcby = f.read_reals(np.float64)
            rotcbz = f.read_reals(np.float64)
        if param['TIDES']:
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
            if param['RHILL_PRESENT']:
                rhill = f.read_reals(np.float64)
            Rpl = f.read_reals(np.float64)
            if param['ROTATION']:
                Ipplx = f.read_reals(np.float64)
                Ipply = f.read_reals(np.float64)
                Ipplz = f.read_reals(np.float64)
                rotplx = f.read_reals(np.float64)
                rotply = f.read_reals(np.float64)
                rotplz = f.read_reals(np.float64)
            if param['TIDES']:
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
        
        clab, plab, tlab, infolab_float, infolab_int, infolab_str = make_swiftest_labels(param)

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
        if param['RHILL_PRESENT']:
           if npl > 0:
              pvec = np.vstack([pvec, rhill])
        if param['ROTATION']:
            cvec = np.vstack([cvec, Ipcbx, Ipcby, Ipcbz, rotcbx, rotcby, rotcbz])
            if npl > 0:
                pvec = np.vstack([pvec, Ipplx, Ipply, Ipplz, rotplx, rotply, rotplz])
        if param['TIDES']:
            cvec = np.vstack([cvec, k2cb, Qcb])
            if npl > 0:
                pvec = np.vstack([pvec, k2pl, Qpl])
        yield t, cbid, cbnames, cvec.T, clab, \
              npl, plid, plnames, pvec.T, plab, \
              ntp, tpid, tpnames, tvec.T, tlab


def swifter2xr(param, verbose=True):
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
    dims = ['time', 'id','vec']
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
        if verbose: print('\nCreating Dataset')
        ds = xr.combine_by_coords([plds, tpds])
        if verbose: print(f"Successfully converted {ds.sizes['time']} output frames.")
    return ds


def swiftest2xr(param, verbose=True):
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

    if ((param['OUT_TYPE'] == 'NETCDF_DOUBLE') or (param['OUT_TYPE'] == 'NETCDF_FLOAT')):
        if verbose: print('\nCreating Dataset from NetCDF file')
        ds = xr.open_dataset(param['BIN_OUT'], mask_and_scale=False)
        if param['OUT_TYPE'] == "NETCDF_DOUBLE":
            ds = fix_types(ds,ftype=np.float64)
        elif param['OUT_TYPE'] == "NETCDF_FLOAT":
            ds = fix_types(ds,ftype=np.float32)
        # Check if the name variable contains unique values. If so, make name the dimension instead of id
        if len(np.unique(ds['name'])) == len(ds['name']):
           ds = ds.swap_dims({"id" : "name"})
           ds = ds.reset_coords("id")
    else:
        print(f"Error encountered. OUT_TYPE {param['OUT_TYPE']} not recognized.")
        return None
    if verbose: print(f"Successfully converted {ds.sizes['time']} output frames.")

    return ds

def xstrip(a):
    """
    Cleans up the string values in the DataSet to remove extra white space

    Parameters
    ----------
    a    : xarray dataset
 
    Returns
    -------
    da : xarray dataset with the strings cleaned up
    """
    func = lambda x: np.char.strip(x)
    return xr.apply_ufunc(func, a.str.decode(encoding='utf-8'),dask='parallelized')

def string_converter(da):
    """
    Converts a string to a unicode string 

    Parameters
    ----------
    da    : xarray dataset

    Returns
    -------
    da : xarray dataset with the strings cleaned up
    """
    if da.dtype == np.dtype(object):
       da = da.astype('<U32')
    elif type(da.values[0]) != np.str_:
       da = xstrip(da)
    return da

def char_converter(da):
    """
    Converts a string to a unicode string

    Parameters
    ----------
    da    : xarray dataset

    Returns
    -------
    da : xarray dataset with the strings cleaned up
    """
    if da.dtype == np.dtype(object):
        da = da.astype('<U1')
    elif type(da.values[0]) != np.str_:
        da = xstrip(da)
    return da

def clean_string_values(ds):
    """
    Cleans up the string values in the DataSet that have artifacts as a result of coming from NetCDF Fortran

    Parameters
    ----------
    ds    : xarray dataset
 
    Returns
    -------
    ds : xarray dataset with the strings cleaned up
    """

    for n in string_varnames:
        if n in ds:
           ds[n] = string_converter(ds[n])

    for n in char_varnames:
        if n in ds:
            ds[n] = char_converter(ds[n])
    return ds


def unclean_string_values(ds):
    """
    Returns strings back to a format readable to NetCDF Fortran

    Parameters
    ----------
    ds    : xarray dataset

    Returns
    -------
    ds : xarray dataset with the strings cleaned up
    """

    for c in string_varnames:
        if c in ds:
            n = string_converter(ds[c])
            ds[c] = n.str.ljust(32).str.encode('utf-8')

    for c in char_varnames:
        if c in ds:
            n = string_converter(ds[c])
            ds[c] = n.str.ljust(1).str.encode('utf-8')
    return ds

def fix_types(ds,itype=np.int64,ftype=np.float64):

    ds = clean_string_values(ds)
    for intvar in int_varnames:
        if intvar in ds:
            ds[intvar] = ds[intvar].astype(itype)

    float_varnames = [x for x in list(ds.keys()) if x not in string_varnames + int_varnames + char_varnames]

    for floatvar in float_varnames:
        ds[floatvar] = ds[floatvar].astype(ftype)

    float_coordnames = [x for x in list(ds.coords) if x not in string_varnames + int_varnames + char_varnames]
    for floatcoord in float_coordnames:
        ds[floatcoord] = ds[floatcoord].astype(np.float64)


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
    """
    Reads in the Swiftest SyMBA-generated PARTICLE_OUT  and converts it to an xarray Dataset

    Parameters
    ----------
    param : dict
        Swiftest parameters

    Returns
    -------
    infoxr : xarray dataset
    """
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

def select_active_from_frame(ds, param, framenum=-1):
    """
    Selects a particular frame from a DataSet and returns only the active particles in that frame

    Parameters
    ----------
    ds : xarray dataset
        Dataset containing Swiftest n-body data
    param : dict
        Swiftest input parameters. This method uses the names of the cb, pl, and tp files from the input
    framenum : int (default=-1)
        Time frame to extract. If this argument is not passed, the default is to use the last frame in the dataset.

    Returns
    -------
    Dataset containing only the active particles at frame, framenum
    """


    # Select the input time frame from the Dataset
    frame = ds.isel(time=[framenum])
    iframe = frame.isel(time=0)

    if "name" in ds.dims:
        count_dim = "name"
    elif "id" in ds.dims:
        count_dim = "id"

    # Select only the active particles at this time step
    # Remove the inactive particles
    if param['OUT_FORM'] == 'XV' or param['OUT_FORM'] == 'XVEL':
        iactive = iframe[count_dim].where((~np.isnan(iframe['Gmass'])) | (~np.isnan(iframe['xhx'])), drop=True)[count_dim]
    else:
        iactive = iframe[count_dim].where((~np.isnan(iframe['Gmass'])) | (~np.isnan(iframe['a'])), drop = True)[count_dim]
    if count_dim == "id":
        frame = frame.sel(id=iactive.values)
    elif count_dim == "name":
        frame = frame.sel(name=iactive.values)


    return frame

def swiftest_xr2infile(ds, param, in_type="NETCDF_DOUBLE", infile_name=None,framenum=-1,verbose=True):
    """
    Writes a set of Swiftest input files from a single frame of a Swiftest xarray dataset

    Parameters
    ----------
    ds : xarray dataset
        Dataset containing Swiftest n-body data in XV format
    param : dict
        Swiftest input parameters. This method uses the names of the cb, pl, and tp files from the input
    framenum : int (default=-1)
        Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.

    Returns
    -------
    A set of input files for a new Swiftest run
    """
    param_tmp = param.copy()
    param_tmp['OUT_FORM'] = param['IN_FORM']
    frame = select_active_from_frame(ds, param_tmp, framenum)

    if "name" in frame.dims:
        frame = frame.swap_dims({"name" : "id"})
        frame = frame.reset_coords("name")

    if in_type == "NETCDF_DOUBLE" or in_type == "NETCDF_FLOAT":
        # Convert strings back to byte form and save the NetCDF file
        # Note: xarray will call the character array dimension string32. The Fortran code
        # will rename this after reading

        if infile_name is None:
            infile_name = param['NC_IN']
        frame = unclean_string_values(frame)
        if verbose:
            print(f"Writing initial conditions to file {infile_name}")
        frame = reorder_dims(frame)
        frame.to_netcdf(path=infile_name)
        return frame

    # All other file types need seperate files for each of the inputs
    cb = frame.where(frame.id == 0, drop=True)
    pl = frame.where(frame.id > 0, drop=True)
    pl = pl.where(np.invert(np.isnan(pl['Gmass'])), drop=True).drop_vars(['j2rp2', 'j2rp2'],errors="ignore")
    tp = frame.where(np.isnan(frame['Gmass']), drop=True).drop_vars(['Gmass', 'radius', 'j2rp2', 'j4rp4'],errors="ignore")
    
    GMSun = np.double(cb['Gmass'])
    if param['CHK_CLOSE']:
       RSun = np.double(cb['radius'])
    else:
       RSun = param['CHK_RMIN']
    J2 = np.double(cb['j2rp2'])
    J4 = np.double(cb['j4rp4'])
    cbname = cb['name'].values[0]
    if param['ROTATION']:
        Ip1cb = np.double(cb['Ip1'])
        Ip2cb = np.double(cb['Ip2'])
        Ip3cb = np.double(cb['Ip3'])
        rotxcb = np.double(cb['rotx'])
        rotycb = np.double(cb['roty'])
        rotzcb = np.double(cb['rotz'])
    cbid = int(0)
    
    if in_type == 'ASCII':
        # Swiftest Central body file
        cbfile = open(param['CB_IN'], 'w')
        print(cbname, file=cbfile)
        print(GMSun, file=cbfile)
        print(RSun, file=cbfile)
        print(J2, file=cbfile)
        print(J4, file=cbfile)
        if param['ROTATION']:
            print(Ip1cb, Ip2cb, Ip3cb, file=cbfile)
            print(rotxcb, rotycb, rotzcb, file=cbfile)
        cbfile.close()
        
        plfile = open(param['PL_IN'], 'w')
        print(pl.id.count().values, file=plfile)
        for i in pl.id:
            pli = pl.sel(id=i)
            if param['RHILL_PRESENT']:
               print(pli['name'].values[0], pli['Gmass'].values[0], pli['rhill'].values[0], file=plfile)
            else:
               print(pli['name'].values[0], pli['Gmass'].values[0], file=plfile)
            if param['CHK_CLOSE']:
               print(pli['radius'].values[0], file=plfile)
            if param['IN_FORM'] == 'XV':
                print(pli['xhx'].values[0], pli['xhy'].values[0], pli['xhz'].values[0], file=plfile)
                print(pli['vhx'].values[0], pli['vhy'].values[0], pli['vhz'].values[0], file=plfile)
            elif param['IN_FORM'] == 'EL':
                print(pli['a'].values[0], pli['e'].values[0], pli['inc'].values[0], file=plfile)
                print(pli['capom'].values[0], pli['omega'].values[0], pli['capm'].values[0], file=plfile)
            else:
                print(f"{param['IN_FORM']} is not a valid input format type.")
            if param['ROTATION']:
                print(pli['Ip1'].values[0], pli['Ip2'].values[0], pli['Ip3'].values[0], file=plfile)
                print(pli['rotx'].values[0], pli['roty'].values[0], pli['rotz'].values[0], file=plfile)
        plfile.close()
        
        # TP file
        tpfile = open(param['TP_IN'], 'w')
        print(tp.id.count().values, file=tpfile)
        for i in tp.id:
            tpi = tp.sel(id=i)
            print(tpi['name'].values[0], file=tpfile)
            if param['IN_FORM'] == 'XV':
                print(tpi['xhx'].values[0], tpi['xhy'].values[0], tpi['xhz'].values[0], file=tpfile)
                print(tpi['vhx'].values[0], tpi['vhy'].values[0], tpi['vhz'].values[0], file=tpfile)
            elif param['IN_FORM'] == 'EL':
                print(tpi['a'].values[0], tpi['e'].values[0], tpi['inc'].values[0], file=tpfile)
                print(tpi['capom'].values[0], tpi['omega'].values[0], tpi['capm'].values[0], file=tpfile)
            else:
                print(f"{param['IN_FORM']} is not a valid input format type.")
        tpfile.close()
    else:
        print(f"{in_type} is an unknown file type")


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

    if "name" in frame.dims:
        frame = frame.swap_dims({"name" : "id"})
        frame = frame.reset_coords("name")

    cb = frame.where(frame.id == 0, drop=True)
    pl = frame.where(frame.id > 0, drop=True)
    pl = pl.where(np.invert(np.isnan(pl['Gmass'])), drop=True).drop_vars(['j2rp2', 'j4rp4'])
    tp = frame.where(np.isnan(frame['Gmass']), drop=True).drop_vars(['Gmass', 'radius', 'j2rp2', 'j4rp4'])
    
    GMSun = np.double(cb['Gmass'])
    if param['CHK_CLOSE']:
       RSun = np.double(cb['radius'])
    else:
       RSun = param['CHK_RMIN']
    param['J2'] = np.double(cb['j2rp2'])
    param['J4'] = np.double(cb['j4rp4'])
    
    if param['IN_TYPE'] == 'ASCII':
        # Swiftest Central body file
        plfile = open(param['PL_IN'], 'w')
        print(pl.id.count().values + 1, file=plfile)
        print(cb.id.values[0], GMSun, file=plfile)
        print('0.0 0.0 0.0', file=plfile)
        print('0.0 0.0 0.0', file=plfile)
        for i in pl.id:
            pli = pl.sel(id=i)
            if param['RHILL_PRESENT']:
                print(i.values, pli['Gmass'].values, pli['rhill'].values, file=plfile)
            else:
                print(i.values, pli['Gmass'].values, file=plfile)
            if param['CHK_CLOSE']:
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

    return

def swift2swifter(swift_param, plname="", tpname="", conversion_questions={}):
    """
    Converts from a Swift run to a Swifter run

    Parameters
    ----------
    swift_param : dictionary
       Swift input parameters. 
    plname : string
        Name of massive body input file
    tpname : string
        Name of test particle input file
    conversion_questions : dictronary
        Dictionary of additional parameters required to convert between formats

    Returns
    -------
    swifter_param : A set of input files for a new Swifter run
    """

    swifter_param = {}
    intxt = conversion_questions.get('RHILL', None)
    if not intxt:
        intxt = input("Is this a SyMBA input file with RHILL values in pl.in? (y/N)> ")
    if intxt.upper() == 'Y':
        isSyMBA = True
        swifter_param['RHILL_PRESENT'] = True
    else:
        isSyMBA = False
        swifter_param['RHILL_PRESENT'] = False
        
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
        swifter_param['CHK_CLOSE'] = True
    else:
        swifter_param['CHK_CLOSE'] = False
        
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
        swifter_param['EXTRA_FORCE'] = True
    else:
        swifter_param['EXTRA_FORCE'] = False

    intxt = conversion_questions.get('BIG_DISCARD', None)
    if not intxt:
        intxt = input("BIG_DISCARD: include data for all bodies > GMTINY for each discard record? (y/N)> ")
    if intxt.upper() == 'Y':
        swifter_param['BIG_DISCARD'] = True
    else:
        swifter_param['BIG_DISCARD'] = False
    
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
            i_list = [i for i in re.split('  +|\t',line) if i.strip()]
            npl = int(i_list[0])
            print(npl, file=plnew)
            line = plold.readline()
            i_list = [i for i in re.split('  +|\t',line) if i.strip()]
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
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                GMpl = real2float(i_list[0])
                if isSyMBA:
                    rhill = real2float(i_list[1])
                    if swift_param['LCLOSE'] == "T":
                        plrad = real2float(i_list[2])
                else:
                    if swift_param['LCLOSE'] == "T":
                        plrad = real2float(i_list[1])
                if swifter_param['RHILL_PRESENT']:
                    print(n + 1, GMpl, rhill, file=plnew)
                else:
                    print(n + 1, GMpl, file=plnew)
                if swifter_param['CHK_CLOSE']:
                    print(plrad, file=plnew)
                line = plold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                xh = real2float(i_list[0])
                yh = real2float(i_list[1])
                zh = real2float(i_list[2])
                print(xh, yh, zh, file=plnew)
                line = plold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
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
            i_list = [i for i in re.split('  +|\t',line) if i.strip()]
            ntp = int(i_list[0])
            print(ntp, file=tpnew)
            for n in range(0, ntp):  # Loop over test particles
                print(npl + n + 1, file=tpnew)
                line = tpold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                xh = real2float(i_list[0])
                yh = real2float(i_list[1])
                zh = real2float(i_list[2])
                print(xh, yh, zh, file=tpnew)
                line = tpold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
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
    """
    Converts from a Swifter run to a Swiftest run

    Parameters
    ----------
    swifter_param : dictionary
       Swifter input parameters. 
    plname : string
        Name of massive body input file
    tpname : string
        Name of test particle input file
    cbname : string
        Name of central body input file
    conversion_questions : dictronary
        Dictionary of additional parameters required to convert between formats

    Returns
    -------
    swiftest_param : A set of input files for a new Swiftest run
    """

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
            i_list = [i for i in re.split('  +|\t',line) if i.strip()]
            npl = int(i_list[0])
            print(npl - 1, file=plnew)
            line = plold.readline()
            i_list = [i for i in re.split('  +|\t',line) if i.strip()]
            GMcb = real2float(i_list[1])  # Store central body GM for later
            line = plold.readline()  # Ignore the two zero vector lines
            line = plold.readline()
            for n in range(1, npl):  # Loop over planets
                line = plold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                idnum = int(i_list[0])
                GMpl = real2float(i_list[1])
                if swifter_param['RHILL_PRESENT']:
                   rhill = real2float(i_list[2])
                   print(idnum, GMpl, rhill, file=plnew)
                else:
                   print(idnum, GMpl, file=plnew)
                if swifter_param['CHK_CLOSE']:
                    line = plold.readline()
                    i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                    plrad = real2float(i_list[0])
                    print(plrad, file=plnew)
                line = plold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                xh = real2float(i_list[0])
                yh = real2float(i_list[1])
                zh = real2float(i_list[2])
                print(xh, yh, zh, file=plnew)
                line = plold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
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
            i_list = [i for i in re.split('  +|\t',line) if i.strip()]
            ntp = int(i_list[0])
            print(ntp, file=tpnew)
            for n in range(0, ntp):  # Loop over test particles
                line = tpold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                name = int(i_list[0])
                print(name, file=tpnew)
                line = tpold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                xh = real2float(i_list[0])
                yh = real2float(i_list[1])
                zh = real2float(i_list[2])
                print(xh, yh, zh, file=tpnew)
                line = tpold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
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
    cbname = conversion_questions.get('CNAME', None)
    if not cbname:
        print("Set central body name:")
        cbname = input("> ")

    print(f'Writing out new CB file: {swiftest_param["CB_IN"]}')
    # Write out new central body file
    try:
        cbnew = open(swiftest_param['CB_IN'], 'w')
        print(cbname, file=cbnew)
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
        swiftest_param['GR'] = True
        swiftest_param.pop('C', None)
    swiftest_param.pop('J2', None)
    swiftest_param.pop('J4', None)
    swiftest_param['IN_FORM'] = "XV"

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
    """
    Converts from a Swift run to a Swiftest run

    Parameters
    ----------
    swift_param : dictionary
       Swift input parameters. 
    plname : string
        Name of massive body input file
    tpname : string
        Name of test particle input file
    cbname : string
        Name of the central body input file
    conversion_questions : dictronary
        Dictionary of additional parameters required to convert between formats

    Returns
    -------
    swiftest_param : A set of input files for a new Swiftest run
    """
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
    """
    Converts from a Swiftest run to a Swifter run

    Parameters
    ----------
    swiftest_param : dictionary
        Swiftest input parameters. 
    J2 : float
        Central body oblateness term. Default spherical.
    J4 : float
        Central body oblateness term. Default spherical.
    Returns
    -------
    swifter_param : A set of input files for a new Swifter run
    """
    swifter_param = swiftest_param
    CBIN = swifter_param.pop("CB_IN", None)
    GMTINY = swifter_param.pop("GMTINY", None)
    DISCARD_OUT = swifter_param.pop("DISCARD_OUT", None)
    MU2KG = swifter_param.pop("MU2KG", 1.0)
    DU2M = swifter_param.pop("DU2M", 1.0)
    TU2S = swifter_param.pop("TU2S", 1.0)
    GR = swifter_param.pop("GR", None)
    if GR is not None:
        if GR:
           swifter_param['C'] =  swiftest.einsteinC * np.longdouble(TU2S) / np.longdouble(DU2M)
    for key in newfeaturelist:
       tmp = swifter_param.pop(key, None)
    swifter_param['J2'] = J2
    swifter_param['J4'] = J4
    if swifter_param['OUT_STAT'] == "REPLACE":
        swifter_param['OUT_STAT'] = "UNKNOWN"
    if swifter_param['OUT_TYPE'] == 'NETCDF_DOUBLE':
       swifter_param['OUT_TYPE'] = 'REAL8'
    elif swifter_param['OUT_TYPE'] == 'NETCDF_FLOAT':
       swifter_param['OUT_TYPE'] = 'REAL4'
    if swifter_param['OUT_FORM'] == 'XVEL':
       swifter_param['OUT_FORM'] = 'XV'
    IN_FORM = swifter_param.pop("IN_FORM", None)
    INTERACTION_LOOPS = swifter_param.pop("INTERACTION_LOOPS", None)
    ENCOUNTER_CHECK = swifter_param.pop("ENCOUNTER_CHECK", None)
    ENCOUNTER_CHECK_PLPL = swifter_param.pop("ENCOUNTER_CHECK_PLPL", None)
    ENCOUNTER_CHECK_PLTP = swifter_param.pop("ENCOUNTER_CHECK_PLTP", None)
    swifter_param['! VERSION'] = "Swifter parameter file converted from Swiftest"

    return swifter_param


def swifter2swift_param(swifter_param, J2=0.0, J4=0.0):
    """
    Converts from a Swifter run to a Swift run

    Parameters
    ----------
    swifter_param : dictionary
        Swifter input parameters. 
    J2 : float
        Central body oblateness term. Default spherical.
    J4 : float
        Central body oblateness term. Default spherical.
    Returns
    -------
    swift_param : A set of input files for a new Swift run
    """
    swift_param = {
        '! VERSION': f"Swift parameter input file converted from Swifter",
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
        'L5': "T",
        'L6': "F",
        'RMIN': -1,
        'RMAX': -1,
        'RMAXU': -1,
        'QMIN': -1,
        'LCLOSE': "F",
        'BINARY_OUTPUTFILE': "bin.dat",
        'STATUS_FLAG_FOR_OPEN_STATEMENTS': "NEW",
    }
    
    swift_param['T0'] = swifter_param['T0']
    swift_param['TSTOP'] = swifter_param['TSTOP']
    swift_param['DT'] = swifter_param['DT']
    # Convert the parameter file values
    swift_param['DTOUT'] = swifter_param['ISTEP_OUT'] * swifter_param['DT']
    swift_param['DTDUMP'] = swifter_param['ISTEP_DUMP'] * swifter_param['DT']
    swift_param['BINARY_OUTPUTFILE'] = swifter_param['BIN_OUT']
    swift_param['STATUS_FLAG_FOR_OPEN_STATEMENTS'] = swifter_param['OUT_STAT']

    if swifter_param['CHK_CLOSE']:
        swift_param['LCLOSE'] = "T"
    else:
        swift_param['LCLOSE'] = "F"

    swift_param['RMIN'] = swifter_param['CHK_RMIN']
    swift_param['RMAX'] = swifter_param['CHK_RMAX']
    swift_param['QMIN'] = swifter_param['CHK_QMIN']
    
    if swift_param['RMIN'] > 0 or swift_param['RMAX'] > 0 or swift_param['QMIN'] > 0:
        swift_param['L2'] = 'T'
        

    return swift_param
