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

from.constants import MSun, AU2M, YR2S, JD2S, GC
from .data import SwiftestDataArray, SwiftestDataset
import numpy as np
import xarray as xr
import sys
import tempfile
import re
import os
import warnings

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
                  "NC_IN",
                  "SEED",
                  "INTERACTION_LOOPS",
                  "ENCOUNTER_CHECK",
                  "ENCOUNTER_CHECK_PLPL",
                  "ENCOUNTER_CHECK_PLTP",
                  "TSTART",
                  "DUMP_CADENCE",
                  "ENCOUNTER_SAVE",
                  "MIN_GMFRAG",
                  "NFRAG_REDUCTION",
                  "COLLISION_MODEL",
                  "COARRAY")

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
              "YORP",
              "COARRAY"]

int_param = ["ISTEP_OUT", "DUMP_CADENCE"]
float_param = ["T0", "TSTART", "TSTOP", "DT", "CHK_RMIN", "CHK_RMAX", "CHK_EJECT", "CHK_QMIN", "DU2M", "MU2KG",
               "TU2S", "MIN_GMFRAG", "GMTINY"]

upper_str_param = ["OUT_TYPE","OUT_FORM","OUT_STAT","IN_TYPE","IN_FORM","ENCOUNTER_SAVE", "CHK_QMIN_COORD"]
lower_str_param = ["NC_IN", "PL_IN", "TP_IN", "CB_IN", "CHK_QMIN_RANGE"]

param_keys = ['! VERSION'] + int_param + float_param + upper_str_param + lower_str_param+ bool_param

# This defines SwiftestDataset variables that are strings, which must be processed due to quirks in how NetCDF-Fortran
# handles strings differently than Python's Xarray.
string_varnames = ["name", "particle_type", "origin_type", "stage", "regime"]
char_varnames = ["space"]
int_varnames = ["id", "ntp", "npl", "nplm", "discard_body_id", "collision_id", "status"]

def _bool2yesno(boolval):
    """
    Converts a boolean into a string of either "YES" or "NO".

    Parameters
    ----------
    boolval : bool
        Input value

    Returns
    -------
    str 
        {"YES","NO"}
    """
    if boolval:
        return "YES"
    else:
        return "NO"
 
 
def _bool2tf(boolval):
    """
    Converts a boolean into a string of either "T" or "F".

    Parameters
    ----------
    boolval : bool
        Input value

    Returns
    -------
    str
        {"T","F"}

    """
    if boolval:
        return "T"
    else:
        return "F"


def _str2bool(input_str):
    """
    Converts a string into an equivalent boolean.

    Parameters
    ----------
    input_str : {"YES", "Y", "T", "TRUE", ".TRUE.", "NO", "N", "F", "FALSE", ".FALSE."}
        Input string. Input is case-insensitive.

    Returns
    -------
    bool 

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


def _real2float(realstr):
    """
    Converts a Fortran-generated ASCII string of a real value into a numpy float type. Handles cases where double precision
    numbers in exponential notation use 'd' or 'D' instead of 'e' or 'E'
    
    Parameters
    ----------
    realstr : str
        Fortran-generated ASCII string of a real value.

    Returns
    -------
    float
        The converted floating point value of the input string
    """
    return float(realstr.replace('d', 'E').replace('D', 'E'))


def read_swiftest_param(param_file_name: os.PathLike, param: dict, verbose: bool=True) -> dict:
    """
    Reads in a Swiftest param.in file and saves it as a dictionary

    Parameters
    ----------
    param_file_name : PathLike
        File name of the input parameter file
    param : dict
        Dictionary to store the entries in the user parameter file
    verbose : bool, default True
        Print out information about the file being read

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
                        alo = _real2float(fields[1])
                        ahi = _real2float(fields[2])
                        param['CHK_QMIN_RANGE'] = f"{alo} {ahi}"

        for uc in upper_str_param:
            if uc in param:
                param[uc] = param[uc].upper()

        for i in int_param:
            if i in param and type(param[i]) != int:
                param[i] = int(float(param[i]))

        for f in float_param:
            if f in param and type(param[f]) is str:
                param[f] = _real2float(param[f])

        for b in bool_param:
            if b in param:
                param[b] = _str2bool(param[b])
    except IOError:
        print(f"{param_file_name} not found.")
    return param


def reorder_dims(ds: SwiftestDataset) -> SwiftestDataset:
    """
    Re-order dimension coordinates so that they are in the same order as the Fortran side
    
    Parameters
    ----------
    ds : SwiftestDataset
    
    Returns
    -------
    ds : SwiftestDataset with the dimensions re-ordered
    """
    idx = ds.indexes
    if "id" in idx:
        dim_order = ["time", "id", "space"]
    elif "name" in idx:
        dim_order = ["time", "name", "space"]
    else:
        dim_order = idx
    idx = {index_name: idx[index_name] for index_name in dim_order}
    ds = ds.reindex(idx)
    return ds


def read_swifter_param(param_file_name: os.PathLike, 
                       verbose: bool=True) -> dict:
    """
    Reads in a Swifter param.in file and saves it as a dictionary

    Parameters
    ----------
    param_file_name : PathLike
        File name of the input parameter file
    verbose : bool, default True
        Print out information about the file being read
    Returns
    -------
    param : dict
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
                        alo = _real2float(fields[1])
                        ahi = _real2float(fields[2])
                        param['CHK_QMIN_RANGE'] = f"{alo} {ahi}"
    
        param['ISTEP_OUT'] = int(param['ISTEP_OUT'])
        param['ISTEP_DUMP'] = int(param['ISTEP_DUMP'])
        param['T0'] = _real2float(param['T0'])
        param['TSTOP'] = _real2float(param['TSTOP'])
        param['DT'] = _real2float(param['DT'])
        param['J2'] = _real2float(param['J2'])
        param['J4'] = _real2float(param['J4'])
        param['CHK_RMIN'] = _real2float(param['CHK_RMIN'])
        param['CHK_RMAX'] = _real2float(param['CHK_RMAX'])
        param['CHK_EJECT'] = _real2float(param['CHK_EJECT'])
        param['CHK_QMIN'] = _real2float(param['CHK_QMIN'])
        for b in bool_param:
            if b in param:
                param[b] = _str2bool(param[b])
        if param['C'] != '-1.0':
            param['C'] = _real2float(param['C'])
        else:
            param.pop('C', None)
    except IOError:
        print(f"{param_file_name} not found.")

    return param


def read_swift_param(param_file_name: os.PathLike, 
                     startfile: os.PathLike="swift.in", 
                     verbose: bool=True) -> dict:
    """
    Reads in a Swift param.in file and saves it as a dictionary

    Parameters
    ----------
    param_file_name : string
        File name of the input parameter file
    startfile : string, default "swift.in"
        File name of the input start file
    verbose : bool, default True
        Print out information about the file being read 
        
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
            param['T0'] = _real2float(line[0])
            param['TSTOP'] = _real2float(line[1])
            param['DT'] = _real2float(line[2])
            line = f.readline().split()
            for i, l in enumerate(line):
                line[i] = l
            param['DTOUT'] = _real2float(line[0])
            param['DTDUMP'] = _real2float(line[1])
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
                param['RMIN'] = _real2float(line[0])
                param['RMAX'] = _real2float(line[1])
                param['RMAXU'] = _real2float(line[2])
                param['QMIN'] = _real2float(line[3])
                param['LCLOSE'] = line[4].upper()
            param['BINARY_OUTPUTFILE'] = f.readline().strip()
            param['STATUS_FLAG_FOR_OPEN_STATEMENTS'] = f.readline().strip().upper()
    except IOError:
       print(f"{param_file_name} not found.")
    
    return param


def write_swift_param(param: dict, param_file_name: os.PathLike) -> None:
    """
    Writes a Swift param.in file.

    Parameters
    ----------
    param : dictionary
        The entries in the user parameter file
    param_file_name : PathLike
        File name of the input parameter file

    Returns
    -------
    None
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


def write_labeled_param(param: dict, param_file_name: os.PathLike) -> None:
    """
    Writes a Swifter/Swiftest param.in file.

    Parameters
    ----------
    param : dictionary
        The entries in the user parameter file
    param_file_name : PathLike
        File name of the input parameter file

    Returns
    -------
    None
    """
    outfile = open(param_file_name, 'w')
    ptmp = param.copy()
    # Print the list of key/value pairs in the preferred order
    for key in param_keys:
        val = ptmp.pop(key, None)
        if val is not None:
            if type(val) is bool:
                print(f"{key:<32} {_bool2yesno(val)}", file=outfile)
            else:
                print(f"{key:<32} {val}", file=outfile)
    # Print the remaining key/value pairs in whatever order
    for key, val in ptmp.items():
        if val != "":
            if type(val) is bool:
                print(f"{key:<32} {_bool2yesno(val)}", file=outfile)
            else:
                print(f"{key:<32} {val}", file=outfile)
    outfile.close()
    return


def _swifter_stream(f, param):
    """
    Reads in a Swifter bin.dat file and returns a single frame of data as a datastream

    Parameters
    ----------
    f : file object
        File object of the Swifter bin.dat file
    param : dict
        Swifter parameters

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
       
        if param['OUT_FORM'] == 'XV':
            rpvec = np.empty((npl,3))
            vpvec = np.empty((npl,3))
            rtvec = np.empty((ntp,3))
            vtvec = np.empty((ntp,3))
        elif param['OUT_FORM'] == 'EL': 
            elempvec = np.empty((npl,6))
            elemtvec = np.empty((ntp,6))
        plid = np.empty(npl, dtype='int')
        tpid = np.empty(ntp, dtype='int')
        if npl > 0:
            GMpl = np.empty(npl)
            Rpl = np.empty(npl)
            for i in range(npl):
                # Read single-line pl frame for
                if param['OUT_FORM'] == 'XV':
                    record = f.read_record('<i4', '<f8', '<f8', '(3,)<f8', '(3,)<f8')
                elif param['OUT_FORM'] == 'EL':
                    record = f.read_record('<i4', '<f8', '<f8', '(6,)<f8')
                plid[i] = record[0]
                GMpl[i] = record[1]
                Rpl[i] = record[2]
                if param['OUT_FORM'] == 'XV':
                    rpvec[i,:] = record[3]
                    vpvec[i,:] = record[4]
                elif param['OUT_FORM'] == 'EL':
                    elempvec[i,:] = record[3]
        if ntp > 0:
            for i in range(ntp):
                if param['OUT_FORM'] == 'XV':
                    record = f.read_record('<i4', '(3,)<f8', '(3,)<f8')
                elif param['OUT_FORM'] == 'EL':
                    record = f.read_record('<i4', '(6,)<f8')
                tpid[i] = record[0]
                if param['OUT_FORM'] == 'XV':
                    rtvec[i,:] = record[1]
                    vtvec[i,:] = record[2]
                elif param['OUT_FORM'] == 'EL':
                    elemtvec[i,:] = record[1]
        
        tlab = ['id']
        if param['OUT_FORM'] == 'XV': 
            tlab.append('rh')
            tlab.append('vh')
        elif param['OUT_FORM'] == 'EL':
            tlab.append('a')
            tlab.append('e')
            tlab.append('inc')
            tlab.append('capom')
            tlab.append('omega')
            tlab.append('capm')
        plab = tlab.copy()
        plab.append('Gmass')
        plab.append('radius')
        
        if param['OUT_FORM'] == 'XV':
            yield t, npl, [plid, rpvec, vpvec, GMpl, Rpl], plab, \
                  ntp, [tpid, rtvec, vtvec], tlab
        elif param['OUT_FORM'] == 'EL':
            yield t, npl, [plid, elempvec, GMpl, Rpl], plab, \
                  ntp, [tpid, elemtvec], tlab


# def swifter2xr(param, verbose=True):
#     """
#     Converts a Swifter binary data file into an SwiftestDataset.

#     Parameters
#     ----------
#     param : dict
#         Swifter parameters
#     verbose : bool, default True
#         Print out information about the file being read
#     Returns
#     -------
#     SwiftestDataset
#     """
#     dims = ['time', 'id','vec']
#     pl = []
#     tp = []
#     with FortranFile(param['BIN_OUT'], 'r') as f:
#         for t, npl, pvec, plab, \
#             ntp, tvec, tlab in _swifter_stream(f, param):
            
#             sys.stdout.write('\r' + f"Reading in time {t[0]:.3e}")
#             sys.stdout.flush()
            
#             pvec_args = dict(zip(plab,pvec))
#             pl.append(vec2xr(param,time=t,**pvec_args))
#             if ntp > 0:
#                 tvec_args = dict(zip(tlab,tvec))
#                 tp.append(vec2xr(param,time=t,**tvec_args))
            
#         plds = xr.concat(pl, dim='time')
#         if len(tp) > 0: 
#             tpds = xr.concat(tp, dim='time')
        
#         if verbose: print('\nCreating Dataset')
#         if len(tp) > 0:
#             ds = xr.combine_by_coords([plds, tpds])
#         if verbose: print(f"Successfully converted {ds.sizes['time']} output frames.")
#     return SwiftestDataset(ds)


def process_netcdf_input(ds: xr.Dataset, param: dict) -> SwiftestDataset:
    """
    Performs several tasks to convert raw NetCDF files output by the Fortran program into a form that
    is used by the Python side. These include:
    - Ensuring all types are correct

    Parameters
    ----------
    ds : Xarray dataset

    Returns
    -------
    ds : SwiftestDataset
    """
    
    if not isinstance(ds, SwiftestDataset):
        ds = SwiftestDataset(ds)

    if param['OUT_TYPE'] == "NETCDF_DOUBLE":
        ds = fix_types(ds,ftype=np.float64)
    elif param['OUT_TYPE'] == "NETCDF_FLOAT":
        ds = fix_types(ds,ftype=np.float32)

    return ds


def swiftest2xr(param: dict, verbose: bool=True, dask: bool=False) -> SwiftestDataset:
    """
    Converts a Swiftest binary data file into an xarray DataSet.

    Parameters
    ----------
    param : dict
        Swiftest parameters
    verbose : bool, default True
        Print out information about the file being read
    dask : bool, default False
        Use Dask to lazily load data (useful for very large datasets)

    Returns
    -------
    SwiftestDataset
    """


    if ((param['OUT_TYPE'] == 'NETCDF_DOUBLE') or (param['OUT_TYPE'] == 'NETCDF_FLOAT')):
        if verbose: print('\nCreating Dataset from NetCDF file')
        if dask:
            ds = xr.open_mfdataset(param['BIN_OUT'], engine='h5netcdf', mask_and_scale=False)
        else:
            ds = xr.open_dataset(param['BIN_OUT'], mask_and_scale=False)
        
        ds = process_netcdf_input(ds, param)
        ds.close()
    else:
        print(f"Error encountered. OUT_TYPE {param['OUT_TYPE']} not recognized.")
        return None
    if verbose: print(f"Successfully converted {ds.sizes['time']} output frames.")

    return ds


def _xstrip_nonstr(da: SwiftestDataArray) -> SwiftestDataArray:
    """
    Cleans up the string values in the DataArray to remove extra white space

    Parameters
    ----------
    da  : SwiftestDataArray
        Input dataset

    Returns
    -------
    SwiftestDataset with the strings cleaned up
    """
    func = lambda x: np.char.strip(x)
    return xr.apply_ufunc(func, da.str.decode(encoding='utf-8'),dask='parallelized')


def _xstrip_str(da: SwiftestDataset) -> SwiftestDataset:
    """
    Cleans up the string values in the DataArray to remove extra white space

    Parameters
    ----------
    da    : SwiftestDataArray
        Input DataArray padded with spaces.

    Returns
    -------
    SwiftestDataset with the strings cleaned up
    """
    func = lambda x: np.char.strip(x)
    return xr.apply_ufunc(func, da,dask='parallelized')


def _string_converter(da: SwiftestDataArray) -> SwiftestDataArray:
    """
    Converts a string to a unicode string 

    Parameters
    ----------
    da    : SwiftestDataArray
        Input DataArray with non-unicode strings

    Returns
    -------
    da : SwiftestDataArray with the strings cleaned up
    """

    da = da.astype('<U32')
    da = _xstrip_str(da)

    return da


def _char_converter(da: SwiftestDataArray) -> SwiftestDataArray:
    """`
    Converts a string to a unicode string

    Parameters
    ----------
    da    : SwiftestDataArray

    Returns
    -------
    da : xarray dataset with the strings cleaned up
    """
    if da.dtype == np.dtype(object):
        da = da.astype('<U1')
    elif type(da.values[0]) != np.str_:
        da = _xstrip_nonstr(da)
    return da


def _clean_string_values(ds: SwiftestDataset) -> SwiftestDataset:
    """
    Cleans up the string values in the DataSet that have artifacts as a result of coming from NetCDF Fortran

    Parameters
    ----------
    ds    : SwiftestDataset
 
    Returns
    -------
    ds : SwiftestDataset with the strings cleaned up
    """

    for n in string_varnames:
        if n in ds:
           ds[n] = _string_converter(ds[n])

    for n in char_varnames:
        if n in ds:
            ds[n] = _char_converter(ds[n])
    return ds


def _unclean_string_values(ds: SwiftestDataset) -> SwiftestDataset:
    """
    Returns strings back to a format readable to NetCDF Fortran

    Parameters
    ----------
    ds    : SwiftestDataset

    Returns
    -------
    ds : SwiftestDataset with the strings cleaned up
    """

    for c in string_varnames:
        if c in ds:
            n = _string_converter(ds[c])
            ds[c] = n.str.ljust(32).str.encode('utf-8')

    for c in char_varnames:
        if c in ds:
            n = _string_converter(ds[c])
            ds[c] = n.str.ljust(1).str.encode('utf-8')
    return ds


def fix_types(ds: SwiftestDataset,
              itype: np.dtype=np.int64,
              ftype: np.dtype=np.float64) -> SwiftestDataset:
    """
    Converts all variables in the dataset to the specified type
    
    Parameters
    ----------
    ds : SwiftestDataset
        Input dataset to convert
    itype : numpy type, default np.int64
        Integer type to convert to
    ftype : numpy type, default np.float64
        Floating point type to convert to 
        
    Returns
    -------
    ds : SwiftestDataset with the types converted
    """
    ds = _clean_string_values(ds)
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


def select_active_from_frame(ds: SwiftestDataset, 
                             param: dict, 
                             framenum: int=-1) -> SwiftestDataset:
    """
    Selects a particular frame from a DataSet and returns only the active particles in that frame

    Parameters
    ----------
    ds : SwiftestDataset
        Dataset containing Swiftest n-body data
    param : dict
        Swiftest input parameters. This method uses the names of the cb, pl, and tp files from the input
    framenum : int, default=-1
        Time frame to extract. If this argument is not passed, the default is to use the last frame in the dataset.

    Returns
    -------
    xarray dataset
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
        if 'rh' in iframe:
            iactive = iframe[count_dim].where((~np.isnan(iframe['Gmass'])) | (~np.isnan(iframe['rh'].isel(space=0))), drop=True)[count_dim]
        else:
            iactive = iframe[count_dim].where(~np.isnan(iframe['Gmass']))
    else:
        if 'a' in iframe:    
            iactive = iframe[count_dim].where((~np.isnan(iframe['Gmass'])) | (~np.isnan(iframe['a'])), drop = True)[count_dim]
        else:
            iactive = iframe[count_dim].where(~np.isnan(iframe['Gmass']))
    if count_dim == "id":
        frame = frame.sel(id=iactive.values)
    elif count_dim == "name":
        frame = frame.sel(name=iactive.values)

    return frame


def swiftest_xr2infile(ds: SwiftestDataset, 
                       param: dict, 
                       in_type: str="NETCDF_DOUBLE", 
                       infile_name: os.PathLike=None,
                       framenum: int=-1,
                       verbose: bool=True) -> None:
    """
    Writes a set of Swiftest input files from a single frame of a Swiftest xarray dataset

    Parameters
    ----------
    ds : SwiftestDataset
        Dataset containing Swiftest n-body data in XV format
    param : dict
        Swiftest input parameters. This method uses the names of the cb, pl, and tp files from the input
    in_type : str, default="NETCDF_DOUBLE"
        Type of input file to write. Options are "NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"
    infile_name : PathLike, default=None
        Name of the input file to write. If not passed, the default is to use the name from the input parameters
    framenum : int (default=-1)
        Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.
    verbose : bool, default True
        Print out information about the file being written.

    Returns
    -------
    None
    """
    param_tmp = param.copy()
    param_tmp['OUT_FORM'] = param['IN_FORM']
    frame = select_active_from_frame(ds, param_tmp, framenum)

    if in_type == "NETCDF_DOUBLE" or in_type == "NETCDF_FLOAT":
        # Convert strings back to byte form and save the NetCDF file
        # Note: xarray will call the character array dimension string32. The Fortran code
        # will rename this after reading

        if infile_name is None:
            infile_name = param['NC_IN']
        frame = _unclean_string_values(frame)
        if verbose:
            print(f"Writing initial conditions to file {infile_name}")
        frame = reorder_dims(frame)
        # This suppresses this warning: RuntimeWarning: numpy.ndarray size changed, may indicate binary incompatibility. Expected 16 from C header, got 96 from PyObject
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            frame.to_netcdf(path=infile_name)
        frame.close()
        return frame

    # All other file types need seperate files for each of the inputs
    cb = frame.isel(name=0)
    pl = frame.where(frame['name'] != cb.name)
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
        Ip1cb = np.double(cb['Ip'].values[0])
        Ip2cb = np.double(cb['Ip'].values[1])
        Ip3cb = np.double(cb['Ip'].values[2])
        rotxcb = np.double(cb['rot'].values[0])
        rotycb = np.double(cb['rot'].values[1])
        rotzcb = np.double(cb['rot'].values[2])
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
                print(pli['rh'].values[0,0], pli['rh'].values[0,1], pli['rh'].values[0,2], file=plfile)
                print(pli['vh'].values[0,0], pli['vh'].values[0,1], pli['vh'].values[0,2], file=plfile)
            elif param['IN_FORM'] == 'EL':
                print(pli['a'].values[0], pli['e'].values[0], pli['inc'].values[0], file=plfile)
                print(pli['capom'].values[0], pli['omega'].values[0], pli['capm'].values[0], file=plfile)
            else:
                print(f"{param['IN_FORM']} is not a valid input format type.")
            if param['ROTATION']:
                print(pli['Ip'].values[0,0], pli['Ip'].values[0,1], pli['Ip'].values[0,2], file=plfile)
                print(pli['rot'].values[0,0], pli['rot'].values[0,1], pli['rot'].values[0,2], file=plfile)
        plfile.close()
        
        # TP file
        tpfile = open(param['TP_IN'], 'w')
        print(tp.id.count().values, file=tpfile)
        for i in tp.id:
            tpi = tp.sel(id=i)
            print(tpi['name'].values[0], file=tpfile)
            if param['IN_FORM'] == 'XV':
                print(tpi['rh'].values[0,0], tpi['rh'].values[0,1], tpi['rh'].values[0,2], file=tpfile)
                print(tpi['vh'].values[0,0], tpi['vh'].values[0,1], tpi['vh'].values[0,2], file=tpfile)
            elif param['IN_FORM'] == 'EL':
                print(tpi['a'].values[0], tpi['e'].values[0], tpi['inc'].values[0], file=tpfile)
                print(tpi['capom'].values[0], tpi['omega'].values[0], tpi['capm'].values[0], file=tpfile)
            else:
                print(f"{param['IN_FORM']} is not a valid input format type.")
        tpfile.close()
    else:
        print(f"{in_type} is an unknown file type")
    return


def swifter_xr2infile(ds: SwiftestDataset, 
                      param: dict, 
                      simdir: os.PathLike=os.getcwd, 
                      framenum: int=-1) -> None:
    """
    Writes a set of Swifter input files from a single frame of a Swiftest xarray dataset

    Parameters
    ----------
    ds : SwiftestDataset
        Dataset containing Swifter n-body data in XV format
    framenum : int
        Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.
    param : dict
        Swifter input parameters. This method uses the names of the pl and tp files from the input

    Returns
    -------
    None
    """

    frame = ds.isel(time=framenum)

    if "name" in frame.dims:
        frame = frame.swap_dims({"name" : "id"})
        frame = frame.reset_coords("name")

    cb = frame.where(frame.id == 0, drop=True)
    pl = frame.where(frame.id > 0, drop=True)
    pl = pl.where(np.invert(np.isnan(pl['Gmass'])), drop=True).drop_vars(['j2rp2', 'j4rp4','c_lm','sign','l','m'],errors="ignore")
    tp = frame.where(np.isnan(frame['Gmass']), drop=True).drop_vars(['Gmass', 'radius', 'j2rp2', 'j4rp4','c_lm','sign','l','m'],errors="ignore")
    
    GMSun = np.double(cb['Gmass'])
    if param['CHK_CLOSE']:
       RSun = np.double(cb['radius'])
    else:
       RSun = param['CHK_RMIN']
    param['J2'] = np.double(cb['j2rp2'])
    param['J4'] = np.double(cb['j4rp4'])
    
    if param['IN_TYPE'] == 'ASCII':
        # Swiftest Central body file
        plfile = open(os.path.join(simdir,param['PL_IN']), 'w')
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
            print(pli['rh'].values[0], pli['rh'].values[1], pli['rh'].values[2], file=plfile)
            print(pli['vh'].values[0], pli['vh'].values[1], pli['vh'].values[2], file=plfile)
        plfile.close()
        
        # TP file
        tpfile = open(os.path.join(simdir,param['TP_IN']), 'w')
        print(tp.id.count().values, file=tpfile)
        for i in tp.id:
            tpi = tp.sel(id=i)
            print(i.values, file=tpfile)
            print(tpi['rh'].values[0], tpi['rh'].values[1], tpi['rh'].values[2], file=tpfile)
            print(tpi['vh'].values[0], tpi['vh'].values[1], tpi['vh'].values[2], file=tpfile)
        tpfile.close()
    else:
        print(f"{param['IN_TYPE']} is an unknown input file type")

    return


def swift2swifter(swift_param: dict, 
                  plname: os.PathLike="", 
                  tpname: os.PathLike="", 
                  conversion_questions: dict={}) -> dict:
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
            GMcb = _real2float(i_list[0])
            if swift_param['L1'] == "T":
                swifter_param['J2'] = _real2float(i_list[1])
                swifter_param['J4'] = _real2float(i_list[2])
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
                GMpl = _real2float(i_list[0])
                if isSyMBA:
                    rhill = _real2float(i_list[1])
                    if swift_param['LCLOSE'] == "T":
                        plrad = _real2float(i_list[2])
                else:
                    if swift_param['LCLOSE'] == "T":
                        plrad = _real2float(i_list[1])
                if swifter_param['RHILL_PRESENT']:
                    print(n + 1, GMpl, rhill, file=plnew)
                else:
                    print(n + 1, GMpl, file=plnew)
                if swifter_param['CHK_CLOSE']:
                    print(plrad, file=plnew)
                line = plold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                rh = _real2float(i_list[0])
                yh = _real2float(i_list[1])
                zh = _real2float(i_list[2])
                print(rh, yh, zh, file=plnew)
                line = plold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                vhx = _real2float(i_list[0])
                vhy = _real2float(i_list[1])
                vhz = _real2float(i_list[2])
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
                rh = _real2float(i_list[0])
                yh = _real2float(i_list[1])
                zh = _real2float(i_list[2])
                print(rh, yh, zh, file=tpnew)
                line = tpold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                vhx = _real2float(i_list[0])
                vhy = _real2float(i_list[1])
                vhz = _real2float(i_list[2])
                print(vhx, vhy, vhz, file=tpnew)
                # Ignore STAT lines
                line = tpold.readline()
                line = tpold.readline()
    except IOError:
        print(f"Error converting TP file")
    swifter_param['! VERSION'] = "Swifter parameter file converted from Swift"
    
    return swifter_param


def swifter2swiftest(swifter_param: dict, 
                     plname: os.PathLike="", 
                     tpname: os.PathLike="", 
                     cbname: os.PathLike="", 
                     conversion_questions: dict={}) -> dict:
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
            GMcb = _real2float(i_list[1])  # Store central body GM for later
            line = plold.readline()  # Ignore the two zero vector lines
            line = plold.readline()
            for n in range(1, npl):  # Loop over planets
                line = plold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                idnum = int(i_list[0])
                GMpl = _real2float(i_list[1])
                if swifter_param['RHILL_PRESENT']:
                   rhill = _real2float(i_list[2])
                   print(idnum, GMpl, rhill, file=plnew)
                else:
                   print(idnum, GMpl, file=plnew)
                if swifter_param['CHK_CLOSE']:
                    line = plold.readline()
                    i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                    plrad = _real2float(i_list[0])
                    print(plrad, file=plnew)
                line = plold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                rh = _real2float(i_list[0])
                yh = _real2float(i_list[1])
                zh = _real2float(i_list[2])
                print(rh, yh, zh, file=plnew)
                line = plold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                vhx = _real2float(i_list[0])
                vhy = _real2float(i_list[1])
                vhz = _real2float(i_list[2])
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
                rh = _real2float(i_list[0])
                yh = _real2float(i_list[1])
                zh = _real2float(i_list[2])
                print(rh, yh, zh, file=tpnew)
                line = tpold.readline()
                i_list = [i for i in re.split('  +|\t',line) if i.strip()]
                vhx = _real2float(i_list[0])
                vhy = _real2float(i_list[1])
                vhz = _real2float(i_list[2])
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
        swiftest_param['MU2KG'] = MSun
        swiftest_param['DU2M'] = AU2M
        swiftest_param['TU2S'] = YR2S
    elif unit_type == 2 or unit_system.upper() == 'MSUN-AU-DAY':
        print("Unit system is MSun-AU-day")
        swiftest_param['MU2KG'] = MSun
        swiftest_param['DU2M'] = AU2M
        swiftest_param['TU2S'] = JD2S
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
    GU = GC / (swiftest_param['DU2M'] ** 3 / (swiftest_param['MU2KG'] * swiftest_param['TU2S'] ** 2))
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
            cbrad = _real2float(cbrad.strip())
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
    if GMTINY != '' and _real2float(GMTINY.strip()) > 0:
        swiftest_param['GMTINY'] = _real2float(GMTINY.strip())
        
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


def swift2swiftest(swift_param: dict, 
                   plname: os.PathLike="", 
                   tpname: os.PathLike="", 
                   cbname: os.PathLike="", 
                   conversion_questions: dict={}) -> dict:
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
    # if GR is not None:
    #     if GR:
    #        swifter_param['C'] =  swiftest.einsteinC * np.longdouble(TU2S) / np.longdouble(DU2M)
    for key in newfeaturelist:
       tmp = swifter_param.pop(key, None)
    if "ISTEP_DUMP" not in swifter_param:
        swifter_param["ISTEP_DUMP"] = swifter_param["ISTEP_OUT"]
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
