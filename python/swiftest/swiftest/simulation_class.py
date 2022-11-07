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

from swiftest import io
from swiftest import init_cond
from swiftest import tool
from swiftest import constants
from datetime import date
import xarray as xr
import numpy as np
import os
import shutil

class Simulation:
    """
    This is a class that defines the basic Swift/Swifter/Swiftest simulation object
    """
    def __init__(self, codename="Swiftest", param_file="param.in", readbin=False, verbose=True):
        self.ds = xr.Dataset()
        self.param = {
            '! VERSION': f"Swiftest parameter input",
            'T0': "0.0",
            'TSTOP': "0.0",
            'DT': "0.0",
            'IN_FORM': "XV",
            'IN_TYPE': "NETCDF_DOUBLE",
            'NC_IN' : "init_cond.nc",
            'CB_IN' : "cb.in",
            'PL_IN' : "pl.in",
            'TP_IN' : "tp.in",
            'ISTEP_OUT': "1",
            'ISTEP_DUMP': "1",
            'BIN_OUT': "bin.nc",
            'OUT_TYPE': 'NETCDF_DOUBLE',
            'OUT_FORM': "XVEL",
            'OUT_STAT': "REPLACE",
            'CHK_RMAX': "-1.0",
            'CHK_EJECT': "-1.0",
            'CHK_RMIN': "-1.0",
            'CHK_QMIN': "-1.0",
            'CHK_QMIN_COORD': "HELIO",
            'CHK_QMIN_RANGE': "-1.0 -1.0",
            'ENC_OUT': "",
            'MU2KG': constants.MSun,
            'TU2S': constants.JD2S,
            'DU2M': constants.AU2M,
            'EXTRA_FORCE': "NO",
            'DISCARD_OUT': "",
            'PARTICLE_OUT' : "",
            'BIG_DISCARD': "NO",
            'CHK_CLOSE': "YES",
            'RHILL_PRESENT': "YES",
            'FRAGMENTATION': "NO",
            'ROTATION': "NO",
            'TIDES': "NO",
            'ENERGY': "NO",
            'GR': "YES",
            'INTERACTION_LOOPS': "TRIANGULAR",
            'ENCOUNTER_CHECK': "TRIANGULAR"
        }
        self.codename = codename
        self.verbose = verbose
        self.set_unit_system()
        if param_file != "" :
            dir_path = os.path.dirname(os.path.realpath(param_file))
            self.read_param(param_file, codename=codename, verbose=self.verbose)
            if readbin:
                binpath = os.path.join(dir_path,self.param['BIN_OUT'])
                if os.path.exists(binpath):
                    self.param['BIN_OUT'] = binpath
                    self.bin2xr()
                else:
                    print(f"BIN_OUT file {self.param['BIN_OUT']} not found.")
        return


    def set_unit_system(self,MU="Msun",DU="AU",TU="YR",MU2KG=None,DU2M=None,TU2S=None):
        """
        Setter for setting the unit conversion between one of the standard sets.

        The units can be set one of two ways:
        1) The user can supply string values to the arguments MU, DU, and TU to select between common systems
        2) The user can supply float values to the arguments MU2KG, DU2M, and TU2S to manually set the conversion
           factor between the desired unit and the SI unit (kg-m-s).

        The two sets of arguments are mutually exclusive. Any values passed to MU2KG, DU2M, or TU2S will override any
        specified in MU, DU, or TU, respectively. The default system is Msun-AU-YR. MU, DU, and TU are case-insenstive

        Parameters
        ----------
        MU : str, default "Msun"
           The mass unit system to use. Case-insensitive valid options are:
           "Msun"   : Solar mass
           "Mearth" : Earth mass
           "kg"     : kilograms
           "g"      : grams

        DU : str, default "AU"
            The distance unit system to use. Case-insensitive valid options are:
            "AU"     : Astronomical Unit
            "Rearth" : Earth radius
            "m"      : meter
            "cm"     : centimeter

        TU : str, default "YR"
            The time unit system to use. Case-insensitive valid options are:
            "YR"     : Year
            "DAY"    : Julian day
            "d"      : Julian day
            "JD"     : Julian day
            "s"      : second

        MU2KG : float, default `None`
            The conversion factor to multiply by the mass unit that would convert it to kilogram. Setting this overrides MU

        DU2M  : float, default `None`
            The conversion factor to multiply by the distance unit that would convert it to meter. Setting this overrides DU

        TU2S  : float, default `None`
            The conversion factor to multiply by the time unit that would convert it to seconds. Setting this overrides TU

        Returns
        ----------
        Sets the values of MU2KG, DU2M, and TU2S in the param dictionary to the appropriate units. Also computes the gravitational constant, GU
        """

        if MU2KG is not None:
            self.param['MU2KG'] = MU2KG
            self.MU_name = None
        else:
            if MU.upper() == "MSUN":
                self.param['MU2KG'] = constants.MSun
                self.MU_name = "MSun"
            elif MU.upper() == "MEARTH":
                self.param['MU2KG'] = constants.MEarth
                self.MU_name = "MEarth"
            elif MU.upper() == "KG":
                self.param['MU2KG'] = 1.0
                self.MU_name = "kg"
            elif MU.upper() == "G":
                self.param['MU2KG'] = 1000.0
                self.MU_name = "g"
            else:
                print(f"{MU} not a recognized unit system. Using MSun as a default.")
                self.param['MU2KG'] = constants.MSun
                self.MU_name = "MSun"

        if DU2M is not None:
            self.param['DU2M'] = DU2M
            self.DU_name = None
        else:
            if DU.upper() == "AU":
                self.param['DU2M'] = constants.AU2M
                self.DU_name = "AU"
            elif DU.upper() == "REARTH":
                self.param['DU2M'] = constants.REarth
                self.DU_name = "REarth"
            elif DU.upper() == "M":
                self.param['DU2M'] = 1.0
                self.DU_name = "m"
            elif DU.upper() == "CM":
                self.param['DU2M'] = 100.0
                self.DU_name = "cm"
            else:
                print(f"{DU} not a recognized unit system. Using AU as a default.")
                self.param['DU2M'] = constants.AU2M
                self.DU_name = "AU"

        if TU2S is not None:
            self.param['TU2S'] = TU2S
            self.TU_name = None
        else:
            if TU.upper() == "YR":
                self.param['TU2S'] = constants.YR2S
                self.TU_name = "y"
            elif TU.upper() == "DAY" or TU.upper() == "D" or TU.upper() == "JD":
                self.param['TU2S'] = constants.JD2S
                self.TU_name = "Day"
            elif TU.upper() == "S":
                self.param['TU2S'] = 1.0
                self.TU_name = "s"
            else:
                print(f"{TU} not a recognized unit system. Using YR as a default.")
                self.param['TU2S'] = constants.YR2S
                self.TU_name = "y"

        self.VU_name = f"{self.DU_name}/{self.TU_name}"
        self.GC = constants.GC * self.param['TU2S']**2 * self.param['MU2KG'] / self.param['DU2M']**3

        if self.verbose:
            if self.MU_name is None or self.DU_name is None or self.TU_name is None:
                print(f"Custom units set.")
                print(f"MU2KG: {self.param['MU2KG']}")
                print(f"DU2M : {self.param['DU2M']}")
                print(f"TU2S : {self.param['TU2S']}")
            else:
                print(f"Units set to {self.MU_name}-{self.DU_name}-{self.TU_name}")

        return
    
    def add(self, plname, date=date.today().isoformat(), idval=None):
        """
        Adds a solar system body to an existing simulation DataSet.
        
        Parameters
        ----------
           plname : string
                Name of planet to add (e.g. "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"
           date : string
                 Date to use when obtaining the ephemerides in the format YYYY-MM-DD. Defaults to "today"
        Returns
        -------
            self.ds : xarray dataset
        """
        #self.ds = init_cond.solar_system_horizons(plname, idval, self.param, date, self.ds)
        self.addp(*init_cond.solar_system_horizons(plname, idval, self.param, date, self.ds))
        return
    
    
    def addp(self, idvals, namevals, v1, v2, v3, v4, v5, v6, GMpl=None, Rpl=None, rhill=None, Ip1=None, Ip2=None, Ip3=None, rotx=None, roty=None, rotz=None, J2=None,J4=None,t=None):
        """
        Adds a body (test particle or massive body) to the internal DataSet given a set up 6 vectors (orbital elements
        or cartesian state vectors, depending on the value of self.param). Input all angles in degress

        Parameters
        ----------
            idvals : integer 
                Array of body index values.
            v1     : float
                xh for param['IN_FORM'] == "XV"; a for param['IN_FORM'] == "EL"
            v2     : float
                yh for param['IN_FORM'] == "XV"; e for param['IN_FORM'] == "EL"
            v3     : float
                zh for param['IN_FORM'] == "XV"; inc for param['IN_FORM'] == "EL"
            v4     : float
                vhxh for param['IN_FORM'] == "XV"; capom for param['IN_FORM'] == "EL"
            v5     : float
                vhyh for param['IN_FORM'] == "XV"; omega for param['IN_FORM'] == "EL"
            v6     : float
                vhzh for param['IN_FORM'] == "XV"; capm for param['IN_FORM'] == "EL"
            Gmass  : float
                Optional: Array of G*mass values if these are massive bodies
            radius : float
                Optional: Array radius values if these are massive bodies
            rhill  : float
                Optional: Array rhill values if these are massive bodies
            Ip1,y,z : float
                Optional: Principal axes moments of inertia
            rotx,y,z: float
                Optional: Rotation rate vector components
            t      :  float
                Optional: Time at start of simulation
        Returns
        -------
            self.ds : xarray dataset
        """
        if t is None:
            t = self.param['T0']

        dsnew = init_cond.vec2xr(self.param, idvals, namevals, v1, v2, v3, v4, v5, v6, GMpl, Rpl, rhill, Ip1, Ip2, Ip3, rotx, roty, rotz, J2, J4, t)
        if dsnew is not None:
            self.ds = xr.combine_by_coords([self.ds, dsnew])
        self.ds['ntp'] = self.ds['id'].where(np.isnan(self.ds['Gmass'])).count(dim="id")
        self.ds['npl'] = self.ds['id'].where(np.invert(np.isnan(self.ds['Gmass']))).count(dim="id") - 1

        return
    
    
    def read_param(self, param_file, codename="Swiftest", verbose=True):
        """
        Reads in a param.in file and determines whether it is a Swift/Swifter/Swiftest parameter file.
        
        Parameters
        ----------
           param_file : string
                File name of the input parameter file
           codename : string
                 Type of parameter file, either "Swift", "Swifter", or "Swiftest"
        Returns
        -------
            self.ds : xarray dataset
        """
        if codename == "Swiftest":
            self.param = io.read_swiftest_param(param_file, self.param, verbose=verbose)
            self.codename = "Swiftest"
        elif codename == "Swifter":
            self.param = io.read_swifter_param(param_file, verbose=verbose)
            self.codename = "Swifter"
        elif codename == "Swift":
            self.param = io.read_swift_param(param_file, verbose=verbose)
            self.codename = "Swift"
        else:
            print(f'{codename} is not a recognized code name. Valid options are "Swiftest", "Swifter", or "Swift".')
            self.codename = "Unknown"
        return
    
    
    def write_param(self, param_file, param=None):
        """
        Writes to a param.in file and determines whether the output format needs to be converted between Swift/Swifter/Swiftest.
        
        Parameters
        ----------
           param_file : string
                File name of the input parameter file
        Returns
        -------
            self.ds : xarray dataset
        """
        if param is None:
            param = self.param
        # Check to see if the parameter type matches the output type. If not, we need to convert
        codename = param['! VERSION'].split()[0]
        if codename == "Swifter" or codename == "Swiftest":
            if param['IN_TYPE'] == "ASCII":
                param.pop("NC_IN", None)
            else:
                param.pop("CB_IN",None)
                param.pop("PL_IN",None)
                param.pop("TP_IN",None)
            io.write_labeled_param(param, param_file)
        elif codename == "Swift":
            io.write_swift_param(param, param_file)
        else:
            print('Cannot process unknown code type. Call the read_param method with a valid code name. Valid options are "Swiftest", "Swifter", or "Swift".')
        return
    
    
    def convert(self, param_file, newcodename="Swiftest", plname="pl.swiftest.in", tpname="tp.swiftest.in", cbname="cb.swiftest.in", conversion_questions={}):
        """
        Converts simulation input files from one format to another (Swift, Swifter, or Swiftest). 

        Parameters
        ----------
           param_file : string
                File name of the input parameter file
            newcodename : string
                Name of the desired format (Swift/Swifter/Swiftest)
            plname : string
                File name of the massive body input file
            tpname : string
                File name of the test particle input file
            cbname : string
                File name of the central body input file
            conversion_questions : dictronary
                Dictionary of additional parameters required to convert between formats

        Returns
        -------
            oldparam : xarray dataset
                The old parameter configuration.
        """
        oldparam = self.param
        if self.codename == newcodename:
            print(f"This parameter configuration is already in {newcodename} format")
            return oldparam
        if newcodename != "Swift" and newcodename != "Swifter" and newcodename != "Swiftest":
            print(f'{newcodename} is an invalid code type. Valid options are "Swiftest", "Swifter", or "Swift".')
            return oldparam
        goodconversion = True
        if self.codename == "Swifter":
            if newcodename == "Swiftest":
                self.param = io.swifter2swiftest(self.param, plname, tpname, cbname, conversion_questions)
            else:
                goodconversion = False
        elif self.codename == "Swift":
            if newcodename == "Swifter":
                self.param = io.swift2swifter(self.param, plname, tpname, conversion_questions)
            elif newcodename == "Swiftest":
                self.param = io.swift2swiftest(self.param, plname, tpname, cbname, conversion_questions)
            else:
                goodconversion = False
        else:
            goodconversion = False
            
        if goodconversion:
            self.write_param(param_file)
        else:
            print(f"Conversion from {self.codename} to {newcodename} is not supported.")
        return oldparam
    
    
    def bin2xr(self):
        """
        Converts simulation output files from a flat binary file to a xarray dataset. 

        Parameters
        ----------

        Returns
        -------
            self.ds : xarray dataset
        """
        if self.codename == "Swiftest":
            self.ds = io.swiftest2xr(self.param, verbose=self.verbose)
            if self.verbose: print('Swiftest simulation data stored as xarray DataSet .ds')
        elif self.codename == "Swifter":
            self.ds = io.swifter2xr(self.param, verbose=self.verbose)
            if self.verbose: print('Swifter simulation data stored as xarray DataSet .ds')
        elif self.codename == "Swift":
            print("Reading Swift simulation data is not implemented yet")
        else:
            print('Cannot process unknown code type. Call the read_param method with a valid code name. Valid options are "Swiftest", "Swifter", or "Swift".')
        return
    
    
    def follow(self, codestyle="Swifter"):
        """
        An implementation of the Swift tool_follow algorithm. Under development. Currently only for Swift simulations. 

        Parameters
        ----------
            codestyle : string
                Name of the desired format (Swift/Swifter/Swiftest)

        Returns
        -------
            fol : xarray dataset
        """
        if self.ds is None:
            self.bin2xr()
        if codestyle == "Swift":
            try:
                with open('follow.in', 'r') as f:
                    line = f.readline() # Parameter file (ignored because bin2xr already takes care of it
                    line = f.readline() # PL file (ignored)
                    line = f.readline() # TP file (ignored)
                    line = f.readline() # ifol
                    i_list = [i for i in line.split(" ") if i.strip()]
                    ifol = int(i_list[0])
                    line = f.readline()  # nskp
                    i_list = [i for i in line.split(" ") if i.strip()]
                    nskp = int(i_list[0])
            except IOError:
                print('No follow.in file found')
                ifol = None
                nskp = None
            fol = tool.follow_swift(self.ds, ifol=ifol, nskp=nskp)
        else:
            fol = None
        
        if self.verbose: print('follow.out written')
        return fol
    
    
    def save(self, param_file, framenum=-1, codename="Swiftest"):
        """
        Saves an xarray dataset to a set of input files.

        Parameters
        ----------
            param_file : string
                Name of the parameter input file
            framenum : integer (default=-1)
                Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.
            codename : string
                Name of the desired format (Swift/Swifter/Swiftest)

        Returns
        -------
            self.ds : xarray dataset
        """

        if codename == "Swiftest":
            io.swiftest_xr2infile(ds=self.ds, param=self.param, in_type=self.param['IN_TYPE'], framenum=framenum,infile_name=self.param['NC_IN'])
            self.write_param(param_file)
        elif codename == "Swifter":
            if self.codename == "Swiftest":
                swifter_param = io.swiftest2swifter_param(self.param)
            else:
                swifter_param = self.param
            io.swifter_xr2infile(self.ds, swifter_param, framenum)
            self.write_param(param_file, param=swifter_param)
        else:
            print(f'Saving to {codename} not supported')

        return

    def initial_conditions_from_bin(self, framenum=-1, new_param=None, new_param_file="param.new.in", new_initial_conditions_file="bin_in.nc",  restart=False, codename="Swiftest"):
        """
        Generates a set of input files from a old output file.

        Parameters
        ----------
            framenum : integer (default=-1)
                Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.
            new_param : string
                File to copy parameters from. Default is the old parameter file.
            new_param_file : string
                Name of the new parameter file.
            new_initial_conditions_file : string
                Name of the new NetCDF file containing the new initial conditions.
            restart : True or False
                If True, overwrite the old output file. If False, generate a new output file.
            codename : string
                Name of the desired format (Swift/Swifter/Swiftest)

        Returns
        -------
            frame : NetCDF dataset 
        """


        if codename != "Swiftest":
            self.save(new_param_file, framenum, codename)
            return
        if new_param is None:
            new_param = self.param.copy()

        if codename == "Swiftest":
            if restart:
                new_param['T0'] = self.ds.time.values[framenum]
            if self.param['OUT_TYPE'] == 'NETCDF_DOUBLE' or self.param['OUT_TYPE'] == 'REAL8':
                new_param['IN_TYPE'] = 'NETCDF_DOUBLE'
            elif self.param['OUT_TYPE'] == 'NETCDF_FLOAT' or self.param['OUT_TYPE'] == 'REAL4':
                new_param['IN_TYPE'] = 'NETCDF_FLOAT'
            else:
                print(f"{self.param['OUT_TYPE']} is an invalid OUT_TYPE file")
                return
            if self.param['BIN_OUT'] != new_param['BIN_OUT'] and restart:
               print(f"Restart run with new output file. Copying {self.param['BIN_OUT']} to {new_param['BIN_OUT']}")
               shutil.copy2(self.param['BIN_OUT'],new_param['BIN_OUT'])
            new_param['IN_FORM'] = 'XV'
            if restart:
                new_param['OUT_STAT'] = 'APPEND'
            new_param['FIRSTKICK'] = 'T'
            new_param['NC_IN'] = new_initial_conditions_file
            new_param.pop('PL_IN', None)
            new_param.pop('TP_IN', None)
            new_param.pop('CB_IN', None)
            print(f"Extracting data from dataset at time frame number {framenum} and saving it to {new_param['NC_IN']}")
            frame = io.swiftest_xr2infile(self.ds, self.param, infile_name=new_param['NC_IN'],framenum=framenum)
            print(f"Saving parameter configuration file to {new_param_file}")
            self.write_param(new_param_file, param=new_param)

        return frame
