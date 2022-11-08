"""
                    self.param['BIN_OUT'] = binpath
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

from swiftest import io
from swiftest import init_cond
from swiftest import tool
from swiftest import constants
from datetime import date
import xarray as xr
import numpy as np
import os
import shutil
from typing import (
    Literal,
    Dict,
    List
)


class Simulation:
    """
    This is a class that defines the basic Swift/Swifter/Swiftest simulation object
    """
    def __init__(self,
                 codename: Literal["Swiftest", "Swifter", "Swift"] = "Swiftest",
                 param_file: os.PathLike | str ="param.in",
                 read_param: bool = False,
                 init_cond_file_type: Literal["NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"] = "NETCDF_DOUBLE",
                 init_cond_file_name: str | os.PathLike | Dict[str, str] | Dict[str, os.PathLike] | None = None,
                 init_cond_format: Literal["EL", "XV"] = "EL",
                 output_file_type: Literal["NETCDF_DOUBLE","NETCDF_FLOAT","REAL4","REAL8","XDR4","XDR8"] = "NETCDF_DOUBLE",
                 output_file_name: os.PathLike | str | None = None,
                 output_format: Literal["XV","XVEL"] = "XVEL",
                 read_old_output_file: bool = False,
                 MU: str = "MSUN",
                 DU: str = "AU",
                 TU: str = "Y",
                 MU2KG: float | None = None,
                 DU2M: float | None = None,
                 TU2S: float | None = None,
                 rmin: float = constants.RSun / constants.AU2M,
                 rmax: float = 10000.0,
                 close_encounter_check: bool = True,
                 general_relativity: bool = True,
                 fragmentation: bool = True,
                 rotation: bool = True,
                 compute_conservation_values: bool = False,
                 extra_force: bool = False,
                 big_discard: bool = False,
                 rhill_present: bool = False,
                 interaction_loops: Literal["TRIANGULAR","FLAT","ADAPTIVE"] = "TRIANGULAR",
                 encounter_check_loops: Literal["TRIANGULAR","SORTSWEEP","ADAPTIVE"] = "TRIANGULAR",
                 verbose: bool = True
                ):
        """

        Parameters
        ----------
        codename : {"Swiftest", "Swifter", "Swift"}, default "Swiftest"
            Name of the n-body integrator that will be used.
        param_file : str, path-like, or file-lke, default "param.in"
            Name of the parameter input file that will be passed to the integrator.
        read_param : bool, default False
            If true, read in a pre-existing parameter input file given by the argument `param_file`. Otherwise, create
            a new parameter file using the arguments passed to Simulation.
            > *Note:* If set to true, the parameters defined in the input file will override any passed into the
            > arguments of Simulation.
        init_cond_file_type : {"NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"}, default "NETCDF_DOUBLE"
            The file type containing initial conditions for the simulation:
            * NETCDF_DOUBLE: A single initial conditions input file in NetCDF file format of type NETCDF_DOUBLE
            * NETCDF_FLOAT: A single initial conditions input file in NetCDF file format of type NETCDF_FLOAT
            * ASCII : Three initial conditions files in ASCII format. The individual files define the central body,
            massive body, and test particle initial conditions.
        init_cond_file_name : str, path-like, or dict, optional
            Name of the input initial condition file or files. Whether to pass a single file name or a dictionary
            depends on the argument passed to `init_cond_file_type`: If `init_cond_file_type={"NETCDF_DOUBLE","NETCDF_FLOAT"}`,
            then this will be a single file name. If `init_cond_file_type="ASCII"` then this must be a dictionary where:
            ```init_cond_file_name = {
                                      "CB" : *path to central body initial conditions file* (Swiftest only),
                                      "PL" : *path to massive body initial conditions file*,
                                      "TP" : *path to test particle initial conditions file*
                                      }```
            If no file name is provided then the following default file names will be used.
            * NETCDF_DOUBLE, NETCDF_FLOAT: `init_cond_file_name = "init_cond.nc"`
            * ASCII: `init_cond_file_name = {"CB" : "cb.in", "PL" : "pl.in", "TP" : "tp.in"}`
        init_cond_format : {"EL", "XV"}, default "EL"
            Indicates whether the input initial conditions are given as orbital elements or cartesian position and
            velocity vectors.
            > *Note:* If `codename` is "Swift" or "Swifter", EL initial conditions are converted to XV.
        output_file_type : {"NETCDF_DOUBLE", "NETCDF_FLOAT","REAL4","REAL8","XDR4","XDR8"}, default "NETCDF_DOUBLE"
            The file type for the outputs of the simulation. Compatible file types depend on the `codename` argument.
            * Swiftest: Only "NETCDF_DOUBLE" or "NETCDF_FLOAT" supported.
            * Swifter: Only "REAL4","REAL8","XDR4" or "XDR8"  supported.
            * Swift: Only "REAL4" supported.
        output_file_name : str or path-like, optional
            Name of output file to generate. If not supplied, then one of the default file names are used, depending on
            the value passed to `output_file_type`. If one of the NetCDF types are used, the default is "bin.nc".
            Otherwise, the default is "bin.dat".
        output_format : {"XV","XVEL"}, default "XVEL"
            Specifies the format for the data saved to the output file. If "XV" then cartesian position and velocity
            vectors for all bodies are stored. If "XVEL" then the orbital elements are also stored.
        read_old_output_file : bool, default False
            If true, read in a pre-existing binary input file given by the argument `output_file_name`.
        MU : str, default "MSUN"
           The mass unit system to use. Case-insensitive valid options are:
           * "Msun"   : Solar mass
           * "Mearth" : Earth mass
           * "kg"     : kilograms
           * "g"      : grams
        DU : str, optional
            The distance unit system to use. Case-insensitive valid options are:
            * "AU"     : Astronomical Unit
            * "Rearth" : Earth radius
            * "m"      : meter
            * "cm"     : centimeter
        TU : str, optional
            The time unit system to use. Case-insensitive valid options are:
            * "YR"     : Year
            * "DAY"    : Julian day
            * "d"      : Julian day
            * "JD"     : Julian day
            * "s"      : second
        MU2KG: float, optional
            The conversion factor to multiply by the mass unit that would convert it to kilogram.
            Setting this overrides MU
        DU2M : float, optional
            The conversion factor to multiply by the distance unit that would convert it to meter.
            Setting this overrides DU
        TU2S : float, optional
            The conversion factor to multiply by the time unit that would convert it to seconds.
            Setting this overrides TU
        rmin : float, default value is the radius of the Sun in the unit system defined by the unit input arguments.
            Minimum distance of the simulation (CHK_QMIN, CHK_RMIN, CHK_QMIN_RANGE[0])
        rmax : float, default value is 10000 AU in the unit system defined by the unit input arguments.
            Maximum distance of the simulation (CHK_RMAX, CHK_QMIN_RANGE[1])
        close_encounter_check : bool, default True
            Check for close encounters between bodies. If set to True, then the radii of massive bodies must be included
            in initial conditions.
        general_relativity : bool, default True
            Include the post-Newtonian correction in acceleration calculations.
        fragmentation : bool, default True
            If set to True, this turns on the Fraggle fragment generation code and `rotation` must also be True.
            This argument only applies to Swiftest-SyMBA simulations. It will be ignored otherwise.
        rotation : bool, default True
            If set to True, this turns on rotation tracking and radius, rotation vector, and moments of inertia values
            must be included in the initial conditions.
            This argument only applies to Swiftest-SyMBA simulations. It will be ignored otherwise.
        compute_conservation_values : bool, default False
            Turns on the computation of energy, angular momentum, and mass conservation and reports the values
            every output step of a running simulation.
        extra_force: bool, default False
            Turns on user-defined force function.
        big_discard: bool, default False
            Includes big bodies when performing a discard (Swifter only)
        rhill_present: bool, default False
            Include the Hill's radius with the input files.
        interaction_loops : {"TRIANGULAR","FLAT","ADAPTIVE"}, default "TRIANGULAR"
            *Swiftest Experimental feature*
            Specifies which algorithm to use for the computation of body-body gravitational forces.
            * "TRIANGULAR" : Upper-triangular double-loops .
            * "FLAT" : Body-body interation pairs are flattened into a 1-D array.
            * "ADAPTIVE" : Periodically times the TRIANGULAR and FLAT methods and determines which one to use based on
              the wall time to complete the loop. *Note:* Using ADAPTIVE means that bit-identical repeatability cannot
              be assured, as the choice of algorithm depends on possible external factors that can influence the wall
              time calculation. The exact floating-point results of the interaction will be different between the two
              algorithm types.
        encounter_check_loops : {"TRIANGULAR","SORTSWEEP","ADAPTIVE"}, default "TRIANGULAR"
            *Swiftest Experimental feature*
            Specifies which algorithm to use for checking whether bodies are in a close encounter state or not.
            * "TRIANGULAR" : Upper-triangular double-loops.
            * "SORTSWEEP" : A Sort-Sweep algorithm is used to reduce the population of potential close encounter bodies.
              This algorithm is still in development, and does not necessarily speed up the encounter checking.
              Use with caution.
            * "ADAPTIVE" : Periodically times the TRIANGULAR and SORTSWEEP methods and determines which one to use based
              on the wall time to complete the loop. *Note:* Using ADAPTIVE means that bit-identical repeatability cannot
              be assured, as the choice of algorithm depends on possible external factors that can influence the wall
              time calculation. The exact floating-point results of the interaction will be different between the two
              algorithm types.
        verbose : bool, default True
            If set to True, then more information is printed by Simulation methods as they are executed. Setting to
            False suppresses most messages other than errors.
        """
        self.ds = xr.Dataset()
        self.param = {
            '! VERSION': f"Swiftest parameter input",
            'T0': 0.0,
            'TSTART' : 0.0,
            'TSTOP': 0.0,
            'DT': 0.0,
            'ISTEP_OUT': 1,
            'ISTEP_DUMP': 1,
            'CHK_QMIN_COORD': "HELIO",
            'INTERACTION_LOOPS': interaction_loops,
            'ENCOUNTER_CHECK': encounter_check_loops
        }
        self.codename = codename
        self.verbose = verbose

        self.set_distance_range(rmin=rmin, rmax=rmax)

        self.set_unit_system(MU=MU, DU=DU, TU=TU,
                             MU2KG=MU2KG, DU2M=DU2M, TU2S=TU2S,
                             recompute_values=True)

        self.set_init_cond_files(init_cond_file_type=init_cond_file_type,
                                 init_cond_file_name=init_cond_file_name,
                                 init_cond_format=init_cond_format)

        self.set_output_files(output_file_type=output_file_type,
                              output_file_name=output_file_name,
                              output_format=output_format)

        self.set_feature(close_encounter_check=close_encounter_check,
                          general_relativity=general_relativity,
                          fragmentation=fragmentation,
                          rotation=rotation,
                          compute_conservation_values=compute_conservation_values,
                          extra_force=extra_force,
                          big_discard=big_discard,
                          rhill_present=rhill_present,
                          verbose = False)

        # If the parameter file is in a different location than the current working directory, we will need
        # to use it to properly open bin files
        self.sim_dir = os.path.dirname(os.path.realpath(param_file))
        if read_param:
            if os.path.exists(param_file):
                self.read_param(param_file, codename=codename, verbose=self.verbose)
            else:
                print(f"{param_file} not found.")

        if read_old_output_file:
            binpath = os.path.join(self.sim_dir,self.param['BIN_OUT'])
            if os.path.exists(binpath):
                self.bin2xr()
            else:
                print(f"BIN_OUT file {binpath} not found.")
        return

    def set_feature(self,
                     close_encounter_check: bool | None = None,
                     general_relativity: bool | None = None,
                     fragmentation: bool | None = None,
                     rotation: bool | None = None,
                     compute_conservation_values: bool | None=None,
                     extra_force: bool | None = None,
                     big_discard: bool | None = None,
                     rhill_present: bool | None = None,
                     tides: bool | None = None,
                     verbose: bool | None = None,
                     ):
        """
        Turns on or off various features of a simulation.

        Parameters
        ----------
        close_encounter_check : bool, optional
            Check for close encounters between bodies. If set to True, then the radii of massive bodies must be included
            in initial conditions.
        general_relativity : bool, optional
            Include the post-Newtonian correction in acceleration calculations.
        fragmentation : bool, optional
            If set to True, this turns on the Fraggle fragment generation code and `rotation` must also be True.
            This argument only applies to Swiftest-SyMBA simulations. It will be ignored otherwise.
        rotation : bool, optional
            If set to True, this turns on rotation tracking and radius, rotation vector, and moments of inertia values
            must be included in the initial conditions.
            This argument only applies to Swiftest-SyMBA simulations. It will be ignored otherwise.
        compute_conservation_values : bool, optional
            Turns on the computation of energy, angular momentum, and mass conservation and reports the values
            every output step of a running simulation.
        extra_force: bool, optional
            Turns on user-defined force function.
        big_discard: bool, optional
            Includes big bodies when performing a discard (Swifter only)
        rhill_present: bool, optional
            Include the Hill's radius with the input files.
        tides: bool, optional
            Turns on tidal model (IN DEVELOPMENT - IGNORED)
        Yarkovsky: bool, optional
            Turns on Yarkovsky model (IN DEVELOPMENT - IGNORED)
        YORP: bool, optional
            Turns on YORP model (IN DEVELOPMENT - IGNORED)
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag

        Returns
        -------
        Sets the appropriate parameters in the self.param dictionary

        """

        update_list = []
        if close_encounter_check is not None:
            self.param["CHK_CLOSE"] = close_encounter_check
            update_list.append("close_encounter_check")

        if general_relativity is not None:
            self.param["GR"] = general_relativity
            update_list.append("general_relativity")

        if fragmentation is not None:
            self.param['FRAGMENTATION'] = fragmentation
            update_list.append("fragmentation")

        if rotation is not None:
            self.param['ROTATION'] = rotation
            update_list.append("rotation")

        if compute_conservation_values is not None:
            self.param["ENERGY"] = compute_conservation_values
            update_list.append("compute_conservation_values")

        if extra_force is not None:
            self.param["EXTRA_FORCE"] = extra_force
            update_list.append("extra_force")

        if big_discard is not None:
            if self.codename != "Swifter":
                self.param["BIG_DISCARD"] = False
            else:
                self.param["BIG_DISCARD"] = big_discard
                update_list.append("big_discard")

        if rhill_present is not None:
            self.param["RHILL_PRESENT"] = rhill_present
            update_list.append("rhill_present")

        if verbose is None:
            verbose = self.verbose
        if verbose:
            if len(update_list) > 0:
                feature_dict = self.get_feature(update_list)
        return


    def get_feature(self,feature: str | List[str] | None = None):
        """

        Returns a subset of the parameter dictionary containing the current value of the feature boolean values.
        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        feature : str | List[str], optional
            A single string or list of strings containing the names of the features to extract. Default is all of:
            ["close_encounter_check", "general_relativity", "fragmentation", "rotation", "compute_conservation_values"]

        Returns
        -------
        feature_dict : dict
           A dictionary containing the requested features.

        """

        valid_features = {"close_encounter_check": "CHK_CLOSE",
                          "general_relativity": "GR",
                          "fragmentation": "FRAGMENTATION",
                          "rotation": "ROTATION",
                          "compute_conservation_values": "ENERGY",
                          "extra_force": "EXTRA_FORCE",
                          "big_discard": "BIG_DISCARD",
                          "rhill_present": "RHILL_PRESENT"
                         }

        if feature is None:
            feature_list = valid_features.keys()
        elif type(feature) is str:
            feature_list = [feature]
        else:
            # Only allow features to be checked if they are valid. Otherwise ignore.
            feature_list = [feat for feat in feature if feat in set(valid_features.keys())]

        # Extract the feature dictionary
        feature_dict = {valid_features[feat]:self.param[valid_features[feat]] for feat in feature_list}

        if self.verbose:
            print("Simulation feature parameters:")
            for key,val in feature_dict.items():
                print(f"{key:<16} {val}")

        return feature_dict

    def set_init_cond_files(self,
        init_cond_file_type: Literal["NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"],
        init_cond_file_name: str | os.PathLike | Dict[str, str] | Dict[str, os.PathLike] | None,
        init_cond_format: Literal["EL", "XV"]):
        """
        Sets the initial condition file parameters in the parameters dictionary.

        Parameters
        ----------
        init_cond_file_type : {"NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"}
            The file type containing initial conditions for the simulation:
            * NETCDF_DOUBLE: A single initial conditions input file in NetCDF file format of type NETCDF_DOUBLE
            * NETCDF_FLOAT: A single initial conditions input file in NetCDF file format of type NETCDF_FLOAT
            * ASCII : Three initial conditions files in ASCII format. The individual files define the central body,
            massive body, and test particle initial conditions.
        init_cond_file_name : str, path-like, or dict, optional
            Name of the input initial condition file or files. Whether to pass a single file name or a dictionary
            depends on the argument passed to `init_cond_file_type`: If `init_cond_file_type={"NETCDF_DOUBLE","NETCDF_FLOAT"}`,
            then this will be a single file name. If `init_cond_file_type="ASCII"` then this must be a dictionary where:
            ```init_cond_file_name = {
                                      "CB" : *path to central body initial conditions file* (Swiftest only),
                                      "PL" : *path to massive body initial conditions file*,
                                      "TP" : *path to test particle initial conditions file*
                                      }```
            If no file name is provided then the following default file names will be used.
            * NETCDF_DOUBLE, NETCDF_FLOAT: `init_cond_file_name = "init_cond.nc"`
            * ASCII: `init_cond_file_name = {"CB" : "cb.in", "PL" : "pl.in", "TP" : "tp.in"}`
        init_cond_format : {"EL", "XV"}
            Indicates whether the input initial conditions are given as orbital elements or cartesian position and
            velocity vectors.
            > *Note:* If `codename` is "Swift" or "Swifter", EL initial conditions are converted to XV.

        Returns
        -------
        Sets the appropriate initial conditions file parameters inside the self.param dictionary.

        """

        def ascii_file_input_error_msg(codename):
            print(f"Error in set_init_cond_files: init_cond_file_name must be a dictionary of the form: ")
            print('{')
            if codename == "Swiftest":
                print('"CB" : *path to central body initial conditions file*,')
            print('"PL" : *path to massive body initial conditions file*,')
            print('"TP" : *path to test particle initial conditions file*')
            print('}')
            return

        if self.codename == "Swiftest":
            init_cond_keys = ["CB", "PL", "TP"]
        else:
            init_cond_keys = ["PL", "TP"]
            if init_cond_file_type != "ASCII":
                print(f"{init_cond_file_type} is not supported by {self.codename}. Using ASCII instead")
                init_cond_file_type="ASCII"
            if init_cond_format != "XV":
                print(f"{init_cond_format} is not supported by {self.codename}. Using XV instead")
                init_cond_format = "XV"

        self.param["IN_TYPE"] = init_cond_file_type
        self.param["IN_FORM"] = init_cond_format

        if init_cond_file_type == "ASCII":
            if init_cond_file_name is None:
                # No file names passed, so we will just use the defaults.
                for key in init_cond_keys:
                    self.param[f"{key}_IN"] = f"{key.lower()}.in"
            elif type(init_cond_file_name) is not dict:
                # Oops, accidentally passed a single string or path-like instead of the expected dictionary for ASCII
                # input type.
                ascii_file_input_error_msg(self.codename)
            elif not all(key in init_cond_file_name for key in init_cond_keys):
                # This is the case where the dictionary doesn't have all the keys we expect. Print an error message.
                ascii_file_input_error_msg(self.codename)
            else:
                # A valid initial conditions file dictionary was passed.
                for key in init_cond_keys:
                    self.param[f"{key}_IN"] = init_cond_file_name[key]
        else:
            if init_cond_file_name is None:
                # No file names passed, so we will just use the default.
                self.param["NC_IN"] = "init_cond.nc"
            elif type(init_cond_file_name) is dict:
                # Oops, accidentally passed a dictionary instead of the expected single string or path-like for NetCDF
                # input type.
                print(f"Only a single input file is used for NetCDF files")
            else:
                self.param["NC_IN"] = init_cond_file_name

        return

    def set_output_files(self,
                         output_file_type: Literal[ "NETCDF_DOUBLE", "NETCDF_FLOAT", "REAL4", "REAL8", "XDR4", "XDR8"],
                         output_file_name: os.PathLike | str | None,
                         output_format: Literal["XV", "XVEL"]
                         ):
        """
        Sets the output file parameters in the parameters dictionary.

        Parameters
        ----------
        output_file_type : {"NETCDF_DOUBLE", "NETCDF_FLOAT","REAL4","REAL8","XDR4","XDR8"}, default "NETCDF_DOUBLE"
            The file type for the outputs of the simulation. Compatible file types depend on the `codename` argument.
            * Swiftest: Only "NETCDF_DOUBLE" or "NETCDF_FLOAT" supported.
            * Swifter: Only "REAL4","REAL8","XDR4" or "XDR8"  supported.
            * Swift: Only "REAL4" supported.
        output_file_name : str or path-like, optional
            Name of output file to generate. If not supplied, then one of the default file names are used, depending on
            the value passed to `output_file_type`. If one of the NetCDF types are used, the default is "bin.nc".
            Otherwise, the default is "bin.dat".
        output_format : {"XV","XVEL"}, default "XVEL"
            Specifies the format for the data saved to the output file. If "XV" then cartesian position and velocity
            vectors for all bodies are stored. If "XVEL" then the orbital elements are also stored.

        Returns
        -------
        Sets the appropriate initial conditions file parameters inside the self.param dictionary.

        """

        if self.codename == "Swiftest":
            if output_file_type not in ["NETCDF_DOUBLE", "NETCDF_FLOAT"]:
                print(f"{output_file_type} is not compatible with Swiftest. Setting to NETCDF_DOUBLE")
                output_file_type = "NETCDF_DOUBLE"
        elif self.codename == "Swifter":
            if output_file_type not in ["REAL4","REAL8","XDR4","XDR8"]:
                print(f"{output_file_type} is not compatible with Swifter. Setting to REAL8")
                output_file_type = "REAL8"
        elif self.codename == "Swift":
            if output_file_type not in ["REAL4"]:
                print(f"{output_file_type} is not compatible with Swift. Setting to REAL4")
                output_file_type = "REAL4"

        self.param['OUT_TYPE'] = output_file_type
        if output_file_name is None:
            if output_file_type in ["NETCDF_DOUBLE", "NETCDF_FLOAT"]:
                self.param['BIN_OUT'] = "bin.nc"
            else:
                self.param['BIN_OUT'] = "bin.dat"
        else:
            self.param['BIN_OUT'] = output_file_name

        if output_format != "XV" and self.codename != "Swiftest":
            print(f"{output_format} is not compatible with {self.codename}. Setting to XV")
            output_format = "XV"
        self.param['OUT_FORM'] = output_format

        return

    def set_unit_system(self,
                        MU: str | None = None,
                        DU: str | None = None,
                        TU: str | None = None,
                        MU2KG: float | None = None,
                        DU2M: float | None = None,
                        TU2S: float | None = None,
                        recompute_values: bool = False):
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
        MU : str, optional
           The mass unit system to use. Case-insensitive valid options are:
           "Msun"   : Solar mass
           "Mearth" : Earth mass
           "kg"     : kilograms
           "g"      : grams
        DU : str, optional
            The distance unit system to use. Case-insensitive valid options are:
            "AU"     : Astronomical Unit
            "Rearth" : Earth radius
            "m"      : meter
            "cm"     : centimeter
        TU : str, optional
            The time unit system to use. Case-insensitive valid options are:
            "YR"     : Year
            "DAY"    : Julian day
            "d"      : Julian day
            "JD"     : Julian day
            "s"      : second
        MU2KG : float, optional
            The conversion factor to multiply by the mass unit that would convert it to kilogram.
            Setting this overrides MU
        DU2M  : float, optional
            The conversion factor to multiply by the distance unit that would convert it to meter.
            Setting this overrides DU
        TU2S  : float, optional
            The conversion factor to multiply by the time unit that would convert it to seconds.
            Setting this overrides TU
        recompute_values : bool, default False
            Recompute all values into the new unit system.
            >*Note:* This is a destructive operation, however if not executed then the values contained in the parameter
            > file and input/output data files computed previously may not be consistent with the new unit conversion
            > factors.

        Returns
        ----------
        Sets the values of MU2KG, DU2M, and TU2S in the param dictionary to the appropriate units. Also computes the
        gravitational constant, self.GU and recomput
        """

        MU2KG_old = None
        DU2M_old = None
        TU2S_old = None

        if MU2KG is not None or MU is not None:
            MU2KG_old = self.param.pop('MU2KG',None)
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

        if DU2M is not None or DU is not None:
            DU2M_old  = self.param.pop('DU2M',None)
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

        if TU2S is not None or TU is not None:
            TU2S_old  = self.param.pop('TU2S',None)
            if TU2S is not None:
                self.param['TU2S'] = TU2S
                self.TU_name = None
            else:
                if TU.upper() == "YR" or TU.upper() == "Y" or TU.upper() == "YEAR" or TU.upper() == "YEARS":
                    self.param['TU2S'] = constants.YR2S
                    self.TU_name = "y"
                elif TU.upper() == "DAY" or TU.upper() == "D" or TU.upper() == "JD" or TU.upper() == "DAYS":
                    self.param['TU2S'] = constants.JD2S
                    self.TU_name = "Day"
                elif TU.upper() == "S" or TU.upper() == "SECONDS" or TU.upper() == "SEC":
                    self.param['TU2S'] = 1.0
                    self.TU_name = "s"
                else:
                    print(f"{TU} not a recognized unit system. Using YR as a default.")
                    self.param['TU2S'] = constants.YR2S
                    self.TU_name = "y"

        self.VU_name = f"{self.DU_name}/{self.TU_name}"
        self.GU = constants.GC * self.param['TU2S']**2 * self.param['MU2KG'] / self.param['DU2M']**3

        if recompute_values:
            self.update_param_units(MU2KG_old, DU2M_old, TU2S_old)

        if self.verbose:
            if self.MU_name is None or self.DU_name is None or self.TU_name is None:
                print(f"Custom units set.")
                print(f"MU2KG: {self.param['MU2KG']}")
                print(f"DU2M : {self.param['DU2M']}")
                print(f"TU2S : {self.param['TU2S']}")
            else:
                print(f"Units set to {self.MU_name}-{self.DU_name}-{self.TU_name}")

        return

    def update_param_units(self,MU2KG_old,DU2M_old,TU2S_old):
        """
        Updates the values of parameters that have units when the units have changed.

        Parameters
        ----------
        MU2KG_old : Old value of the mass unit conversion factor
        DU2M_old :  Old value of the distance unit conversion factor
        TU2S_old :  Old value of the time unit conversion factor

        Returns
        -------
        Updates the set of param dictionary values for the new unit system

        """

        mass_keys = ['GMTINY', 'MIN_GMFRAG']
        distance_keys = ['CHK_QMIN','CHK_RMIN','CHK_RMAX', 'CHK_EJECT']
        time_keys = ['T0','TSTOP','DT']

        if MU2KG_old is not None:
            for k in mass_keys:
                if k in self.param:
                    print(f"param['{k}']: {self.param[k]}")
                    self.param[k] *= self.param['MU2KG'] / MU2KG_old

        if DU2M_old is not None:
            for k in distance_keys:
                if k in self.param:
                    self.param[k] *= self.param['DU2M'] / DU2M_old

            CHK_QMIN_RANGE = self.param.pop('CHK_QMIN_RANGE', None)
            if CHK_QMIN_RANGE is not None:
                CHK_QMIN_RANGE = CHK_QMIN_RANGE.split(" ")
                for i, v in enumerate(CHK_QMIN_RANGE):
                    CHK_QMIN_RANGE[i] = float(CHK_QMIN_RANGE[i]) * self.param['DU2M'] / DU2M_old
                self.param['CHK_QMIN_RANGE'] = f"{CHK_QMIN_RANGE[0]} {CHK_QMIN_RANGE[1]}"

        if TU2S_old is not None:
            for k in time_keys:
                if k in self.param:
                    self.param[k] *= self.param['TU2S'] / TU2S_old

        return


    def set_distance_range(self,rmin=None,rmax=None):
        """
        Sets the minimum and maximum distances of the simulation.

        Parameters
        ----------
        rmin : float
            Minimum distance of the simulation (CHK_QMIN, CHK_RMIN, CHK_QMIN_RANGE[0])
        rmax : float
            Maximum distance of the simulation (CHK_RMAX, CHK_QMIN_RANGE[1])

        Returns
        -------
        Sets the appropriate param dictionary values.

        """
        CHK_QMIN_RANGE = self.param.pop('CHK_QMIN_RANGE',None)
        if CHK_QMIN_RANGE is None:
            CHK_QMIN_RANGE = [-1, -1]
        else:
            CHK_QMIN_RANGE = CHK_QMIN_RANGE.split(" ")
        if rmin is not None:
            self.param['CHK_QMIN'] = rmin
            self.param['CHK_RMIN'] = rmin
            CHK_QMIN_RANGE[0] = rmin
        if rmax is not None:
            self.param['CHK_RMAX'] = rmax
            self.param['CHK_EJECT'] = rmax
            CHK_QMIN_RANGE[1] = rmax

        self.param['CHK_QMIN_RANGE'] =f"{CHK_QMIN_RANGE[0]} {CHK_QMIN_RANGE[1]}"

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

        # Make a temporary copy of the parameter dictionary so we can supply the absolute path of the binary file
        # This is done to handle cases where the method is called from a different working directory than the simulation
        # results
        param_tmp = self.param.copy()
        param_tmp['BIN_OUT'] = os.path.join(self.dir_path,self.param['BIN_OUT'])
        if self.codename == "Swiftest":
            self.ds = io.swiftest2xr(param_tmp, verbose=self.verbose)
            if self.verbose: print('Swiftest simulation data stored as xarray DataSet .ds')
        elif self.codename == "Swifter":
            self.ds = io.swifter2xr(param_tmp, verbose=self.verbose)
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
