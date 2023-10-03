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

from swiftest import io
from swiftest import init_cond
from swiftest import tool
from swiftest import constants
from swiftest._bindings import driver
from swiftest import __file__ as _pyfile
import json
import os
from glob import glob
from pathlib import Path
import datetime
import xarray as xr
import numpy as np
import numpy.typing as npt
import shutil
import warnings
import contextlib
from typing import (
    Literal,
    Dict,
    List,
    Tuple,
    Any
)
from cython import nogil

@contextlib.contextmanager
def _cwd(newdir):
    olddir = os.getcwd()
    os.chdir(newdir)
    try:
        yield
    finally:
        os.chdir(olddir)

class Simulation(object):
    """
    This is a class that defines the basic Swift/Swifter/Swiftest simulation object
    """

    def __init__(self,read_param: bool = False, 
                 read_data: bool = False, 
                 read_collisions: bool | None = None, 
                 read_encounters: bool | None = None,
                 simdir: os.PathLike | str = "simdata", 
                 dask: bool = False,
                 **kwargs: Any):
        """

        Parameters
        ----------
        read_param : bool, default True
            If true, read in a pre-existing parameter input file given by the argument `param_file` if it exists.
            Otherwise, create a new parameter file using the arguments passed to Simulation or defaults

            Parameters for a given Simulation object can be set a number of different ways, including via a parameter input
            file, arguments to Simulation, the general `set_parameter` method, or the specific setters for groups of
            similar parameters (e.g. set_init_cond_files, set_simulation_time, etc.). Each parameter has a default value
            that can be overridden by an argument to Simulation(). Some argument parameters have equivalent values that
            are passed to the Swiftest driver Fortran function via a parameter input file. When declaring a new
            Simulation object, parameters are chosen in the following way, from highest to lowest priority"
            1. Arguments to Simulation()
            2. The parameter input file given by `param_file` under the following conditions:
                - `read_param` is set to True (default behavior).
                - The file given by `param_file` exists. The default file is `param.in` located in the `simdata` directory
                  inside the current working directory, which can be changed by passing `param_file` as an argument.
                - The argument has an equivalent parameter or set of parameters in the parameter input file.
            3. Default values (see below)
        read_data : bool, default False
            If True, read in a pre-existing binary input file given by the argument `output_file_name` if it exists.
            Parameter input file equivalent: None
        read_collisions : bool, default None
            If True, read in a pre-existing collision file `collisions.nc`. If None, then it will take the value of `read_data`. 
        read_encounters : bool, default None
            If True, read in a pre-existing encounter file `encounters.nc`. If None, then it will take the value of `read_data`. 
        simdir : PathLike, default `"simdir"`
            Directory where simulation data will be stored, including the parameter file, initial conditions file, output file,
            dump files, and log files.

        **kwargs : See list of valid parameters and their defaults below

        codename : {"Swiftest", "Swifter", "Swift"}, default "Swiftest"
            Name of the n-body code that will be used.
            Parameter input file equivalent: None
        integrator : {"symba","rmvs","whm","helio"}, default "symba"
            Name of the n-body integrator that will be used when executing a run.
            Parameter input file equivalent: None
        read_param : bool, default False
            Read the parameter file given by `param_file`.
        param_file : str, path-like, or file-lke, default "param.in"
            Name of the parameter input file that will be passed to the integrator.
            Parameter input file equivalent: None
        t0 : float, default 0.0
            The reference time for the start of the simulation. Defaults is 0.0.
            Parameter input file equivalent: `T0`
        tstart : float, default 0.0
            The start time for a restarted simulation. For a new simulation, tstart will be set to t0 automatically.
            Parameter input file equivalent: `TSTART`
        tstop : float, optional
            The stopping time for a simulation. `tstop` must be greater than `tstart`.
            Parameter input file equivalent: `TSTOP`
        dt : float, optional
            The step size of the simulation. `dt` must be less than or equal to `tstop-tstart`.
            Parameter input file equivalent: `DT`
        istep_out : int, optional
            The number of time steps between output saves to file. *Note*: only `istep_out` or `tstep_out` can be set.
            Parameter input file equivalent: `ISTEP_OUT`
        tstep_out : float, optional
            The approximate time between when outputs are written to file. Passing this computes
            `istep_out = floor(tstep_out/dt)`. *Note*: only `istep_out` or `tstep_out` can be set. `tstep_out` must be less than `tstop`
            Parameter input file equivalent: None 
        nstep_out : int, optional
            The total number of times that outputs are written to file. Passing this allows for a geometric progression of output steps:
            `TSTART, f**0 * TSTEP_OUT, f**1 * TSTEP_OUT, f**2 * TSTEP_OUT, ..., f**(nstep_out-1) * TSTEP_OUT`, 
            where `f` is a factor that can stretch (or shrink) the time between outputs. Setting `nstep_out = int((tstart - tstop) / (tstep_out))` is
            equivalent to the standard linear output (i.e. `f==1`) and is the same as not passing anything for this argument. 
            *Note*: Passing `nstep_out` requires passing either `istep_out` or `tstep_out` as well.     
        dump_cadence : int, optional
            The number of output steps (given by `istep_out`) between when the saved data is dumped to a file. Setting it to 0
            is equivalent to only dumping data to file at the end of the simulation. Default value is 10.
            Parameter input file equivalent: `DUMP_CADENCE`
        init_cond_file_type : {"NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"}, default "NETCDF_DOUBLE"
            The file type containing initial conditions for the simulation:
            * NETCDF_DOUBLE: A single initial conditions input file in NetCDF file format of type NETCDF_DOUBLE.
            * NETCDF_FLOAT: A single initial conditions input file in NetCDF file format of type NETCDF_FLOAT.
            * ASCII : Three initial conditions files in ASCII format. The individual files define the central body,
            massive body, and test particle initial conditions.
            Parameter input file equivalent: `IN_TYPE`
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
            Parameter input file equivalent: `NC_IN`, `CB_IN`, `PL_IN`, `TP_IN`
        init_cond_format : {"EL", "XV"}, default "EL"
            Indicates whether the input initial conditions are given as orbital elements or cartesian position and
            velocity vectors.
            > *Note:* If `codename` is "Swift" or "Swifter", EL initial conditions are converted to XV.
            Parameter input file equivalent: `IN_FORM`
        output_file_type : {"NETCDF_DOUBLE", "NETCDF_FLOAT","REAL4","REAL8","XDR4","XDR8"}, default "NETCDF_DOUBLE"
            The file type for the outputs of the simulation. Compatible file types depend on the `codename` argument.
            * Swiftest: Only "NETCDF_DOUBLE" or "NETCDF_FLOAT" supported.
            * Swifter: Only "REAL4","REAL8","XDR4" or "XDR8"  supported.
            * Swift: Only "REAL4" supported.
            Parameter input file equivalent: `OUT_TYPE`
        output_file_name : str or path-like, optional
            Name of output file to generate. If not supplied, then one of the default file names are used, depending on
            the value passed to `output_file_type`. The default is "data.nc".
            Parameter input file equivalent: `BIN_OUT`
        output_format : {"XV","XVEL"}, default "XVEL"
            Specifies the format for the data saved to the output file. If "XV" then cartesian position and velocity
            vectors for all bodies are stored. If "XVEL" then the orbital elements are also stored.
            Parameter input file equivalent: `OUT_FORM`
        MU : str, default "MSUN"
            The mass unit system to use. Case-insensitive valid options are:
            * "Msun"   : Solar mass
            * "Mearth" : Earth mass
            * "kg"     : kilograms
            * "g"      : grams
            Parameter input file equivalent: None
        DU : str, optional
            The distance unit system to use. Case-insensitive valid options are:
            * "AU"     : Astronomical Unit
            * "Rearth" : Earth radius
            * "m"      : meter
            * "cm"     : centimeter
            Parameter input file equivalent: None
        TU : str, optional
            The time unit system to use. Case-insensitive valid options are:
            * "YR"     : Year
            * "DAY"    : Julian day
            * "d"      : Julian day
            * "JD"     : Julian day
            * "s"      : second
            Parameter input file equivalent: None
        MU2KG: float, optional
            The conversion factor to multiply by the mass unit that would convert it to kilogram.
            Setting this overrides MU
            Parameter input file equivalent: `MU2KG`
        DU2M : float, optional
            The conversion factor to multiply by the distance unit that would convert it to meter.
            Setting this overrides DU
            Parameter input file equivalent: `DU2M`
        TU2S : float, optional
            The conversion factor to multiply by the time unit that would convert it to seconds.
            Setting this overrides TU
            Parameter input file equivalent: `TU2S`
        MU_name : str, optional
            The name of the mass unit. When setting one of the standard units via `MU` a name will be
            automatically set for the unit, so this argument will override the automatic name.
            Parameter input file equivalent: None
        DU_name : str, optional
            The name of the distance unit. When setting one of the standard units via `DU` a name will be
            automatically set for the unit, so this argument will override the automatic name.
            Parameter input file equivalent: None
        TU_name : str, optional
            The name of the time unit. When setting one of the standard units via `TU` a name will be
            automatically set for the unit, so this argument will override the automatic name.
            Parameter input file equivalent: None
        rmin : float, default value is the radius of the Sun in the unit system defined by the unit input arguments.
            Minimum distance of the simulation
            Parameter input file equivalent: `CHK_QMIN`, `CHK_RMIN`, `CHK_QMIN_RANGE[0]`
        rmax : float, default value is 10000 AU in the unit system defined by the unit input arguments.
            Maximum distance of the simulation (CHK_RMAX, CHK_QMIN_RANGE[1])
            Parameter input file equivalent: `CHK_RMAX`, `CHK_QMIN_RANGE[1]`
        qmin_coord : str, {"HELIO", "BARY"}, default "HELIO"
            coordinate frame to use for checking the minimum periapsis distance
            Parameter input file equivalent: `QMIN_COORD`
        mtiny : float, optional
            The minimum mass of fully interacting bodies. Bodies below this mass interact with the larger bodies,
            but not each other (SyMBA only). *Note.* Only mtiny or gmtiny is accepted, not both.
            Parameter input file equivalent: None
        gmtiny : float, optional
            The minimum G*mass of fully interacting bodies. Bodies below this mass interact with the larger bodies,
            but not each other (SyMBA only). *Note.* Only mtiny or gmtiny is accepted, not both.
            Parameter input file equivalent: `GMTINY`
        close_encounter_check : bool, default True
            Check for close encounters between bodies. If set to True, then the radii of massive bodies must be included
            in initial conditions.
            Parameter input file equivalent: `CHK_CLOSE`
        encounter_save : {"NONE","TRAJECTORY","CLOSEST", "BOTH"}, default "NONE"
            Indicate if and how encounter data should be saved. If set to "TRAJECTORY", the position and velocity vectors
            of all bodies undergoing close encounters are saved at each intermediate step to the encounter files.
            If set to "CLOSEST", the position  and velocities at the point of closest approach between pairs of bodies are 
            computed and stored to the encounter files. If set to "BOTH", then this stores the values that would be computed
            in "TRAJECTORY" and "CLOSEST". If set to "NONE" no trajectory information is saved.
            *WARNING*: Enabling this feature could lead to very large files.
        general_relativity : bool, default True
            Include the post-Newtonian correction in acceleration calculations.
            Parameter input file equivalent: `GR`
        collision_model: {"MERGE","BOUNCE","FRAGGLE"}, default "MERGE"
            This is used to set the collision/fragmentation model. [TODO: DESCRIBE THESE] 
            This argument only applies to Swiftest-SyMBA simulations. It will be ignored otherwise.
            Parameter input file equivalent: `COLLISION_MODEL`
        minimum_fragment_gmass : float, optional
            If fragmentation is turned on, this sets the mimimum G*mass of a collisional fragment that can be generated if a 
            fragmentation model is enabled. Ignored otherwise.
            *Note.* Only set one of minimum_fragment_gmass or minimum_fragment_mass
            Parameter input file equivalent: None
        minimum_fragment_mass : float, optional
            If fragmentation is turned on, this sets the mimimum mass of a collisional fragment that can be generated. if a 
            fragmentation model is enabled. Ignored otherwise
            *Note.* Only set one of minimum_fragment_gmass or minimum_fragment_mass
            Parameter input file equivalent: `MIN_GMFRAG`
        nfrag_reduction : float, optional
            If fragmentation is turne don, this is a reduction factor used to limit the number of fragments generated in a collision.
            For instance, if the SFD of the collision would generated 300 fragments above the `minimum_fragment_mass`, then a value
            of `nfrag_reduction = 30.0` would reduce it to 10.  
            *Note.* Currently only used by the Fraggle collision model.
        rotation : bool, default False
            If set to True, this turns on rotation tracking and radius, rotation vector, and moments of inertia values
            must be included in the initial conditions.
            This argument only applies to Swiftest-SyMBA simulations. It will be ignored otherwise.
            Parameter input file equivalent: `ROTATION`
        compute_conservation_values : bool, default False
            Turns on the computation of energy, angular momentum, and mass conservation and reports the values
            every output step of a running simulation.
            Parameter input file equivalent: `ENERGY`
        extra_force: bool, default False
            Turns on user-defined force function.
            Parameter input file equivalent: `EXTRA_FORCE`
        big_discard: bool, default False
            Includes big bodies when performing a discard (Swifter only)
            Parameter input file equivalent: `BIG_DISCARD`
        rhill_present: bool, default False
            Include the Hill's radius with the input files .
            Parameter input file equivalent: `RHILL_PRESENT`
        restart : bool, default False
            If true, will restart an old run. The file given by `output_file_name` must exist before the run can
            execute. If false, will start a new run. If the file given by `output_file_name` exists, it will be replaced
            when the run is executed.
            Parameter input file equivalent: `OUT_STAT`
        interaction_loops : {"TRIANGULAR","FLAT"}, default "TRIANGULAR"
            > *Swiftest Experimental feature*
            Specifies which algorithm to use for the computation of body-body gravitational forces.
            * "TRIANGULAR" : Upper-triangular double-loops .
            * "FLAT" : Body-body interation pairs are flattened into a 1-D array.
            Parameter input file equivalent: `INTERACTION_LOOPS`
        encounter_check_loops : {"TRIANGULAR","SORTSWEEP"}, default "TRIANGULAR"
            > *Swiftest Experimental feature*
            Specifies which algorithm to use for checking whether bodies are in a close encounter state or not.
            * "TRIANGULAR" : Upper-triangular double-loops.
            * "SORTSWEEP" : A Sort-Sweep algorithm is used to reduce the population of potential close encounter bodies.
              This algorithm is still in development, and does not necessarily speed up the encounter checking.
              Use with caution.
            Parameter input file equivalent: `ENCOUNTER_CHECK`
        dask : bool, default False
            Use Dask to lazily load data (useful for very large datasets)
        coarray : bool, default False
            If true, will employ Coarrays on test particle structures to run in single program/multiple data parallel mode. 
            *Note" In order to use this capability, Swiftest must be compiled for Coarray support. Only certain integrators
            can use Coarrays: RMVS, WHM, Helio are all compatible, but SyMBA is not, due to the way tp-pl close encounters 
            are handeled.
        verbose : bool, default True
            If set to True, then more information is printed by Simulation methods as they are executed. Setting to
            False suppresses most messages other than errors.
            Parameter input file equivalent: None
        """

        self._getter_column_width = 32
        self.verbose = kwargs.pop("verbose",True)

        self.param = {}
        self.data = xr.Dataset()
        self.init_cond = xr.Dataset()
        self.encounters = xr.Dataset()
        self.collisions = xr.Dataset()

        # Set the location of the parameter input file, choosing the default if it isn't specified.
        self.simdir = Path.cwd() / Path(simdir)
        param_file = Path(kwargs.pop("param_file", "param.in"))

        # Parameters are set in reverse priority order. First the defaults, then values from a pre-existing input file,
        # then using the arguments passed via **kwargs.
        #--------------------------
        # Lowest Priority: Defaults:
        #--------------------------
        # Quietly set all parameters to their defaults.
        self.set_parameter(verbose=False,param_file=param_file)

        #-----------------------------------------------------------------
        # Higher Priority: Values from a file (if requested and it exists)
        #-----------------------------------------------------------------

        # If the user asks to read in an old parameter file or output file, override any default parameters with values from the file
        # If the file doesn't exist, flag it for now so we know to create it
        param_file_found = False
        
        if read_collisions is None:
            self.read_collisions = read_data
        else:
            self.read_collisions = read_collisions
        
        if read_encounters is None:
            self.read_encounters = read_data
        else:
            self.read_encounters = read_encounters
            
        if read_param or read_data:
            if self.read_param(read_init_cond = True, dask=dask, **kwargs):
                # We will add the parameter file to the kwarg list. This will keep the set_parameter method from
                # overriding everything with defaults when there are no arguments passed to Simulation()
                kwargs['param_file'] = self.param_file
                param_file_found = True
            else:
                param_file_found = False

        # -----------------------------------------------------------------
        # Highest Priority: Values from arguments passed to Simulation()
        # -----------------------------------------------------------------
        if len(kwargs) > 0 and "param_file" not in kwargs or len(kwargs) > 1 and "param_file" in kwargs:
            self.set_parameter(verbose=False, **kwargs)

        # Let the user know that there was a problem reading an old parameter file and we're going to create a new one
        if read_param and not param_file_found:
            warnings.warn(f"{self.param_file} not found. Creating a new file using default values for parameters not passed to Simulation().",stacklevel=2)
            self.write_param()

        # Read in an old simulation file if requested
        if read_data:
            binpath = os.path.join(self.simdir, self.param['BIN_OUT'])
            if os.path.exists(binpath):
                self.read_output_file(dask=dask)
            else:
                raise FileNotFoundError(f"BIN_OUT file {binpath} not found.")
        return

    def _run_swiftest_driver(self):
        """
        Internal callable function that executes the swiftest_driver run
        """

        with _cwd(self.simdir):
            driver(self.integrator,str(self.param_file), "progress")

        return

    def run(self,dask: bool = False, **kwargs):
        """
        Runs a Swiftest integration. Uses the parameters set by the `param` dictionary unless overridden by keyword
        arguments. Accepts any keyword arguments that can be passed to `set_parameter`.

        Parameters
        ----------
        **kwargs : Any valid keyword arguments accepted by `set_parameter`

        Returns
        -------
        None

        """

        if len(kwargs) > 0:
            self.set_parameter(**kwargs)

        if self.codename != "Swiftest":
            warnings.warn(f"Running an integration is not yet supported for {self.codename}",stacklevel=2)
            return

        # Save initial conditions
        if not self.restart:
            self.clean()
            
        # Write out the current parameter set before executing run
        self.write_param(verbose=False,**kwargs)

        print(f"Running a {self.codename} {self.integrator} run from tstart={self.param['TSTART']} {self.TU_name} to tstop={self.param['TSTOP']} {self.TU_name}")

        self._run_swiftest_driver()

        # Read in new data
        self.read_encounters = True
        self.read_collisions = True
        self.read_output_file(dask=dask)

        return

    def _get_valid_arg_list(self, arg_list: str | List[str] | None = None, valid_var: Dict | None = None):
        """
        Internal function for getters that extracts subset of arguments that is contained in the dictionary of valid
        argument/parameter variable pairs.

        Parameters
        ----------
        arg_list : str | List[str], optional
            A single string or list of strings containing the Simulation argument. If none are supplied,
            then it will create the arg_list out of all keys in the valid_var dictionary.
        valid_var : valid_var: Dict
            A dictionary where the key is the argument name and the value is the equivalent key in the Simulation
            parameter dictionary (i.e. the left-hand column of a param.in file)

        Returns
        -------
        valid_arg : [str]
            The list of valid arguments that match the keys in valid_var
        param : dict
            A dictionary that is the subset of the Simulation parameter dictionary corresponding to the arguments listed
            in arg_list.

        """

        if arg_list is None:
            valid_arg = None
        else:
            valid_arg = arg_list.copy()

        if valid_arg is None:
            valid_arg = list(valid_var.keys())
        elif type(arg_list) is str:
            valid_arg = [arg_list]
        else:
            # Only allow arg_lists to be checked if they are valid. Otherwise ignore.
            valid_arg = [k for k in arg_list if k in list(valid_var.keys())]

        # Extract the arg_list dictionary
        param = {valid_var[arg]: self.param[valid_var[arg]] for arg in valid_arg if valid_var[arg] in self.param}

        return valid_arg, param

    def set_simulation_time(self,
                            t0: float | None = None,
                            tstart: float | None = None,
                            tstop: float | None = None,
                            dt: float | None = None,
                            istep_out: int | None = None,
                            tstep_out: float | None = None,
                            nstep_out: int | None = None,
                            dump_cadence: int | None = None,
                            verbose: bool | None = None,
                            **kwargs: Any
                            ):
        """

        Parameters
        ----------
        t0 : float, optional
            The reference time for the start of the simulation. Defaults is 0.0
        tstart : float, optional
            The start time for a restarted simulation. For a new simulation, tstart will be set to t0 automatically.
        tstop : float, optional
            The stopping time for a simulation. `tstop` must be greater than `tstart`.
        dt : float, optional
            The step size of the simulation. `dt` must be less than or equal to `tstop-dstart`.
        istep_out : int, optional
            The number of time steps between output saves to file. *Note*: only `istep_out` or `tstep_out` can be set.
            Parameter input file equivalent: `ISTEP_OUT`
        tstep_out : float, optional
            The approximate time between when outputs are written to file. Passing this computes
            `istep_out = floor(tstep_out/dt)`. *Note*: only `istep_out` or `tstep_out` can be set.
            Parameter input file equivalent: None 
        nstep_out : int, optional
            The total number of times that outputs are written to file. Passing this allows for a geometric progression of output steps:
            `TSTART, f**0 * TSTEP_OUT, f**1 * TSTEP_OUT, f**2 * TSTEP_OUT, ..., f**(nstep_out-1) * TSTEP_OUT`, 
            where `f` is a factor that can stretch (or shrink) the time between outputs. Setting `nstep_out = int((tstart - tstop) / (tstep_out))` is
            equivalent to the standard linear output (i.e. `f==1`) and is the same as not passing anything for this argument. 
            *Note*: Passing `nstep_out` requires passing either `istep_out` or `tstep_out` as well.     
        dump_cadence : int, optional
            The number of output steps (given by `istep_out`) between when the saved data is dumped to a file. Setting it to 0
            is equivalent to only dumping data to file at the end of the simulation. Default value is 10.
            Parameter input file equivalent: `DUMP_CADENCE`
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        time_dict : dict
           A dictionary containing the requested parameters

        """
        if t0 is None and tstart is None and tstop is None and dt is None and istep_out is None and \
                tstep_out is None and dump_cadence is None:
            return {}

        update_list = []

        if t0 is None:
            t0 = self.param.pop("T0", None)
            if t0 is None:
                t0 = 0.0
        else:
            update_list.append("t0")

        if tstart is None:
            tstart = self.param.pop("TSTART", None)
            if tstart is None:
                tstart = t0
        else:
            update_list.append("tstart")

        self.param['T0'] = t0
        self.param['TSTART'] = tstart

        if tstop is None:
            tstop = self.param.pop("TSTOP", None)
        else:
            update_list.append("tstop")

        if tstop is not None:
            if tstop <= tstart:
                warnings.warn("tstop must be greater than tstart.",stacklevel=2)
                return {}

        if tstop is not None:
            self.param['TSTOP'] = tstop

        if dt is None:
            dt = self.param.pop("DT", None)
        else:
            update_list.append("dt")

        if dt is not None and tstop is not None:
            if dt > (tstop - tstart):
                msg = "dt must be smaller than tstop-tstart"
                msg +=f"\nSetting dt = {tstop - tstart} instead of {dt}"
                warnings.warn(msg,stacklevel=2)
                dt = tstop - tstart

        if dt is not None:
            self.param['DT'] = dt

        if istep_out is None and tstep_out is None:
            istep_out = self.param.pop("ISTEP_OUT", None)
        elif istep_out is not None and tstep_out is not None:
            warnings.warn("istep_out and tstep_out cannot both be set",stacklevel=2)
            return {}
        else:
            update_list.append("istep_out")
            
        if tstep_out is not None and tstep_out > tstop:
            warnings.warn("tstep_out must be less than tstop. Setting tstep_out=tstop",stacklevel=2)
            tstep_out = tstop

        if tstep_out is not None and dt is not None:
            istep_out = int(tstep_out / dt)

        if istep_out is not None:
            self.param['ISTEP_OUT'] = int(istep_out)
            
            
        if nstep_out is not None:
            if istep_out is None:
                warnings.warn("nstep_out requires either istep_out or tstep_out to also be set", stacklevel=2)
            else:
                self.param['NSTEP_OUT'] = int(nstep_out)

        if dump_cadence is None:
            dump_cadence = self.param.pop("DUMP_CADENCE", 1)
        else:
            update_list.append("dump_cadence")
        self.param['DUMP_CADENCE'] = dump_cadence
        
        time_dict = self.get_simulation_time(update_list, verbose=verbose)

        return time_dict

    def get_simulation_time(self, arg_list: str | List[str] | None = None, verbose: bool | None = None, **kwargs):
        """

        Returns a subset of the parameter dictionary containing the current simulation time parameters.
        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list : str | List[str], optional
            A single string or list of strings containing the names of the simulation time parameters to extract.
            Default is all of:
            ["t0", "tstart", "tstop", "dt", "istep_out", "tstep_out", "dump_cadence"]
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.


        Returns
        -------
        time_dict : dict
           A dictionary containing the requested parameters

        """

        valid_var = {"t0": "T0",
                     "tstart": "TSTART",
                     "tstop": "TSTOP",
                     "dt": "DT",
                     "istep_out": "ISTEP_OUT",
                     "nstep_out": "NSTEP_OUT",
                     "dump_cadence": "DUMP_CADENCE",
                     }

        units = {"t0": self.TU_name,
                 "tstart": self.TU_name,
                 "tstop": self.TU_name,
                 "dt": self.TU_name,
                 "tstep_out": self.TU_name,
                 "istep_out": "",
                 "nstep_out": "",
                 "dump_cadence": ""}

        tstep_out = None
        if arg_list is None or "tstep_out" in arg_list or "istep_out" in arg_list:
            if "ISTEP_OUT" in self.param and "DT" in self.param:
                istep_out = self.param['ISTEP_OUT']
                dt = self.param['DT']
                tstep_out = istep_out * dt

        valid_arg, time_dict = self._get_valid_arg_list(arg_list, valid_var)

        if verbose is None:
            verbose = self.verbose

        if verbose:
            for arg in valid_arg:
                key = valid_var[arg]
                if key in time_dict:
                    print(f"{arg:<{self._getter_column_width}} {time_dict[key]} {units[arg]}")
                else:
                    print(f"{arg:<{self._getter_column_width}} NOT SET")
            if tstep_out is not None:
                print(f"{'tstep_out':<{self._getter_column_width}} {tstep_out} {units['tstep_out']}")

        return time_dict

    def set_parameter(self, verbose: bool = True, **kwargs):
        """
        Setter for all possible parameters. This will call each of the specialized setters using keyword arguments.
        If no arguments are passed, then default values will be used.
        Parameters
        ----------
        **kwargs : Any argument listed listed in the Simulation class definition.

        Returns
        -------
        param : A dictionary of all Simulation parameters that changed

        """

        default_arguments = {
            "codename" : "Swiftest",
            "integrator": "symba",
            "t0": 0.0,
            "tstart": 0.0,
            "tstop": None,
            "dt": None,
            "istep_out": 1,
            "tstep_out": None,
            "nstep_out": None,
            "dump_cadence": 10,
            "init_cond_file_type": "NETCDF_DOUBLE",
            "init_cond_file_name": None,
            "init_cond_format": "EL",
            "read_data": False,
            "output_file_type": "NETCDF_DOUBLE",
            "output_file_name": None,
            "output_format": "XVEL",
            "MU": "MSUN",
            "DU": "AU",
            "TU": "Y",
            "MU2KG": None,
            "DU2M": None,
            "TU2S": None,
            "MU_name": None,
            "DU_name": None,
            "TU_name": None,
            "rmin": constants.RSun / constants.AU2M,
            "rmax": 10000.0,
            "qmin_coord": "HELIO",
            "gmtiny": 0.0,
            "mtiny": None,
            "nfrag_reduction": 30.0,
            "close_encounter_check": True,
            "general_relativity": True,
            "collision_model": "FRAGGLE",
            "minimum_fragment_mass": None,
            "minimum_fragment_gmass": 0.0,
            "rotation": True,
            "compute_conservation_values": False,
            "extra_force": False,
            "big_discard": False,
            "rhill_present": False,
            "interaction_loops": "TRIANGULAR",
            "encounter_check_loops": "TRIANGULAR",
            "ephemeris_date": "MBCL",
            "restart": False,
            "encounter_save" : "NONE",
            "coarray" : False,
            "simdir" : self.simdir,
        }
        param_file = kwargs.pop("param_file",None)

        if param_file is not None:
            self.param_file = param_file

        # If no arguments (other than, possibly, verbose) are requested, use defaults
        if len(kwargs) == 0:
            kwargs = default_arguments

        unrecognized = [k for k,v in kwargs.items() if k not in default_arguments]
        if len(unrecognized) > 0:
            for k in unrecognized:
                warnings.warn(f'Unrecognized argument "{k}"',stacklevel=2)

        # Add the verbose flag to the kwargs for passing down to the individual setters
        kwargs["verbose"] = verbose

        # Setters returning parameter dictionary values
        param_dict = {}
        param_dict.update(self.set_integrator(**kwargs))
        param_dict.update(self.set_unit_system(**kwargs))
        param_dict.update(self.set_simulation_time(**kwargs))
        param_dict.update(self.set_init_cond_files(**kwargs))
        param_dict.update(self.set_output_files(**kwargs))
        param_dict.update(self.set_distance_range(**kwargs))
        param_dict.update(self.set_feature(**kwargs))

        # Non-returning setters
        self.set_ephemeris_date(**kwargs)

        return param_dict

    def get_parameter(self, **kwargs):
        """
        Setter for all possible parameters. Calls each of the specialized setters using keyword arguments
        Parameters
        ----------
        **kwargs : Any of the arguments defined in Simulation. If none provided, it returns all arguments.

        Returns
        -------
        param : A dictionary of all Simulation parameters requested

        """

        # Getters returning parameter dictionary values
        param_dict = {}
        param_dict.update(self.get_integrator(**kwargs))
        param_dict.update(self.get_simulation_time(**kwargs))
        param_dict.update(self.get_init_cond_files(**kwargs))
        param_dict.update(self.get_output_files(**kwargs))
        param_dict.update(self.get_distance_range(**kwargs))
        param_dict.update(self.get_unit_system(**kwargs))
        param_dict.update(self.get_feature(**kwargs))

        self.get_ephemeris_date(**kwargs)

        return param_dict

    def set_integrator(self,
                       codename: None | Literal["Swiftest", "Swifter", "Swift"]  = "Swiftest",
                       integrator: Literal["symba","rmvs","whm","helio"] | None = None,
                       mtiny: float | None = None,
                       gmtiny: float | None = None,
                       verbose: bool | None = None,
                       **kwargs: Any
                       ):
        """

        Parameters
        ----------
        codename : {"swiftest", "swifter", "swift"}, optional
        integrator : {"symba","rmvs","whm","helio"}, optional
            Name of the n-body integrator that will be used when executing a run.
        mtiny : float, optional
            The minimum mass of fully interacting bodies. Bodies below this mass interact with the larger bodies,
            but not each other (SyMBA only). *Note.* Only mtiny or gmtiny is accepted, not both.
        gmtiny : float, optional
            The minimum G*mass of fully interacting bodies. Bodies below this mass interact with the larger bodies,
            but not each other (SyMBA only). *Note.* Only mtiny or gmtiny is accepted, not both.
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        integrator_dict: dict
            A dictionary containing the subset of the parameter dictonary that was updated by this setter

        """
        # TODO: Improve how it finds the executable binary

        update_list = []
        
        if codename is not None:
            valid_codename = ["Swiftest", "Swifter", "Swift"]
            if codename.title() not in valid_codename:
                warnings.warn(f"{codename} is not a valid codename. Valid options are ",",".join(valid_codename),stacklevel=2)
                try:
                    self.codename
                except:
                    self.codename = valid_codename[0]
            else:
                self.codename = codename.title()

            self.param['! VERSION'] = f"{self.codename} input file"
            update_list.append("codename")
        
        if self.codename == "Swifter":
            J2 = self.param.pop("J2",0.0)
            J4 = self.param.pop("J4",0.0)
            self.param = io.swiftest2swifter_param(self.param, J2, J4) 

        if integrator is not None:
            valid_integrator = ["symba","rmvs","whm","helio"]
            if integrator.lower() not in valid_integrator:
                warnings.warn(f"{integrator} is not a valid integrator. Valid options are ",",".join(valid_integrator),stacklevel=2)
                try:
                    self.integrator
                except:
                    self.integrator = valid_integrator[0]
            else:
                self.integrator = integrator.lower()
            update_list.append("integrator")

        if mtiny is not None or gmtiny is not None:
            if self.integrator != "symba":
                warnings.warn("mtiny and gmtiny are only used by SyMBA.",stacklevel=2)
            if mtiny is not None and gmtiny is not None:
                warnings.warn("Only set mtiny or gmtiny, not both.",stacklevel=2)
            elif gmtiny is not None:
                self.param['GMTINY'] = gmtiny
                update_list.append("gmtiny")
            elif mtiny is not None:
                self.param['GMTINY'] = self.GU * mtiny
                update_list.append("gmtiny")

        integrator_dict = self.get_integrator(update_list, verbose)

        return integrator_dict

    def get_integrator(self,arg_list: str | List[str] | None = None, verbose: bool | None = None, **kwargs: Any):
        """

        Returns a subset of the parameter dictionary containing the current values of the distance range parameters.
        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list: str | List[str], optional
            A single string or list of strings containing the names of the features to extract. Default is all of:
            ["integrator"]
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        integrator_dict : dict
            The subset of the dictionary containing the code name if codename is selected
        """

        valid_var = {"gmtiny"  : "GMTINY"}

        valid_instance_vars = {"codename": self.codename,
                               "integrator": self.integrator,
                               "param_file": str(self.param_file),
                              }

        try:
            self.integrator
        except:
            warnings.warn(f"integrator is not set",stacklevel=2)
            return {}

        try:
            self.codename
        except:
            warnings.warn(f"codename is not set",stacklevel=2)
            return {}

        if verbose is None:
            verbose = self.verbose

        if not bool(kwargs) and arg_list is None:
            arg_list = list(valid_instance_vars.keys())
            arg_list.append(*[a for a in valid_var.keys() if a not in valid_instance_vars])

        valid_arg, integrator_dict = self._get_valid_arg_list(arg_list, valid_var)

        if verbose:
            for arg in arg_list:
                if arg in valid_instance_vars:
                    print(f"{arg:<{self._getter_column_width}} {valid_instance_vars[arg]}")
            for arg in valid_arg:
                key = valid_var[arg]
                if key in integrator_dict:
                    if arg == "gmtiny":
                        if self.integrator == "symba":
                           print(f"{arg:<{self._getter_column_width}} {integrator_dict[key]} {self.DU_name}^3 / {self.TU_name}^2 ")
                           print(f"{'mtiny':<{self._getter_column_width}} {integrator_dict[key] / self.GU} {self.MU_name}")
                    else:
                        print(f"{arg:<{self._getter_column_width}} {integrator_dict[key]}")
                else:
                    print(f"{arg:<{self._getter_column_width}} NOT SET")


        return integrator_dict

    def set_feature(self,
                    close_encounter_check: bool | None = None,
                    general_relativity: bool | None = None,
                    collision_model: Literal["MERGE","BOUNCE","FRAGGLE"] | None = None,
                    minimum_fragment_gmass: float | None = None,
                    minimum_fragment_mass: float | None = None,
                    nfrag_reduction: float | None = None,
                    rotation: bool | None = None,
                    compute_conservation_values: bool | None = None,
                    extra_force: bool | None = None,
                    big_discard: bool | None = None,
                    rhill_present: bool | None = None,
                    restart: bool | None = None,
                    tides: bool | None = None,
                    interaction_loops: Literal["TRIANGULAR", "FLAT"] | None = None,
                    encounter_check_loops: Literal["TRIANGULAR", "SORTSWEEP"] | None = None,
                    encounter_save: Literal["NONE", "TRAJECTORY", "CLOSEST", "BOTH"] | None = None,
                    coarray: bool | None = None,
                    verbose: bool | None = None,
                    simdir: str | os.PathLike = None, 
                    **kwargs: Any
                    ):
        """
        Turns on or off various features of a simulation.

        Parameters
        ----------
        close_encounter_check : bool, optional
            Check for close encounters between bodies. If set to True, then the radii of massive bodies must be included
            in initial conditions.
        encounter_save : {"NONE","TRAJECTORY","CLOSEST","BOTH"}, default "NONE"
            Indicate if and how encounter data should be saved. If set to "TRAJECTORY", the position and velocity vectors
            of all bodies undergoing close encounter are saved at each intermediate step to the encounter files.
            If set to "CLOSEST", the position  and velocities at the point of closest approach between pairs of bodies are 
            computed and stored to the encounter files. If set to "BOTH", then this stores the values that would be computed
            in "TRAJECTORY" and "CLOSEST". If set to "NONE" no trajectory information is saved.
            *WARNING*: Enabling this feature could lead to very large files.
        general_relativity : bool, optional
            Include the post-Newtonian correction in acceleration calculations.
        collision_model: {"MERGE","BOUNCE","FRAGGLE"}, default "MERGE"
            This is used to set the collision/fragmentation model. [TODO: DESCRIBE THESE] 
            This argument only applies to Swiftest-SyMBA simulations. It will be ignored otherwise.
            Parameter input file equivalent: `COLLISION_MODEL`
        minimum_fragment_gmass : float, optional
            If fragmentation is turned on, this sets the mimimum G*mass of a collisional fragment that can be generated if a 
            fragmentation model is enabled. Ignored otherwise.
            *Note.* Only set one of minimum_fragment_gmass or minimum_fragment_mass
            Parameter input file equivalent: None
        minimum_fragment_mass : float, optional
            If fragmentation is turned on, this sets the mimimum mass of a collisional fragment that can be generated. if a 
            fragmentation model is enabled. Ignored otherwise
            *Note.* Only set one of minimum_fragment_gmass or minimum_fragment_mass
            Parameter input file equivalent: `MIN_GMFRAG`
        nfrag_reduction : float, optional
            If fragmentation is turne don, this is a reduction factor used to limit the number of fragments generated in a collision.
            For instance, if the SFD of the collision would generated 300 fragments above the `minimum_fragment_mass`, then a value
            of `nfrag_reduction = 30.0` would reduce it to 10.  
            *Note.* Currently only used by the Fraggle collision model. 
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
        interaction_loops : {"TRIANGULAR","FLAT"}, default "TRIANGULAR"
            *Swiftest Experimental feature*
            Specifies which algorithm to use for the computation of body-body gravitational forces.
            * "TRIANGULAR" : Upper-triangular double-loops .
            * "FLAT" : Body-body interation pairs are flattened into a 1-D array.
        encounter_check_loops : {"TRIANGULAR","SORTSWEEP"}, default "TRIANGULAR"
            *Swiftest Experimental feature*
            Specifies which algorithm to use for checking whether bodies are in a close encounter state or not.
            * "TRIANGULAR" : Upper-triangular double-loops.
            * "SORTSWEEP" : A Sort-Sweep algorithm is used to reduce the population of potential close encounter bodies.
              This algorithm is still in development, and does not necessarily speed up the encounter checking.
              Use with caution.
        coarray : bool, default False
            If true, will employ Coarrays on test particle structures to run in single program/multiple data parallel mode. 
            *Note" In order to use this capability, Swiftest must be compiled for Coarray support. Only certain integrators
            can use Coarrays: RMVS, WHM, Helio are all compatible, but SyMBA is not, due to the way tp-pl close encounters 
            are handeled.           
        tides: bool, optional
            Turns on tidal model (IN DEVELOPMENT - IGNORED)
        Yarkovsky: bool, optional
            Turns on Yarkovsky model (IN DEVELOPMENT - IGNORED)
        YORP: bool, optional
            Turns on YORP model (IN DEVELOPMENT - IGNORED)
        restart : bool, default False
            If true, will restart an old run. The file given by `output_file_name` must exist before the run can
            execute. If false, will start a new run. If the file given by `output_file_name` exists, it will be replaced
            when the run is executed.
        simdir: PathLike, optional
            Directory where simulation data will be stored, including the parameter file, initial conditions file, output file,
            dump files, and log files. 
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        feature_dict : dict
            A dictionary containing the requested features.

        """

        update_list = []
        if close_encounter_check is not None:
            self.param["CHK_CLOSE"] = close_encounter_check
            update_list.append("close_encounter_check")
            
        if rhill_present is not None:
            self.param["RHILL_PRESENT"] = rhill_present
            update_list.append("rhill_present")            
            
        if extra_force is not None:
            self.param["EXTRA_FORCE"] = extra_force
            update_list.append("extra_force")

        if big_discard is not None:
            if self.codename != "Swifter":
                self.param["BIG_DISCARD"] = False
            else:
                self.param["BIG_DISCARD"] = big_discard
                update_list.append("big_discard")     
                
        if simdir is not None:
            self.simdir = Path(simdir)
            if self.simdir.exists():
                if not self.simdir.is_dir():
                    msg = f"Cannot create the {self.simdir.resolve()} directory: File exists."
                    msg += "\nDelete the file or change the location of param_file"
                    raise NotADirectoryError(msg)
            self.param_file = Path(kwargs.pop("param_file","param.in"))

        if self.codename == "Swiftest": 
            if general_relativity is not None:
                self.param["GR"] = general_relativity
                update_list.append("general_relativity")

            fragmentation_models = ["FRAGGLE"]
            if collision_model is not None:
                collision_model = collision_model.upper()
                fragmentation = collision_model in fragmentation_models
                if self.codename != "Swiftest" and self.integrator != "symba" and fragmentation:
                    warnings.warn("Fragmentation is only available on Swiftest SyMBA.",stacklevel=2)
                    self.param['COLLISION_MODEL'] = "MERGE"
                else:
                    self.param['COLLISION_MODEL'] = collision_model
                    update_list.append("collision_model")
                    if fragmentation:
                        if "MIN_GMFRAG" not in self.param and minimum_fragment_mass is None and minimum_fragment_gmass is None:
                            warnings.warn("Minimum fragment mass is not set. Set it using minimum_fragment_gmass or minimum_fragment_mass",stacklevel=2)
                        else:
                            update_list.append("minimum_fragment_gmass")

            if minimum_fragment_gmass is not None and minimum_fragment_mass is not None:
                warnings.warn("Only set either minimum_fragment_mass or minimum_fragment_gmass, but not both!",stacklevel=2)

            if minimum_fragment_gmass is not None:
                self.param["MIN_GMFRAG"] = minimum_fragment_gmass
                if "minmum_fragment_gmass" not in update_list:
                    update_list.append("minimum_fragment_gmass")
            elif minimum_fragment_mass is not None:
                self.param["MIN_GMFRAG"] = minimum_fragment_mass * self.GU
                if "minimum_fragment_gmass" not in update_list:
                    update_list.append("minimum_fragment_gmass")
                    
            if nfrag_reduction is not None:
                self.param["NFRAG_REDUCTION"] = nfrag_reduction
                update_list.append("nfrag_reduction")

            if rotation is not None:
                self.param['ROTATION'] = rotation
                update_list.append("rotation")

            if self.param['COLLISION_MODEL'] == "FRAGGLE" and not self.param['ROTATION']:
                self.param['ROTATION'] = True
                update_list.append("rotation")

            if compute_conservation_values is not None:
                self.param["ENERGY"] = compute_conservation_values
                update_list.append("compute_conservation_values")

            if restart is not None:
                self.param["RESTART"] = restart
                update_list.append("restart")

            if interaction_loops is not None:
                valid_vals = ["TRIANGULAR", "FLAT"]
                interaction_loops = interaction_loops.upper()
                if interaction_loops not in valid_vals:
                    msg = f"{interaction_loops} is not a valid option for interaction loops."
                    msg += f"\nMust be one of {valid_vals}"
                    warnings.warn(msg,stacklevel=2)
                    if "INTERACTION_LOOPS" not in self.param:
                        self.param["INTERACTION_LOOPS"] = valid_vals[0]
                else:
                    self.param["INTERACTION_LOOPS"] = interaction_loops
                    update_list.append("interaction_loops")

            if encounter_check_loops is not None:
                valid_vals = ["TRIANGULAR", "SORTSWEEP"]
                encounter_check_loops = encounter_check_loops.upper()
                if encounter_check_loops not in valid_vals:
                    msg = f"{encounter_check_loops} is not a valid option for interaction loops."
                    msg += f"\nMust be one of {valid_vals}"
                    warnings.warn(msg,stacklevel=2)
                    if "ENCOUNTER_CHECK" not in self.param:
                        self.param["ENCOUNTER_CHECK"] = valid_vals[0]
                else:
                    self.param["ENCOUNTER_CHECK"] = encounter_check_loops
                    update_list.append("encounter_check_loops")

            if encounter_save is not None:
                valid_vals = ["NONE", "TRAJECTORY", "CLOSEST", "BOTH"]
                encounter_save = encounter_save.upper()
                if encounter_save not in valid_vals:
                    msg = f"{encounter_save} is not a valid option for encounter_save."
                    msg += f"\nMust be one of {valid_vals}"
                    warnings.warn(msg,stacklevel=2)
                    if "ENCOUNTER_SAVE" not in self.param:
                        self.param["ENCOUNTER_SAVE"] = valid_vals[0]
                else:
                    self.param["ENCOUNTER_SAVE"] = encounter_save
                    update_list.append("encounter_save")
        
            if coarray is not None:
                if self.codename == "Swiftest":
                    self.param["COARRAY"] = coarray
                    update_list.append("coarray")     
                    
            self.param["TIDES"] = False
                
                    
        feature_dict = self.get_feature(update_list, verbose)
        return feature_dict

    def get_feature(self, arg_list: str | List[str] | None = None, verbose: bool | None = None, **kwargs: Any):
        """

        Returns a subset of the parameter dictionary containing the current value of the feature boolean values.
        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list: str | List[str], optional
            A single string or list of strings containing the names of the features to extract. Default is all of:
            ["close_encounter_check", "general_relativity", "collision_model", "rotation", "compute_conservation_values"]
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        feature_dict : dict
           A dictionary containing the requested features.

        """

        valid_var = {"close_encounter_check": "CHK_CLOSE",
                     "collision_model": "COLLISION_MODEL",
                     "encounter_save": "ENCOUNTER_SAVE",
                     "minimum_fragment_gmass": "MIN_GMFRAG",
                     "nfrag_reduction": "NFRAG_REDUCTION",
                     "rotation": "ROTATION",
                     "general_relativity": "GR",
                     "compute_conservation_values": "ENERGY",
                     "rhill_present": "RHILL_PRESENT",
                     "extra_force": "EXTRA_FORCE",
                     "big_discard": "BIG_DISCARD",
                     "interaction_loops": "INTERACTION_LOOPS",
                     "encounter_check_loops": "ENCOUNTER_CHECK",
                     "coarray" : "COARRAY",
                     "restart": "RESTART"
                     }

        valid_arg, feature_dict = self._get_valid_arg_list(arg_list, valid_var)

        if verbose is None:
            verbose = self.verbose

        if verbose:
            for arg in valid_arg:
                key = valid_var[arg]
                if key in feature_dict:
                    if arg == "minimum_fragment_gmass":
                        print(f"{arg:<{self._getter_column_width}} {feature_dict[key]} {self.DU_name}^3 / {self.TU_name}^2 ")
                        print(f"{'minimum_fragment_mass':<{self._getter_column_width}} {feature_dict[key] / self.GU} {self.MU_name}")
                    else:
                        print(f"{arg:<{self._getter_column_width}} {feature_dict[key]}")
                else:
                    print(f"{arg:<{self._getter_column_width}} NOT SET")


        return feature_dict

    def set_init_cond_files(self,
                            init_cond_file_type: Literal["NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"] | None = None,
                            init_cond_file_name: str | os.PathLike | Dict[str, str] |
                                                 Dict[ str, os.PathLike] | None = None,
                            init_cond_format: Literal["EL", "XV"] | None = None,
                            verbose: bool | None = None,
                            **kwargs: Any
                            ):
        """
        Sets the initial condition file parameters in the parameters dictionary.

        Parameters
        ----------
        init_cond_file_type : {"NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"}, optional
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
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        init_cond_file_dict : dict
           A dictionary containing the requested parameters

        """

        update_list = []
        if init_cond_file_name is not None:
            update_list.append("init_cond_file_name")
        if init_cond_file_type is not None:
            update_list.append("init_cond_file_type")
        if init_cond_format is not None:
            update_list.append("init_cond_format")

        if len(update_list) == 0:
            return {}

        def ascii_file_input_error_msg(codename):
            msg = f"in set_init_cond_files: init_cond_file_name must be a dictionary of the form: "
            msg += "\n {"
            if codename == "Swiftest":
                msg += '\n"CB" : *path to central body initial conditions file*,'
            msg += '\n"PL" : *path to massive body initial conditions file*,'
            msg += '\n"TP" : *path to test particle initial conditions file*'
            msg += '\n}'
            warnings.warn(msg,stacklevel=2)
            return {}

        if init_cond_format is None:
            if "IN_FORM" in self.param:
                init_cond_format = self.param['IN_FORM']
            else:
                if self.codename.title() == "Swiftest":
                    init_cond_format = "EL"
                else:
                    init_cond_format = "XV"

        if init_cond_file_type is None:
            if "IN_TYPE" in self.param:
                init_cond_file_type = self.param['IN_TYPE']
            else:
                if self.codename.title() == "Swiftest":
                    init_cond_file_type = "NETCDF_DOUBLE"
                else:
                    init_cond_file_type = "ASCII"

        if self.codename.title() == "Swiftest":
            init_cond_keys = ["CB", "PL", "TP"]
        else:
            init_cond_keys = ["PL", "TP"]
            if init_cond_file_type != "ASCII":
                warnings.warn(f"{init_cond_file_type} is not supported by {self.codename}. Using ASCII instead",stacklevel=2)
                init_cond_file_type = "ASCII"
            if init_cond_format != "XV":
                warnings.warn(f"{init_cond_format} is not supported by {self.codename}. Using XV instead",stacklevel=2)
                init_cond_format = "XV"

        valid_formats = {"EL", "XV"}
        if init_cond_format not in valid_formats:
            warnings.warn(f"{init_cond_format} is not a valid input format",stacklevel=2)
        else:
            self.param['IN_FORM'] = init_cond_format

        valid_types = {"NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"}
        if init_cond_file_type not in valid_types:
            warnings.warn(f"{init_cond_file_type} is not a valid input type",stackevel=2)
        else:
            self.param['IN_TYPE'] = init_cond_file_type

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
                warnings.warn(f"Only a single input file is used for NetCDF files",stacklevel=2)
            else:
                self.param["NC_IN"] = init_cond_file_name

        init_cond_file_dict = self.get_init_cond_files(update_list, verbose)

        return init_cond_file_dict

    def get_init_cond_files(self, arg_list: str | List[str] | None = None, verbose: bool | None = None, **kwargs):
        """

        Returns a subset of the parameter dictionary containing the current initial condition file parameters
        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list : str | List[str], optional
            A single string or list of strings containing the names of the simulation time parameters to extract.
            Default is all of:
            ["init_cond_file_type", "init_cond_file_name", "init_cond_format"]
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.


        Returns
        -------
        init_cond_file_dict : dict
           A dictionary containing the requested parameters

        """

        valid_var = {"init_cond_file_type": "IN_TYPE",
                     "init_cond_format": "IN_FORM",
                     "init_cond_file_name": "NC_IN",
                     "init_cond_file_name['CB']": "CB_IN",
                     "init_cond_file_name['PL']": "PL_IN",
                     "init_cond_file_name['TP']": "TP_IN",
                     }

        three_file_args = ["init_cond_file_name['CB']",
                           "init_cond_file_name['PL']",
                           "init_cond_file_name['TP']"]

        if self.codename == "Swifter":
            three_file_args.remove("init_cond_file_name['CB']")
            valid_var.pop("init_cond_file_name['CB']",None)

        # We have to figure out which initial conditions file model we are using (1 vs. 3 files)
        if arg_list is None:
            valid_arg = list(valid_var.keys())
        elif type(arg_list) is str:
            valid_arg = [arg_list]
        else:
            valid_arg = [k for k in arg_list if k in list(valid_var.keys())]

        # Figure out which input file model we need to use
        if "init_cond_file_name" in valid_arg:
            if self.param["IN_TYPE"] == "ASCII":
                valid_arg.remove("init_cond_file_name")
                for key in three_file_args:
                    if key not in valid_arg:
                        valid_arg.append(key)
            else:
                for key in three_file_args:
                    if key in valid_arg:
                        valid_arg.remove(key)

        valid_arg, init_cond_file_dict = self._get_valid_arg_list(valid_arg, valid_var)

        if verbose is None:
            verbose = self.verbose

        if verbose:
            for arg in valid_arg:
                key = valid_var[arg]
                print(f"{arg:<{self._getter_column_width}} {init_cond_file_dict[key]}")

        return init_cond_file_dict

    def set_output_files(self,
                         output_file_type: Literal[
                                               "NETCDF_DOUBLE", "NETCDF_FLOAT", "REAL4", "REAL8", "XDR4", "XDR8"] | None = None,
                         output_file_name: os.PathLike | str | None = None,
                         output_format: Literal["XV", "XVEL"] | None = None,
                         restart: bool | None = None,
                         verbose: bool | None = None,
                         **kwargs: Any
                         ):
        """
        Sets the output file parameters in the parameter dictionary.

        Parameters
        ----------
        output_file_type : {"NETCDF_DOUBLE", "NETCDF_FLOAT","REAL4","REAL8","XDR4","XDR8"}, optional
            The file type for the outputs of the simulation. Compatible file types depend on the `codename` argument.
            * Swiftest: Only "NETCDF_DOUBLE" or "NETCDF_FLOAT" supported.
            * Swifter: Only "REAL4","REAL8","XDR4" or "XDR8"  supported.
            * Swift: Only "REAL4" supported.
        output_file_name : str or path-like, optional
            Name of output file to generate. If not supplied, then one of the default file names are used, depending on
            the value passed to `output_file_type`. If one of the NetCDF types are used, the default is "data.nc".
            Otherwise, the default is "bin.dat".
        output_format : {"XV","XVEL"}, optional
            Specifies the format for the data saved to the output file. If "XV" then cartesian position and velocity
            vectors for all bodies are stored. If "XVEL" then the orbital elements are also stored.
        restart: bool, optional
            Indicates whether this is a restart of an old run or a new run.
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        output_file_dict : dict
           A dictionary containing the requested parameters

        """
        update_list = []
        if output_file_type is not None:
            update_list.append("output_file_type")
        if output_file_name is not None:
            update_list.append("output_file_name")
        if output_format is not None:
            update_list.append("output_format")
        if restart is not None:
            self.restart = restart
            update_list.append("restart")
        if len(update_list) == 0:
            return {}

        if self.codename == "Swiftest":
            if output_file_type is None:
                output_file_type = self.param.pop("OUT_TYPE", None)
                if output_file_type is None:
                    output_file_type = "NETCDF_DOUBLE"
            elif output_file_type not in ["NETCDF_DOUBLE", "NETCDF_FLOAT"]:
                warnings.warn(f"{output_file_type} is not compatible with Swiftest. Setting to NETCDF_DOUBLE",stacklevel=2)
                output_file_type = "NETCDF_DOUBLE"
        elif self.codename == "Swifter":
            if output_file_type is None:
                output_file_type = self.param.pop("OUT_TYPE", None)
                if output_file_type is None:
                    output_file_type = "REAL8"
            elif output_file_type not in ["REAL4", "REAL8", "XDR4", "XDR8"]:
                warnings.warn(f"{output_file_type} is not compatible with Swifter. Setting to REAL8",stacklevel=2)
                output_file_type = "REAL8"
        elif self.codename == "Swift":
            if output_file_type is None:
                output_file_type = self.param.pop("OUT_TYPE", None)
                if output_file_type is None:
                    output_file_type = "REAL4"
            if output_file_type not in ["REAL4"]:
                warnings.warn(f"{output_file_type} is not compatible with Swift. Setting to REAL4",stacklevel=2)
                output_file_type = "REAL4"

        self.param['OUT_TYPE'] = output_file_type
        if output_file_name is None:
            if output_file_type in ["NETCDF_DOUBLE", "NETCDF_FLOAT"]:
                self.param['BIN_OUT'] = "data.nc"
            else:
                self.param['BIN_OUT'] = "bin.dat"
        else:
            self.param['BIN_OUT'] = output_file_name

        if output_format is not None:
            if output_format != "XV" and self.codename != "Swiftest":
                warnings.warn(f"{output_format} is not compatible with {self.codename}. Setting to XV",stacklevel=2)
                output_format = "XV"
            self.param["OUT_FORM"] = output_format

        if self.restart:
            self.param["OUT_STAT"] = "APPEND"
        else:
            self.param["OUT_STAT"] = "REPLACE"

        output_file_dict = self.get_output_files(update_list, verbose=verbose)

        return output_file_dict

    def get_output_files(self, arg_list: str | List[str] | None = None, verbose: bool | None = None, **kwargs):
        """

        Returns a subset of the parameter dictionary containing the current output file parameters
        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list : str | List[str], optional
            A single string or list of strings containing the names of the simulation time parameters to extract.
            Default is all of:
            ["output_file_type", "output_file_name", "output_format"]
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.


        Returns
        -------
        output_file_dict : dict
           A dictionary containing the requested parameters

        """

        valid_var = {"output_file_type": "OUT_TYPE",
                     "output_file_name": "BIN_OUT",
                     "output_format": "OUT_FORM",
                     "restart": "OUT_STAT"
                     }

        valid_arg, output_file_dict = self._get_valid_arg_list(arg_list, valid_var)

        if verbose is None:
            verbose = self.verbose

        if verbose:
            for arg in valid_arg:
                key = valid_var[arg]
                print(f"{arg:<{self._getter_column_width}} {output_file_dict[key]}")

        return output_file_dict

    def set_unit_system(self,
                        MU: str | None = None,
                        DU: str | None = None,
                        TU: str | None = None,
                        MU2KG: float | None = None,
                        DU2M: float | None = None,
                        TU2S: float | None = None,
                        MU_name: str | None = None,
                        DU_name: str | None = None,
                        TU_name: str | None = None,
                        recompute_unit_values: bool = True,
                        verbose: bool | None = None,
                        **kwargs: Any):
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
        MU_name : str, optional
            The name of the mass unit. When setting one of the standard units via `MU` a name will be
            automatically set for the unit, so this argument will override the automatic name.
        DU_name : str, optional
            The name of the distance unit. When setting one of the standard units via `DU` a name will be
            automatically set for the unit, so this argument will override the automatic name.
        TU_name : str, optional
            The name of the time unit. When setting one of the standard units via `TU` a name will be
            automatically set for the unit, so this argument will override the automatic name.
        recompute_unit_values : bool, default True
            Recompute all values into the new unit system.
            >*Note:* This is a destructive operation, however if not executed then the values contained in the parameter
            > file and input/output data files computed previously may not be consistent with the new unit conversion
            > factors.
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        ----------
        unit_dict : dict
           A dictionary containing the requested unit conversion parameters
        """

        MU2KG_old = None
        DU2M_old = None
        TU2S_old = None

        if "MU_name" not in dir(self):
            self.MU_name = "MU"
        if "DU_name" not in dir(self):
            self.DU_name = "DU"
        if "TU_name" not in dir(self):
            self.TU_name = "TU"

        update_list = []
        if MU is not None or MU2KG is not None:
            update_list.append("MU")
        if DU is not None or DU2M is not None:
            update_list.append("DU")
        if TU is not None or TU2S is not None:
            update_list.append("TU")

        if MU2KG is not None or MU is not None:
            MU2KG_old = self.param.pop('MU2KG', None)
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
                    warnings.warn(f"{MU} not a recognized unit system. Using MSun as a default.",stacklevel=2)
                    self.param['MU2KG'] = constants.MSun
                    self.MU_name = "MSun"
            self.MU2KG = self.param['MU2KG']
            self.KG2MU = 1.0 / self.MU2KG

        if DU2M is not None or DU is not None:
            DU2M_old = self.param.pop('DU2M', None)
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
                    warnings.warn(f"{DU} not a recognized unit system. Using AU as a default.",stacklevel=2)
                    self.param['DU2M'] = constants.AU2M
                    self.DU_name = "AU"
            self.DU2M = self.param['DU2M']
            self.M2DU = 1.0 / self.DU2M

        if TU2S is not None or TU is not None:
            TU2S_old = self.param.pop('TU2S', None)
            if TU2S is not None:
                self.param['TU2S'] = TU2S
                self.TU_name = None
            else:
                if TU.upper() == "YR" or TU.upper() == "Y" or TU.upper() == "YEAR" or TU.upper() == "YEARS":
                    self.param['TU2S'] = constants.YR2S
                    self.TU_name = "y"
                elif TU.upper() == "DAY" or TU.upper() == "D" or TU.upper() == "JD" or TU.upper() == "DAYS":
                    self.param['TU2S'] = constants.JD2S
                    self.TU_name = "d"
                elif TU.upper() == "S" or TU.upper() == "SECONDS" or TU.upper() == "SEC":
                    self.param['TU2S'] = 1.0
                    self.TU_name = "s"
                else:
                    warnings.warn(f"{TU} not a recognized unit system. Using YR as a default.",stacklevel=2)
                    self.param['TU2S'] = constants.YR2S
                    self.TU_name = "y"
            self.TU2S = self.param['TU2S']
            self.S2TU = 1.0 / self.TU2S

        if MU_name is not None:
            self.MU_name = MU_name
        if DU_name is not None:
            self.DU_name = DU_name
        if TU_name is not None:
            self.TU_name = TU_name

        if "DU_name" in dir(self) and "TU_name" in dir(self):
            self.VU_name = f"{self.DU_name}/{self.TU_name}"
        if all(key in self.param for key in ["MU2KG","DU2M","TU2S"]):
            self.GU = constants.GC * self.param["TU2S"] ** 2 * self.param["MU2KG"] / self.param["DU2M"] ** 3

        if "MU2KG" in self.param and "DU2M" in self.param and "TU2S" in self.param:
            if recompute_unit_values and \
                    MU2KG_old != self.param['MU2KG'] or \
                    DU2M_old != self.param['DU2M'] or \
                    TU2S_old != self.param['TU2S']:
                self.update_param_units(MU2KG_old, DU2M_old, TU2S_old)

        unit_dict = self.get_unit_system(update_list, verbose)

        return unit_dict

    def get_unit_system(self, arg_list: str | List[str] | None = None, verbose: bool | None = None, **kwargs):
        """

        Returns a subset of the parameter dictionary containing the current simulation unit system.
        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list : str | List[str], optional
            A single string or list of strings containing the names of the simulation unit system
            Default is all of:
            ["MU", "DU", "TU"]
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        unit_dict : dict
           A dictionary containing the requested unit conversion parameters

        """

        valid_var = {
            "MU": "MU2KG",
            "DU": "DU2M",
            "TU": "TU2S",
        }

        if "MU_name" not in dir(self) or self.MU_name is None:
            MU_name = "mass unit"
        else:
            MU_name = self.MU_name
        if "DU_name" not in dir(self) or self.DU_name is None:
            DU_name = "distance unit"
        else:
            DU_name = self.DU_name
        if "TU_name" not in dir(self) or self.TU_name is None:
            TU_name = "time unit"
        else:
            TU_name = self.TU_name

        units1 = {
            "MU": MU_name,
            "DU": DU_name,
            "TU": TU_name
        }
        units2 = {
            "MU": f"kg / {MU_name}",
            "DU": f"m / {DU_name}",
            "TU": f"s / {TU_name}"
        }

        valid_arg, unit_dict = self._get_valid_arg_list(arg_list, valid_var)

        if verbose is None:
            verbose = self.verbose

        if verbose:
            for arg in valid_arg:
                key = valid_var[arg]
                col_width = str(int(self._getter_column_width) - 4)
                print(f"{arg}: {units1[arg]:<{col_width}} {unit_dict[key]} {units2[arg]}")

        return unit_dict

    def update_param_units(self, MU2KG_old, DU2M_old, TU2S_old):
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
        distance_keys = ['CHK_QMIN', 'CHK_RMIN', 'CHK_RMAX', 'CHK_EJECT']
        time_keys = ['T0', 'TSTOP', 'DT']

        if MU2KG_old is not None:
            for k in mass_keys:
                if k in self.param:
                    self.param[k] *= MU2KG_old / self.param['MU2KG']

        if DU2M_old is not None:
            for k in distance_keys:
                if k in self.param:
                    self.param[k] *= DU2M_old / self.param['DU2M']

            CHK_QMIN_RANGE = self.param.pop('CHK_QMIN_RANGE', None)
            if CHK_QMIN_RANGE is not None:
                CHK_QMIN_RANGE = CHK_QMIN_RANGE.split(" ")
                for i, v in enumerate(CHK_QMIN_RANGE):
                    CHK_QMIN_RANGE[i] = float(CHK_QMIN_RANGE[i]) * self.param['DU2M'] / DU2M_old
                self.param['CHK_QMIN_RANGE'] = f"{CHK_QMIN_RANGE[0]} {CHK_QMIN_RANGE[1]}"

        if TU2S_old is not None:
            for k in time_keys:
                if k in self.param:
                    self.param[k] *= TU2S_old / self.param['TU2S']

        return

    def set_distance_range(self,
                           rmin: float | None = None,
                           rmax: float | None = None,
                           qmin_coord: Literal["HELIO","BARY"] | None = None,
                           verbose: bool | None = None,
                           **kwargs: Any):
        """
        Sets the minimum and maximum distances of the simulation.

        Parameters
        ----------
        rmin : float
            Minimum distance of the simulation (CHK_QMIN, CHK_RMIN, CHK_QMIN_RANGE[0])
        rmax : float
            Maximum distance of the simulation (CHK_RMAX, CHK_QMIN_RANGE[1])
        qmin_coord : str, {"HELIO", "BARY"}
            coordinate frame to use for CHK_QMIN
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        range_dict : dict
           A dictionary containing the requested parameters.

        """
        if rmax is None and rmin is None and qmin_coord is None:
            return {}

        update_list = []
        CHK_QMIN_RANGE = self.param.pop('CHK_QMIN_RANGE', None)
        if CHK_QMIN_RANGE is None:
            CHK_QMIN_RANGE = [-1, -1]
        else:
            CHK_QMIN_RANGE = CHK_QMIN_RANGE.split(" ")
        if rmin is not None:
            self.param['CHK_QMIN'] = rmin
            self.param['CHK_RMIN'] = rmin
            CHK_QMIN_RANGE[0] = rmin
            update_list.append("rmin")
        if rmax is not None:
            self.param['CHK_RMAX'] = rmax
            self.param['CHK_EJECT'] = rmax
            CHK_QMIN_RANGE[1] = rmax
            update_list.append("rmax")
        if qmin_coord is not None:
            valid_qmin_coord = ["HELIO","BARY"]
            if qmin_coord.upper() not in valid_qmin_coord:
                warnings.warn(f"qmin_coord = {qmin_coord} is not a valid option.  Must be one of",','.join(valid_qmin_coord),stacklevel=2)
                self.param['CHK_QMIN_COORD'] = valid_qmin_coord[0]
            else:
                self.param['CHK_QMIN_COORD'] = qmin_coord.upper()
            update_list.append("qmin_coord")

        self.param['CHK_QMIN_RANGE'] = f"{CHK_QMIN_RANGE[0]} {CHK_QMIN_RANGE[1]}"

        range_dict = self.get_distance_range(update_list, verbose=verbose)

        return range_dict

    def get_distance_range(self, arg_list: str | List[str] | None = None, verbose: bool | None = None, **kwargs: Any):
        """

        Returns a subset of the parameter dictionary containing the current values of the distance range parameters.
        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list: str | List[str], optional
            A single string or list of strings containing the names of the features to extract. Default is all of:
            ["rmin", "rmax"]
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        range_dict : dict
           A dictionary containing the requested parameters.

        """

        valid_var = {"rmin": "CHK_RMIN",
                     "rmax": "CHK_RMAX",
                     "qmin_coord": "CHK_QMIN_COORD",
                     "qmin": "CHK_QMIN",
                     "qminR": "CHK_QMIN_RANGE"
                     }

        units = {"rmin": self.DU_name,
                 "rmax": self.DU_name,
                 "qmin": self.DU_name,
                 "qminR": self.DU_name,
                 }

        if type(arg_list) is str:
            arg_list = [arg_list]
        if arg_list is not None:
            if "rmin" in arg_list:
                arg_list.append("qmin")
            if "rmax" in arg_list or "rmin" in arg_list:
                arg_list.append("qminR")

        valid_arg, range_dict = self._get_valid_arg_list(arg_list, valid_var)

        if verbose is None:
            verbose = self.verbose

        if verbose:
            if "rmin" in valid_arg:
                key = valid_var["rmin"]
                print(f"{'rmin':<{self._getter_column_width}} {range_dict[key]} {units['rmin']}")
            if "rmax" in valid_arg:
                key = valid_var["rmax"]
                print(f"{'rmax':<{self._getter_column_width}} {range_dict[key]} {units['rmax']}")
            if "qmin_coord" in valid_arg:
                key = valid_var["qmin_coord"]
                print(f"{'qmin_coord':<{self._getter_column_width}} {range_dict[key]}")

        return range_dict

    def add_solar_system_body(self,
                              name: str | List[str] | None = None,
                              ephemeris_id: int | List[int] | None = None,
                              date: str | None = None,
                              source: str = "HORIZONS", 
                              **kwargs: Any):
        """
        Adds a solar system body to an existing simulation Dataset from the JPL Horizons ephemeris service. The JPL Horizons service
        will be searched for a body matching the string passed by `name`, or alternatively `ephemeris_id` if passed. Bodies will be
        named in the Swiftest initial conditions Dataset using `name`. Use `ephemeris_id` to have finer control over which body is
        searched in Horizons while using a custom name.
        
        If `name` is not passed, then the target name property is used as the name. You must pass either `name` and/or `ephemeris_id`
        
        When passing `name` == "Earth" or `name` == "Pluto", it a body is generated that has initial conditions matching the system
        barycenter and mass equal to the sum of Earth+Moon or Pluto+Charon. 
        
        To obtain initial conditions for either Earth or Pluto alone, pass `ephemeris_id` == "399" for Earth or 
        `ephemeris_id` == "999" for Pluto.  


        Parameters
        ----------
        name : str, optional | List[str]
            Add solar system body by name. This will be the name used in the Swiftest initial conditions Dataset unless not supplied 
        ephemeris_id : int | List[int], optional but must be the same length as `name` if passed.
            In that case, the id is passed to the ephemeris service and the name is used. 
        date : str, optional
            ISO-formatted date sto use when obtaining the ephemerides in the format YYYY-MM-DD. Defaults to value
            set by `set_ephemeris_date`.
        source : str, default "Horizons"
            The source of the ephemerides.
            >*Note.* Currently only the JPL Horizons ephemeris is implemented, so this is ignored.
        **kwargs: Any
            Additional keyword arguments to pass to the query method (i.e. astroquery.Horizons)
        Returns
        -------
        None
            initial conditions data stored as an Xarray Dataset in the init_cond instance variable
        """
        
        if name == None and ephemeris_id == None:
            warnings.warn("Either `name` and/or `ephemeris_id` must be supplied to add_solar_system_body")
            return None

        if ephemeris_id is not None:
            if type(ephemeris_id) is int or type(ephemeris_id) is str:
                ephemeris_id = [ephemeris_id]
            if name is None:
                name = [None] * len(ephemeris_id)
            elif len(ephemeris_id) != len(name):
                warnings.warn(f"The length of ephemeris_id ({len(ephemeris_id)}) does not match the length of name ({len(name)})",stacklevel=2)
                return None
        else:
            if type(name) is str:
                name = [name]
            ephemeris_id = [None] * len(name)

        if self.ephemeris_date is None:
            self.set_ephemeris_date()

        if date is None:
            date = self.ephemeris_date
        try:
            datetime.datetime.fromisoformat(date)
        except:
            warnings.warn(f"{date} is not a valid date format. Must be 'YYYY-MM-DD'. Setting to {self.ephemeris_date}",stacklevel=2)
            date = self.ephemeris_date

        if source.upper() != "HORIZONS":
            warnings.warn("Currently only the JPL Horizons ephemeris service is supported",stacklevel=2)

        body_list = []
        for i,n in enumerate(name):
            body = init_cond.solar_system_horizons(n, self.param, date, ephemeris_id=ephemeris_id[i],**kwargs)
            if body is not None:
                body_list.append(body)

        #Convert the list receieved from the solar_system_horizons output and turn it into arguments to vec2xr
        if len(body_list) == 0:
           print("No valid bodies found")
           return 
        elif len(body_list) == 1:
            values = list(np.hsplit(np.array(body_list[0],dtype=np.dtype(object)),17))
        else:
            values = list(np.squeeze(np.hsplit(np.array(body_list,np.dtype(object)),17)))
        keys = ["id","name","a","e","inc","capom","omega","capm","rh","vh","Gmass","radius","rhill","Ip","rot","j2rp2","j4rp4"]
        kwargs = dict(zip(keys,values))
        scalar_floats = ["a","e","inc","capom","omega","capm","Gmass","radius","rhill","j2rp2","j4rp4"]
        vector_floats = ["rh","vh","Ip","rot"]
        scalar_ints = ["id"]

        for k,v in kwargs.items():
            if k in scalar_ints:
                v[v == None] = -1
                kwargs[k] = v.astype(int)
            elif k in scalar_floats:
                kwargs[k] = v.astype(np.float64)
                if all(np.isnan(kwargs[k])):
                    kwargs[k] = None
            elif k in vector_floats:
                kwargs[k] = np.vstack(v)
                kwargs[k] = kwargs[k].astype(np.float64)
                if np.all(np.isnan(kwargs[k])):
                    kwargs[k] = None

        kwargs['time'] = np.array([self.param['TSTART']])
        
        if len(self.data) == 0:
            maxid = -1
        else:
            maxid = self.data.id.max().values[()]

        nbodies = kwargs["name"].size
        kwargs['id'] = np.where(kwargs['id'] < 0,np.arange(start=maxid+1,stop=maxid+1+nbodies,dtype=int),kwargs['id'])
        
        dsnew = init_cond.vec2xr(self.param,**kwargs)

        dsnew = self._combine_and_fix_dsnew(dsnew)
        if dsnew['id'].max(dim='name') > 0 and dsnew['name'].size > 0:
           self.save(verbose=False)

        self.init_cond = self.data.copy(deep=True)

        return

    def set_ephemeris_date(self,
                          ephemeris_date: str | None = None,
                          verbose: bool | None = None,
                          **kwargs: Any):
        """

        Parameters
        ----------
        ephemeris_date : str, optional
            Date to use when obtaining the ephemerides.
            Valid options are "today", "MBCL", or date in the format YYYY-MM-DD.
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        Sets the `ephemeris_date` instance variable.

        """

        if ephemeris_date is None:
            return

        # The default value is Prof. Minton's Brimley/Cocoon line crossing date (aka MBCL)
        minton_bday = datetime.date.fromisoformat('1976-08-05')
        brimley_cocoon_line = datetime.timedelta(days=18530)
        minton_bcl = (minton_bday + brimley_cocoon_line).isoformat()

        if ephemeris_date is None or ephemeris_date.upper() == "MBCL":
            ephemeris_date = minton_bcl
        elif ephemeris_date.upper() == "TODAY":
            ephemeris_date = datetime.date.today().isoformat()
        else:
            try:
                datetime.datetime.fromisoformat(ephemeris_date)
            except:
                valid_date_args = ['"MBCL"', '"TODAY"', '"YYYY-MM-DD"']
                msg = f"{ephemeris_date} is not a valid format. Valid options include:", ', '.join(valid_date_args)
                msg += "\nUsing MBCL for date."
                warnings.warn(msg,stacklevel=2)
                ephemeris_date = minton_bcl

        self.ephemeris_date = ephemeris_date

        ephemeris_date = self.get_ephemeris_date(verbose=verbose)

        return ephemeris_date

    def get_ephemeris_date(self, arg_list: str | List[str] | None = None, verbose: bool | None = None, **kwargs: Any):
        """

        Prints the current value of the ephemeris date

        Parameters
        ----------
        arg_list: str | List[str], optional
            A single string or list of strings containing the names of the features to extract. Default is all of:
            ["integrator"]
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        ephemeris_date: str
            The ISO-formatted date string for the ephemeris computation

        """

        try:
            self.ephemeris_date
        except:
            warnings.warn(f"ephemeris_date is not set",stacklevel=2)
            return

        valid_arg = {"ephemeris_date": self.ephemeris_date}

        ephemeris_date = self._get_instance_var(arg_list, valid_arg,verbose, **kwargs)

        return ephemeris_date

    def _get_instance_var(self, arg_list: str | List[str], valid_arg: Dict, verbose: bool | None = None, **kwargs: Any):
        """

        Prints the current value of an instance variable.

        Parameters
        ----------
        arg_list: str | List[str]
            A single string or list of strings containing the names of the the instance variable to get.
        valid_arg: dict
            A dictionary where the key is the parameter argument and the value is the equivalent instance variable value.
        verbose: bool, optional
            If passed, it will override the Simulation object's verbose flag
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        Tuple of instance variable values given by the arg_list

        """

        arg_vals = []
        if verbose is None:
            verbose = self.verbose
        if verbose:
            if arg_list is None:
                arg_list = list(valid_arg.keys())
            for arg in arg_list:
                if arg in valid_arg:
                    print(f"{arg:<{self._getter_column_width}} {valid_arg[arg]}")
                    arg_vals.append(valid_arg[arg])

        return tuple(arg_vals)

    def add_body(self,
                 name: str | List[str] | npt.NDArray[np.str_] | None=None,
                 id: int | list[int] | npt.NDArray[np.int_] | None=None,
                 a: float | List[float] | npt.NDArray[np.float_] | None = None,
                 e: float | List[float] | npt.NDArray[np.float_] | None = None,
                 inc: float | List[float] | npt.NDArray[np.float_] | None = None,
                 capom: float | List[float] | npt.NDArray[np.float_] | None = None,
                 omega: float | List[float] | npt.NDArray[np.float_] | None = None,
                 capm: float | List[float] | npt.NDArray[np.float_] | None = None,
                 rh: List[float] | List[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
                 vh: List[float] | List[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
                 mass: float | List[float] | npt.NDArray[np.float_] | None=None,
                 Gmass: float | List[float] | npt.NDArray[np.float_] | None=None,
                 radius: float | List[float] | npt.NDArray[np.float_] | None=None,
                 rhill: float | List[float] | npt.NDArray[np.float_] | None=None,
                 rot: List[float] | List[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None=None,
                 Ip: List[float] | npt.NDArray[np.float_] | None=None,
                 J2: float | List[float] | npt.NDArray[np.float_] | None=None,
                 J4: float | List[float] | npt.NDArray[np.float_] | None=None):
        """
        Adds a body (test particle or massive body) to the internal DataSet given a set up 6 vectors (orbital elements
        or cartesian state vectors, depending on the value of self.param). Input all angles in degress.

        This method will update self.data with the new body or bodies added to the existing Dataset.

        Parameters
        ----------
        name : str or array-like of str, optional
        Name or names of Bodies. If none passed, name will be "Body<idval>"
        id : int or array-like of int, optional
        Unique id values. If not passed, an id will be assigned in ascending order starting from the pre-existing
        Dataset ids.
        a : float or array-like of float, optional
        semimajor axis for param['IN_FORM'] == "EL"
        e : float or array-like of float, optional
        eccentricity  for param['IN_FORM'] == "EL"
        inc : float or array-like of float, optional
        inclination for param['IN_FORM'] == "EL"
        capom : float or array-like of float, optional
        longitude of ascending node for param['IN_FORM'] == "EL"
        omega : float or array-like of float, optional
        argument of periapsis for param['IN_FORM'] == "EL"
        capm : float or array-like of float, optional
        mean anomaly for param['IN_FORM'] == "EL"
        rh : (n,3) array-like of float, optional
        Position vector array.
        vh : (n,3) array-like of float, optional
        Velocity vector array.
        mass : float or array-like of float, optional
        mass values if these are massive bodies (only one of mass or Gmass can be passed)
        Gmass : float or array-like of float, optional
        G*mass values if these are massive bodies (only one of mass or Gmass can be passed)
        radius : float or array-like of float, optional
        Radius values if these are massive bodies
        rhill : float or array-like of float, optional
        Hill's radius values if these are massive bodies
        rot: (3) or (n,3) array-like of float, optional
        Rotation rate vectors if these are massive bodies with rotation enabled.
        Ip: (3) or (n,3) array-like of float, optional
        Principal axes moments of inertia vectors if these are massive bodies with rotation enabled.

        Returns
        -------
        data : Xarray Dataset
        Dasaset containing the body or bodies that were added

        """

        #convert all inputs to numpy arrays
        def input_to_array(val,t,n=None):
            if t == "f":
                t = np.float64
            elif t == "i":
                t = np.int64
            elif t == "s":
                t = str

            if val is None:
                return None, n
            elif isinstance(val, np.ndarray):
                pass
            elif np.isscalar(val):
                val = np.array([val],dtype=t)
            else:
                try:
                    val = np.array(val,dtype=t)
                except:
                    raise ValueError(f"{val} cannot be converted to a numpy array")

            if n is None:
                return val, len(val)
            else:
                if n != len(val):
                    raise ValueError(f"Mismatched array lengths in add_body. Got {len(val)} when expecting {n}")
                return val, n

        def input_to_array_3d(val,n=None):
            if val is None:
                return None, n
            elif isinstance(val, np.ndarray):
                pass
            else:
                try:
                    val = np.array(val,dtype=np.float64)
                except:
                    raise ValueError(f"{val} cannot be converted to a numpy array")
                if n is None:
                    ndims = len(val.shape)
                    if ndims > 2 or ndims == 0:
                        raise ValueError(f"Argument must be an (n,3) or (3,) array. This one is {val.shape}")
                    else:
                        if val.shape[-1] != 3:
                            raise ValueError(f"Argument must be a 3-dimensional vector. This one has {val.shape[0]}!")
                        if val.dim == 1:
                            n = 1
                        else:
                            n = val.shape[0]
                elif n == 1:
                    if val.shape != (1,3) and val.shape != (3,):
                        raise ValueError(f"Argument is an incorrect shape. Expected {(n,3)} or {(3,1)}. Got {val.shape} instead")
                    elif val.shape == (3,):
                        val = np.expand_dims(val,axis=0)
                elif val.shape != (n,3) and val.shape != (3,n):
                    raise ValueError(f"Argument is an incorrect shape. Expected {(n,3)} or {(3,n)}. Got {val.shape} instead")
                elif val.shape == (3,n):
                    val = val.T

            return val, n

        nbodies = None
        name,nbodies = input_to_array(name,"s",nbodies)
        a,nbodies = input_to_array(a,"f",nbodies)
        e,nbodies = input_to_array(e,"f",nbodies)
        inc,nbodies = input_to_array(inc,"f",nbodies)
        capom,nbodies = input_to_array(capom,"f",nbodies)
        omega,nbodies = input_to_array(omega,"f",nbodies)
        capm,nbodies = input_to_array(capm,"f",nbodies)
        id,nbodies = input_to_array(id,"i",nbodies)
        mass,nbodies = input_to_array(mass,"f",nbodies)
        Gmass,nbodies = input_to_array(Gmass,"f",nbodies)
        rhill,nbodies = input_to_array(rhill,"f",nbodies)
        radius,nbodies = input_to_array(radius,"f",nbodies)
        J2,nbodies = input_to_array(J2,"f",nbodies)
        J4,nbodies = input_to_array(J4,"f",nbodies)

        rh,nbodies = input_to_array_3d(rh,nbodies)
        vh,nbodies = input_to_array_3d(vh,nbodies)
        rot,nbodies = input_to_array_3d(rot,nbodies)
        Ip,nbodies = input_to_array_3d(Ip,nbodies)

        if len(self.data) == 0:
            maxid = -1
        else:
            maxid = self.data.id.max().values[()]

        if id is None:
            id = np.arange(start=maxid+1,stop=maxid+1+nbodies,dtype=int)

        if name is None:
            name=np.char.mod(f"Body%d",id)

        if len(self.data) > 0:
            dup_id = np.in1d(id, self.data.id)
            if any(dup_id):
                raise ValueError(f"Duplicate ids detected: ", *id[dup_id])

        time = [self.param['TSTART']]

        if mass is not None:
            if Gmass is not None:
                raise ValueError("Cannot use mass and Gmass inputs simultaneously!")
            else: 
                Gmass = self.GU * mass

        dsnew = init_cond.vec2xr(self.param, name=name, a=a, e=e, inc=inc, capom=capom, omega=omega, capm=capm, id=id,
                                 Gmass=Gmass, radius=radius, rhill=rhill, Ip=Ip, rh=rh, vh=vh,rot=rot, j2rp2=J2, j4rp4=J4, time=time)

        dsnew = self._combine_and_fix_dsnew(dsnew)
        self.save(verbose=False)
        self.init_cond = self.data.copy(deep=True)

        return

    def _combine_and_fix_dsnew(self,dsnew):
        """
        Combines the new Dataset with the old one. Also computes the values of ntp and npl and sets the proper types.
        Parameters
        ----------
        dsnew : xarray Dataset
            Dataset with new bodies

        Returns
        -------
        dsnew : xarray Dataset
            Updated Dataset with ntp, npl values and types fixed.

        """
        if "id" not in self.data.dims:
            if len(np.unique(dsnew['name'])) == len(dsnew['name']):
               dsnew = dsnew.swap_dims({"id" : "name"})
               if "id" in dsnew:
                   dsnew = dsnew.reset_coords("id")
            else:
                msg = "Non-unique names detected for bodies. The Dataset will be dimensioned by integer id instead of name."
                msg +="\nConsider using unique names instead."
                print(msg)

        self.data = xr.combine_by_coords([self.data, dsnew])

        if self.param['OUT_TYPE'] == "NETCDF_DOUBLE":
            dsnew = io.fix_types(dsnew, ftype=np.float64)
            self.data = io.fix_types(self.data, ftype=np.float64)
        elif self.param['OUT_TYPE'] == "NETCDF_FLOAT":
            dsnew = io.fix_types(dsnew, ftype=np.float32)
            self.data = io.fix_types(self.data, ftype=np.float32)

        def get_nvals(ds):
            if "name" in ds.dims:
                count_dim = "name"
            elif "id" in ds.dims:
                count_dim = "id"
            if "Gmass" in ds:
                ds['ntp'] = ds[count_dim].where(np.isnan(ds['Gmass'])).count(dim=count_dim)
                ds['npl'] = ds[count_dim].where(~(np.isnan(ds['Gmass']))).count(dim=count_dim) - 1
                if self.integrator == "symba" and "GMTINY" in self.param and self.param['GMTINY'] is not None:
                    ds['nplm'] = ds[count_dim].where(ds['Gmass'] > self.param['GMTINY']).count(dim=count_dim) - 1
            else:
                ds['ntp'] = ds[count_dim].count(dim=count_dim)
                ds['npl'] = xr.full_like(ds['ntp'],0)
                if self.integrator == "symba" and "GMTINY" in self.param and self.param['GMTINY'] is not None:
                   ds['nplm'] = xr.full_like(ds['ntp'],0)
            return ds

        dsnew = get_nvals(dsnew)
        self.data = get_nvals(self.data)

        self.data = self.data.sortby("id")
        self.data = io.reorder_dims(self.data)

        return dsnew

    def read_param(self,
                   param_file : os.PathLike | str | None = None,
                   codename: Literal["Swiftest", "Swifter", "Swift"] | None = None,
                   read_init_cond : bool | None = None,
                   verbose: bool | None = None,
                   dask: bool = False, 
                   **kwargs : Any):
        """
        Reads in an input parameter file and stores the values in the param dictionary.
        
        Parameters
        ----------
        param_file : str or path-like, default is the value of the Simulation object's internal `param_file`.
            File name of the input parameter file
        codename : {"Swiftest", "Swifter", "Swift"}, default is the value of the Simulation object's internal`codename`
            Type of parameter file, either "Swift", "Swifter", or "Swiftest"
        read_init_cond : bool, optional
            If true, will read in the initial conditions file into the data instance variable. Default True
        verbose : bool, default is the value of the Simulation object's internal `verbose`
            If set to True, then more information is printed by Simulation methods as they are executed. Setting to
            False suppresses most messages other than errors.
        dask : bool, default False
            Use Dask to lazily load data (useful for very large datasets)
            
        Returns
        -------
        True if the parameter file exists and is read correctly. False otherwise.
        """
        if param_file is None:
            param_file = self.simdir / self.param_file
        if read_init_cond is None:
            read_init_cond = True
        if codename is None:
            codename = self.codename

        if verbose is None:
            verbose = self.verbose

        if not os.path.exists(param_file):
            return False

        if codename == "Swiftest":
            self.param = io.read_swiftest_param(param_file, self.param, verbose=verbose)
            if read_init_cond:
                if "NETCDF" in self.param['IN_TYPE']:
                   init_cond_file = self.simdir / self.param['NC_IN']
                   if os.path.exists(init_cond_file):
                       param_tmp = self.param.copy()
                       param_tmp['BIN_OUT'] = init_cond_file
                       self.data = io.swiftest2xr(param_tmp, verbose=self.verbose, dask=dask)
                       self.init_cond = self.data.copy(deep=True)
                   else:
                       warnings.warn(f"Initial conditions file file {init_cond_file} not found.", stacklevel=2)
                else:
                    warnings.warn("Reading in ASCII initial conditions files in Python is not yet supported")
        elif codename == "Swifter":
            self.param = io.read_swifter_param(param_file, verbose=verbose)
        elif codename == "Swift":
            self.param = io.read_swift_param(param_file, verbose=verbose)
        else:
            warnings.warn(f'{codename} is not a recognized code name. Valid options are "Swiftest", "Swifter", or "Swift".',stacklevel=2)
            return False

        return True

    def write_param(self,
                    codename: Literal["Swiftest", "Swifter", "Swift"]  | None = None,
                    param_file: str | os.PathLike | None = None,
                    param: Dict | None = None,
                    **kwargs: Any):
        """
        Writes to a param.in file and determines whether the output format needs to be converted between Swift/Swifter/Swiftest.
        
        Parameters
        ----------
        codename : {"Swiftest", "Swifter", "Swift"}, optional
            Alternative name of the n-body code that the parameter file will be formatted for. Defaults to current instance
            variable codename
        param_file : str or path-like, optional
            Alternative file name of the input parameter file. Defaults to current instance variable self.param_file
        param: Dict, optional
            An alternative parameter dictionary to write out. Defaults to the current instance variable self.param
        **kwargs
            A dictionary of additional keyword argument. These are ignored.

        Returns
        -------
        None
        """

        if codename is None:
            codename = self.codename
        if param_file is None:
            param_file = self.param_file
        if param is None:
            param = self.param

        if "verbose" in kwargs:
            verbose = kwargs['verbose']
        else:
            verbose = self.verbose

        if verbose:
            print(f"Writing parameter inputs to file {str(self.simdir / param_file)}")
        param['! VERSION'] = f"{codename} input file"
        
        self.simdir.mkdir(parents=True, exist_ok=True)
        # Check to see if the parameter type matches the output type. If not, we need to convert
        if codename == "Swifter" or codename == "Swiftest":
            if param['IN_TYPE'] == "ASCII":
                param.pop("NC_IN", None)
            else:
                param.pop("CB_IN", None)
                param.pop("PL_IN", None)
                param.pop("TP_IN", None)
            io.write_labeled_param(param, self.simdir / param_file)
        elif codename == "Swift":
            io.write_swift_param(param, self.simdir / param_file)
        else:
            warnings.warn('Cannot process unknown code type. Call the read_param method with a valid code name. Valid options are "Swiftest", "Swifter", or "Swift".',stacklevel=2)
            
        return

    def convert(self, param_file, newcodename="Swiftest", plname="pl.swiftest.in", tpname="tp.swiftest.in",
                cbname="cb.swiftest.in", conversion_questions={}, dask = False):
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
        dask : bool, default False
            Use Dask to lazily load data (useful for very large datasets)

        Returns
        -------
        oldparam : xarray dataset
        The old parameter configuration.
        """
        oldparam = self.param
        if self.codename == newcodename:
            warnings.warn(f"This parameter configuration is already in {newcodename} format",stacklevel=2)
            return oldparam
        if newcodename != "Swift" and newcodename != "Swifter" and newcodename != "Swiftest":
            warnings.warn(f'{newcodename} is an invalid code type. Valid options are "Swiftest", "Swifter", or "Swift".',stacklevel=2)
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
            warnings.warn(f"Conversion from {self.codename} to {newcodename} is not supported.",stacklevel=2)
        return oldparam

    def read_output_file(self,read_init_cond : bool = True, dask : bool = False):
        """
        Reads in simulation data from an output file and stores it as an Xarray Dataset in the `data` instance variable.

        Parameters
        ----------
        read_init_cond : bool, default True
            Read in an initial conditions file along with the output file. Default is True
        dask : bool, default False
            Use Dask to lazily load data (useful for very large datasets)
        Returns
        -------
        self.data : xarray dataset
        """

        # Make a temporary copy of the parameter dictionary so we can supply the absolute path of the binary file
        # This is done to handle cases where the method is called from a different working directory than the simulation
        # results

        param_tmp = self.param.copy()
        param_tmp['BIN_OUT'] = os.path.join(self.simdir, self.param['BIN_OUT'])
        if self.codename == "Swiftest":
            self.data = io.swiftest2xr(param_tmp, verbose=self.verbose, dask=dask)
            if self.verbose:
                print('Swiftest simulation data stored as xarray DataSet .data')
            if read_init_cond:
                if self.verbose:
                    print("Reading initial conditions file as .init_cond")
                if "NETCDF" in self.param['IN_TYPE']:
                    param_tmp['BIN_OUT'] = self.simdir / self.param['NC_IN']
                    self.init_cond = io.swiftest2xr(param_tmp, verbose=False, dask=dask)
                else:
                    self.init_cond = self.data.isel(time=0)

            if self.read_encounters:
                self.read_encounter_file(dask=dask)
            if self.read_collisions:
                self.read_collision_file(dask=dask)
            if self.verbose:
                print("Finished reading Swiftest dataset files.")

        elif self.codename == "Swifter":
            self.data = io.swifter2xr(param_tmp, verbose=self.verbose)
            if self.verbose: print('Swifter simulation data stored as xarray DataSet .data')
        elif self.codename == "Swift":
            warnings.warn("Reading Swift simulation data is not implemented yet",stacklevel=2)
        else:
            warnings.warn('Cannot process unknown code type. Call the read_param method with a valid code name. Valid options are "Swiftest", "Swifter", or "Swift".',stacklevel=2)
        return

    def read_encounter_file(self, dask=False):
        enc_file = self.simdir / "encounters.nc"
        if not os.path.exists(enc_file):
           return

        if dask:
            self.encounters = xr.open_mfdataset(enc_file,engine='h5netcdf', mask_and_scale=False)
        else:
            self.encounters = xr.open_dataset(enc_file, mask_and_scale=False)
        self.encounters = io.process_netcdf_input(self.encounters, self.param)

        # Remove any overlapping time values
        tgood,tid = np.unique(self.encounters.time,return_index=True)
        self.encounters = self.encounters.isel(time=tid)
        # Remove any NaN values
        tgood=self.encounters.time.where(~np.isnan(self.encounters.time),drop=True)
        self.encounters = self.encounters.sel(time=tgood)

        return

    def read_collision_file(self, dask=False):
       
        col_file = self.simdir / "collisions.nc"
        if not os.path.exists(col_file):
           return

        if self.verbose:
                print("Reading collisions history file as .collisions")

        if dask:
            self.collisions = xr.open_mfdataset(col_file,engine='h5netcdf', mask_and_scale=False)
        else:
            self.collisions = xr.open_dataset(col_file, mask_and_scale=False)
            
        self.collisions = io.process_netcdf_input(self.collisions, self.param)

        return

    def follow(self, codestyle="Swifter", dask=False):
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
        if self.data is None:
            self.read_output_file(dask=dask)
        if codestyle == "Swift":
            try:
                with open('follow.in', 'r') as f:
                    line = f.readline()  # Parameter file (ignored because read_output_file already takes care of it
                    line = f.readline()  # PL file (ignored)
                    line = f.readline()  # TP file (ignored)
                    line = f.readline()  # ifol
                    i_list = [i for i in line.split(" ") if i.strip()]
                    ifol = int(i_list[0])
                    line = f.readline()  # nskp
                    i_list = [i for i in line.split(" ") if i.strip()]
                    nskp = int(i_list[0])
            except IOError:
                warnings.warn('No follow.in file found',stacklevel=2)
                ifol = None
                nskp = None
            fol = tool.follow_swift(self.data, ifol=ifol, nskp=nskp)
        else:
            fol = None

        if self.verbose: print('follow.out written')
        return fol

    def save(self,
             codename: Literal["Swiftest", "Swifter", "Swift"] | None = None,
             param_file: str | os.PathLike | None = None,
             param: Dict | None = None,
             framenum: int = -1,
             **kwargs: Any):
        """
        Saves an xarray dataset to a set of input files.

        Parameters
        ----------
        codename : {"Swiftest", "Swifter", "Swift"}, optional
            Alternative name of the n-body code that the parameter file will be formatted for. Defaults to current instance
            variable self.codename
        param_file : str or path-like, optional
            Alternative file name of the input parameter file. Defaults to current instance variable self.param_file
        param: Dict, optional
            An alternative parameter dictionary to write out. Defaults to the current instance variable self.param
        framenum : int Default=-1
            Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.
        **kwargs
            A dictionary of additional keyword argument. These are ignored.

        Returns
        -------
        None
        """

        if "verbose" in kwargs:
            verbose = kwargs['verbose']
        else:
            verbose = self.verbose
        
        if codename is None:
            codename = self.codename
        if param_file is None:
            param_file = self.param_file
        if param is None:
            param = self.param

        if not self.simdir.exists():
            self.simdir.mkdir(parents=True, exist_ok=True)
        
        if codename == "Swiftest":
            infile_name = Path(self.simdir) / param['NC_IN']
            io.swiftest_xr2infile(ds=self.data, param=param, in_type=self.param['IN_TYPE'], infile_name=infile_name, framenum=framenum, verbose=verbose)
            self.write_param(param_file=param_file,**kwargs)
        elif codename == "Swifter":
            swifter_param = io.swiftest2swifter_param(param)
            if "rhill" in self.data:
                swifter_param['RHILL_PRESENT'] = 'YES'
            io.swifter_xr2infile(self.data, swifter_param, self.simdir, framenum)
            self.write_param(codename=codename,param_file=param_file,param=swifter_param,**kwargs)
        else:
            warnings.warn(f'Saving to {codename} not supported',stacklevel=2)

        return

    def initial_conditions_from_bin(self, framenum=-1, new_param=None, new_param_file="param.new.in",
                                    new_initial_conditions_file="bin_in.nc", restart=False, codename="Swiftest"):
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
            A dataset containing the extracted initial condition data.
        """

        frame = None
        if codename != "Swiftest":
            self.save(new_param_file, framenum, codename)
            return

        if new_param is None:
            new_param = self.param.copy()

        if codename == "Swiftest":
            if restart:
                new_param['T0'] = self.data.time.values[framenum]
            if self.param['OUT_TYPE'] == 'NETCDF_DOUBLE':
                new_param['IN_TYPE'] = 'NETCDF_DOUBLE'
            elif self.param['OUT_TYPE'] == 'NETCDF_FLOAT':
                new_param['IN_TYPE'] = 'NETCDF_FLOAT'
            else:
                warnings.warn(f"{self.param['OUT_TYPE']} is an invalid OUT_TYPE file",stacklevel=2)
                return

            if self.param['BIN_OUT'] != new_param['BIN_OUT'] and restart:
                print(f"Restart run with new output file. Copying {self.param['BIN_OUT']} to {new_param['BIN_OUT']}")
                shutil.copy2(self.param['BIN_OUT'], new_param['BIN_OUT'])

            new_param['IN_FORM'] = 'XV'
            if restart:
                new_param['OUT_STAT'] = 'APPEND'

            new_param['FIRSTKICK'] = 'T'
            new_param['NC_IN'] = new_initial_conditions_file
            new_param.pop('PL_IN', None)
            new_param.pop('TP_IN', None)
            new_param.pop('CB_IN', None)
            print(f"Extracting data from dataset at time frame number {framenum} and saving it to {new_param['NC_IN']}")
            frame = io.swiftest_xr2infile(self.data, self.param, infile_name=new_param['NC_IN'], framenum=framenum)
            print(f"Saving parameter configuration file to {new_param_file}")
            self.write_param(new_param_file, param=new_param)

        return frame

    def clean(self):
        """
        Cleans up simulation directory by deleting old files (dump, logs, and output files).
        """
        old_files = [self.simdir / self.param['BIN_OUT'],
                     self.simdir / "fraggle.log",
                     self.simdir / "swiftest.log",
                     self.simdir / "collisions.log",
                     self.simdir / "collisions.nc",
                     self.simdir / "encounters.nc",
                     ]
        
        glob_files = [self.simdir.glob("**/param.*.in")]

        for f in old_files:
            if f.exists():
                os.remove(f)
                
        for g in glob_files:
            for f in g:
               if f.exists():
                  os.remove(f)        
        return

