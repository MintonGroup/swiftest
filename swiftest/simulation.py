"""
Copyright 2025 - David Minton.

This file is part of Swiftest.
Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Swiftest.
If not, see: https://www.gnu.org/licenses.
"""

from __future__ import annotations

import contextlib
import datetime
import os
import shutil
import warnings
from pathlib import Path
from typing import Any, Literal, Union

import numpy as np
import numpy.typing as npt
import xarray as xr

from . import constants, init_cond, io, tool
from .data import SwiftestDataset

FloatLike = float | int | np.number


@contextlib.contextmanager
def _cwd(newdir: Path):
    olddir = Path.cwd()
    os.chdir(str(newdir))
    try:
        yield
    finally:
        os.chdir(str(olddir))


class Simulation:
    """
    Defines the basic Swift/Swifter/Swiftest simulation object.

    """

    def __init__(
        self,
        simdir: os.PathLike | Path | str = "simdata",
        param_file: os.PathLike | Path | str = "param.in",
        read_param: bool = False,
        read_data: bool = False,
        read_init_cond: bool | None = None,
        read_collisions: bool | None = None,
        read_encounters: bool | None = None,
        clean: bool = False,
        dask: bool = False,
        codename: Literal["Swiftest", "Swifter", "Swift"] = "Swiftest",
        integrator: Literal["symba", "rmvs", "whm", "helio"] = "symba",
        coarray: bool = False,
        verbose: bool = True,
        **kwargs: Any,
    ):
        """
        Set up a new simulation object with the given parameters.

        Parameters for a given Simulation object can be set a number of different ways, including via a parameter input
        file, arguments to Simulation, the general `set_parameter` method, or the specific setters for groups of
        similar parameters (e.g. set_init_cond_files, set_simulation_time, etc.). Each parameter has a default value
        that can be overridden by an argument to Simulation(). Some argument parameters have equivalent values that
        are passed to the Swiftest driver Fortran function via a parameter input file. When declaring a new
        Simulation object, parameters are chosen in the following way, from highest to lowest priority

        #. Arguments to Simulation()
        #. The parameter input file given by `param_file` under the following conditions

            - `read_param` is set to True (default behavior).
            - The file given by `param_file` exists. The default file is `param.in` located in the `simdata` directory
              inside the current working directory, which can be changed by passing `param_file` as an argument.
            - The argument has an equivalent parameter or set of parameters in the parameter input file.

        #. Default values (see below)

        Parameters
        ----------
        simdir : str, path-like, or None, default "simdata"
            The directory where the simulation data will be stored.
        param_file : str, path-like, or None, default "param.in"
            The name of the parameter input file that will be passed to the integrator. The parameter file should be
            located in the directory set by `simdir`. If the file does not exist, it will be created.
        read_param : bool, default True
            If true, read in a pre-existing parameter input file given by the argument `param_file` if it exists.
            Otherwise, create a new parameter file using the arguments passed to Simulation or defaults
        read_data : bool, default False
            If True, read in a pre-existing binary input file given by the argument `output_file_name` if it exists.
            If True, then `read_param` will also be set to True.
        read_init_cond: bool, default None
            If True, read in a pre-existing initial condition file given by the param.in parameter `NC_IN` or the argument
            `init_cond_file_name` if it exists.
            If None, then it will take the value of `read_data`.
            If True, then `read_param` will also be set to True
        read_collisions : bool, default None
            If True, read in a pre-existing collision file `collisions.nc`.
            If None, then it will take the value of `read_data`.
            If True, then `read_param` will also be set to True
        read_encounters : bool, default None
            If True, read in a pre-existing encounter file `encounters.nc`.
            If True, then `read_param` will also be set to True
            If None, then it will take the value of `read_data`.
        clean : bool, default False
            If True, then delete the `simdir` directory before initializing the simulation. If `read_param` is True, then this will
            be ignored.
        dask : bool, default False
            If true, will use Dask to lazily load data (useful for very large datasets).
        coarray : bool, default False
            If true, will employ Coarrays on test particle structures to run in single program/multiple data parallel mode.
            In order to use this capability, Swiftest must be compiled for Coarray support. Only certain integrators can use
            Coarrays. RMVS, WHM, Helio are all compatible, but SyMBA is not, due to the way tp-pl close encounters are handeled.
        verbose : bool, default False
            If set to True, then more information is printed by Simulation methods as they are executed. Setting to
            False suppresses most messages other than errors and some warnings.
        **kwargs : Any
            Any valid keyword arguments accepted by :meth:`~swiftest.Simulation.set_parameter` (see below)
        """
        self.verbose = verbose

        # Define some instance variables
        self._getter_column_width = 32
        self._param = {}
        self._param_file = None
        self._data = SwiftestDataset()
        self._init_cond = SwiftestDataset()
        self._encounters = SwiftestDataset()
        self._collisions = SwiftestDataset()
        self._MU_name = "MU"
        self._DU_name = "DU"
        self._TU_name = "TU"
        self._MU2KG = None
        self._KG2MU = None
        self._TU2S = None
        self._S2TU = None
        self._DU2M = None
        self._M2DU = None
        self._GU = None
        self._ephemeris_date = constants.MINTON_BCL
        self._codename = None
        self._integrator = None
        self._restart = False

        # Define a dictionary that maps parameter names used in `param.in` files to arguments passed to Simulation() and `set_parameter()`
        self._param_to_argument = {
            "T0": "t0",
            "TSTART": "tstart",
            "TSTOP": "tstop",
            "DT": "dt",
            "ISTEP_OUT": "istep_out",
            "NSTEP_OUT": "nstep_out",
            "DUMP_CADENCE": "dump_cadence",
            "GMTINY": "gmtiny",
            "CHK_CLOSE": "close_encounter_check",
            "COLLISION_MODEL": "collision_model",
            "ENCOUNTER_SAVE": "encounter_save",
            "MIN_GMFRAG": "minimum_fragment_gmass",
            "NFRAG_REDUCTION": "nfrag_reduction",
            "ROTATION": "rotation",
            "GR": "general_relativity",
            "ENERGY": "compute_conservation_values",
            "RHILL_PRESENT": "rhill_present",
            "EXTRA_FORCE": "extra_force",
            "BIG_DISCARD": "big_discard",
            "INTERACTION_LOOPS": "interaction_loops",
            "ENCOUNTER_CHECK": "encounter_check_loops",
            "COARRAY": "coarray",
            "RESTART": "restart",
            "IN_TYPE": "init_cond_file_type",
            "IN_FORM": "init_cond_format",
            "NC_IN": "init_cond_file_name",
            "CB_IN": "init_cond_file_name['CB']",
            "PL_IN": "init_cond_file_name['PL']",
            "TP_IN": "init_cond_file_name['TP']",
            "OUT_TYPE": "output_file_type",
            "BIN_OUT": "output_file_name",
            "OUT_FORM": "output_format",
            "OUT_STAT": "restart",
            "MU2KG": "MU2KG",
            "DU2M": "DU2M",
            "TU2S": "TU2S",
            "CHK_RMIN": "rmin",
            "CHK_RMAX": "rmax",
            "CHK_QMIN_COORD": "qmin_coord",
            "CHK_QMIN": "qmin",
            "CHK_QMIN_RANGE": "qminR",
            "SEED": "seed",
        }

        # Define default parameters
        self.MU_name = "Msun"
        self.MU2KG = constants.MSun
        self.DU_name = "AU"
        self.DU2M = constants.AU2M
        self.TU_name = "YR"
        self.TU2S = constants.YR2S
        self.param = {
            "T0": 0.0,
            "TSTART": 0.0,
            "ISTEP_OUT": 1,
            "DUMP_CADENCE": 1,
            "IN_TYPE": "NETCDF_DOUBLE",
            "NC_IN": "init_cond.nc",
            "IN_FORM": "XV",
            "OUT_TYPE": "NETCDF_DOUBLE",
            "BIN_OUT": "data.nc",
            "OUT_FORM": "XVEL",
            "OUT_STAT": "REPLACE",
            "MU2KG": self.MU2KG,
            "DU2M": self.DU2M,
            "TU2S": self.TU2S,
            "CHK_RMIN": constants.RSun * self.M2DU,
            "CHK_RMAX": 10000.0,
            "CHK_EJECT": 10000.0,
            "CHK_QMIN_RANGE": "-1 10000.0",
            "CHK_QMIN_COORD": "HELIO",
            "CHK_CLOSE": True,
            "GR": True,
            "GMTINY": 0.0,
            "COLLISION_MODEL": "FRAGGLE",
            "MIN_GMFRAG": 0.0,
            "NFRAG_REDUCTION": 30.0,
            "ROTATION": True,
            "ENERGY": True,
            "EXTRA_FORCE": False,
            "BIG_DISCARD": False,
            "RHILL_PRESENT": False,
            "INTERACTION_LOOPS": "TRIANGULAR",
            "ENCOUNTER_CHECK": "TRIANGULAR",
            "RESTART": False,
            "ENCOUNTER_SAVE": "NONE",
            "TIDES": False,
        }

        self.codename = codename
        self.integrator = integrator
        self.simdir = simdir

        if read_init_cond is None:
            read_init_cond = read_data

        if read_collisions is None:
            read_collisions = read_data

        if read_encounters is None:
            read_encounters = read_data

        if read_init_cond or read_data or read_collisions or read_encounters:
            read_param = True

        if clean and not read_param:
            self.clean(deep=True)

        self.param_file = param_file

        # Parameters are set in reverse priority order. First the defaults, then values from a pre-existing input file,
        # then using the arguments passed via **kwargs.
        # --------------------------
        # Lowest Priority: Defaults:
        # --------------------------
        # Set all parameters to their defaults.
        self.set_parameter()

        # -----------------------------------------------------------------
        # Higher Priority: Values from a file (if requested and it exists)
        # -----------------------------------------------------------------

        # If the user asks to read in an old parameter file or output file, override any default parameters with values from the file
        # If the file doesn't exist, flag it for now so we know to create it

        if read_param:
            self.read_param(**kwargs)

        # -----------------------------------------------------------------
        # Highest Priority: Values from arguments passed to Simulation()
        # -----------------------------------------------------------------
        if len(kwargs) > 0:
            self.set_parameter(**kwargs)

        # Read in an old simulation file if requested
        if read_init_cond:
            icpath = self.simdir / self.param["NC_IN"]
            if icpath.exists():
                self.read_init_cond(dask=dask)

        if read_data:
            binpath = self.simdir / self.param["BIN_OUT"]
            if binpath.exists():
                self.read_output_file(
                    dask=dask,
                    read_collisions=read_collisions,
                    read_encounters=read_encounters,
                )
        else:
            if read_encounters:
                self.read_encounter_file(dask=dask, verbose=verbose)
            if read_collisions:
                self.read_collision_file(dask=dask, verbose=verbose)
        return

    def _run_swiftest_driver(self, **kwargs: Any) -> None:
        """
        Internal callable function that executes the swiftest_driver run.
        """
        from .core import driver

        verbose = kwargs.pop("verbose", self.verbose)

        if "display_style" in kwargs:
            display_style = kwargs.pop("display_style")
        else:
            if verbose:
                display_style = "progress"
            else:
                display_style = "quiet"

        # Set up signal handling to catch termination signals
        try:
            with _cwd(self.simdir):
                driver(self.integrator, str(self.param_file), display_style)
        except Exception as e:
            warnings.warn(f"Failed to run the simulation with Exception: {e}", stacklevel=2)

        return

    def run(self, dask: bool = False, **kwargs: Any) -> None:
        """
        Runs a Swiftest integration.

        Uses the parameters set by the `param` dictionary unless overridden by keyword arguments. Accepts any keyword arguments that can be passed to `set_parameter`.

        Parameters
        ----------
        dask : bool, default False
            If true, will use Dask to lazily load data (useful for very large datasets)
        **kwargs : Any
            Any valid keyword arguments accepted by :meth:`~swiftest.Simulation.set_parameter`

        Returns
        -------
        None
        """
        verbose = kwargs.pop("verbose", self.verbose)

        if len(kwargs) > 0:
            old_verbose = self.verbose
            self.set_parameter(verbose=verbose, **kwargs)
            self.verbose = old_verbose

        if self.codename != "Swiftest":
            if verbose:
                warnings.warn(
                    f"Running an integration is not yet supported for {self.codename}",
                    stacklevel=2,
                )
            return

        # Save initial conditions
        if not self.restart:
            self.clean(verbose=verbose)
            self.save(framenum=0, verbose=verbose)

        # Write out the current parameter set before executing run
        self.write_param(verbose=verbose, **kwargs)

        if self.param["DT"] > (self.param["TSTOP"] - self.param["TSTART"]) and verbose:
            msg = "dt should be smaller than tstop-tstart"
            warnings.warn(msg, stacklevel=2)

        if verbose:
            print(
                f"Running a {self.codename} {self.integrator} run from tstart={self.param['TSTART']} {self.TU_name} to tstop={self.param['TSTOP']} {self.TU_name}"
            )

        # This is temporary until a better algorithm for setting OMP_NUM_THREADS is implemented
        omp_num_threads_old = None
        if self.init_cond.isel(time=0).npl.values[()] < 100:
            if "OMP_NUM_THREADS" in os.environ:
                omp_num_threads_old = os.environ["OMP_NUM_THREADS"]
            os.environ["OMP_NUM_THREADS"] = "1"

        # Prevent NetCDF file locking
        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

        self._run_swiftest_driver(verbose=verbose)
        if verbose:
            print("\nRun complete.")
        if omp_num_threads_old:
            os.environ["OMP_NUM_THREADS"] = omp_num_threads_old
        else:
            del os.environ["OMP_NUM_THREADS"]

        # Read in new data
        self.read_output_file(dask=dask, read_encounters=True, read_collisions=True, verbose=verbose)

        return

    def _get_valid_arg_list(self, arg_list: str | list[str] | None = None, valid_var: dict | None = None) -> tuple[list[str], dict]:
        """
        Internal function for getters that extracts subset of arguments that is contained in the dictionary of valid argument/parameter variable pairs.

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

    def _create_valid_var(self, valid_arg):
        """
        Internal function for getters that extracts subset of arguments that is contained in the dictionary of valid argument/parameter variable pairs.

        Parameters
        ----------
        valid_arg : List
            A list of argument values

        Returns
        -------
        Dict :
            A dictionary where each key is the argument name and the values are the corresponding parameter name
        """
        valid_var = {arg: key for key, arg in self._param_to_argument.items() if arg in valid_arg}

        return valid_var

    def set_simulation_time(
        self,
        t0: float | None = None,
        tstart: float | None = None,
        tstop: float | None = None,
        dt: float | None = None,
        istep_out: int | None = None,
        tstep_out: float | None = None,
        nstep_out: int | None = None,
        dump_cadence: int | None = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """
        Set the parameters that control how a simulation is run, such as start and stop time, step size, and the cadence of output to both the screen and to file. Returns a dictionary of the parameters that were set.

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
            The number of time steps between output saves to file. Only `istep_out` or `tstep_out` can be set.
            Parameter input file equivalent is `ISTEP_OUT`
        tstep_out : float, optional
            The approximate time between when outputs are written to file. Passing this computes::

                `istep_out = floor(tstep_out/dt)`.

            Only `istep_out` or `tstep_out` can be set.
            Parameter input file equivalent is None
        nstep_out : int, optional
            The total number of times that outputs are written to file. Passing this allows for a geometric progression of output steps:

                TSTART, f**0 * TSTEP_OUT, f**1 * TSTEP_OUT, f**2 * TSTEP_OUT, ..., f**(nstep_out-1) * TSTEP_OUT

            where `f` is a factor that can stretch (or shrink) the time between outputs. Setting::

                nstep_out = int((tstart - tstop) / (tstep_out))

            is equivalent to the standard linear output (i.e. `f==1`) and is the same as not passing anything for this argument.
            Passing `nstep_out` requires passing either `istep_out` or `tstep_out` as well.
        dump_cadence : int, optional
            The number of output steps (given by `istep_out`) between when the saved data is dumped to a file. Setting it to 0
            is equivalent to only dumping data to file at the end of the simulation. Default value is 10.
            Parameter input file equivalent is `DUMP_CADENCE`
        **kwargs
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        time_dict : dict
           A dictionary containing the requested parameters
        """
        if (
            t0 is None
            and tstart is None
            and tstop is None
            and dt is None
            and istep_out is None
            and tstep_out is None
            and dump_cadence is None
        ):
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

        self.param["T0"] = t0
        self.param["TSTART"] = tstart

        if tstop is None:
            tstop = self.param.pop("TSTOP", None)
        else:
            update_list.append("tstop")

        if tstop is not None and tstop < tstart:
            warnings.warn("tstop should be greater than tstart.", stacklevel=2)

        if tstop is not None:
            self.param["TSTOP"] = tstop

        if dt is None:
            dt = self.param.pop("DT", None)
        else:
            update_list.append("dt")

        if dt is not None:
            self.param["DT"] = dt

        if istep_out is None and tstep_out is None:
            istep_out = self.param.pop("ISTEP_OUT", None)
        elif istep_out is not None and tstep_out is not None:
            warnings.warn("istep_out and tstep_out cannot both be set", stacklevel=2)
            return {}
        else:
            update_list.append("istep_out")

        if tstep_out is not None and tstep_out > tstop:
            warnings.warn(
                "tstep_out must be less than tstop. Setting tstep_out=tstop",
                stacklevel=2,
            )
            tstep_out = tstop

        if tstep_out is not None and dt is not None:
            istep_out = int(tstep_out / dt)

        if istep_out is not None:
            self.param["ISTEP_OUT"] = int(istep_out)

        if nstep_out is not None:
            if istep_out is None:
                warnings.warn(
                    "nstep_out requires either istep_out or tstep_out to also be set",
                    stacklevel=2,
                )
            else:
                self.param["NSTEP_OUT"] = int(nstep_out)

        if dump_cadence is None:
            dump_cadence = self.param.pop("DUMP_CADENCE", 1)
        else:
            update_list.append("dump_cadence")
        self.param["DUMP_CADENCE"] = dump_cadence

        time_dict = self.get_simulation_time(update_list)

        return time_dict

    def get_simulation_time(self, arg_list: str | list[str] | None = None, **kwargs: Any) -> dict[str, Any]:
        """
        Returns a subset of the parameter dictionary containing the current simulation time parameters.

        Parameters
        ----------
        arg_list : str | List[str], optional
            A single string or list of strings containing the names of the simulation time parameters to extract.
            Default is all of:
            ["t0", "tstart", "tstop", "dt", "istep_out", "tstep_out", "dump_cadence"]
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        time_dict : dict
           A dictionary containing the requested parameters
        """
        valid_arg = [
            "t0",
            "tstart",
            "tstop",
            "dt",
            "istep_out",
            "tstep_out",
            "nstep_out",
            "dump_cadence",
        ]
        valid_var = self._create_valid_var(valid_arg)
        verbose = kwargs.pop("verbose", self.verbose)

        units = {
            "t0": self.TU_name,
            "tstart": self.TU_name,
            "tstop": self.TU_name,
            "dt": self.TU_name,
            "tstep_out": self.TU_name,
            "istep_out": "",
            "nstep_out": "",
            "dump_cadence": "",
        }

        tstep_out = None
        if (
            arg_list is None
            or "tstep_out" in arg_list
            or "istep_out" in arg_list
            and ("ISTEP_OUT" in self.param and "DT" in self.param)
        ):
            istep_out = self.param["ISTEP_OUT"]
            dt = self.param["DT"]
            tstep_out = istep_out * dt

        valid_arg, time_dict = self._get_valid_arg_list(arg_list, valid_var)

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

    def set_parameter(self, **kwargs: Any) -> dict[str, Any]:
        """
        Setter for all possible parameters. This will call each of the specialized setters using keyword arguments.

        If no arguments are passed, then either values from the param_file or default values will be used.

        Parameters
        ----------
        **kwargs : dict
            See list of valid parameters and their defaults below
        param_file : str, path-like, or file-like, default None
            Name of the parameter input file that will be passed to the integrator.
            If passed, then the file will be read in and the parameters set from them, but can be overriden by any explicit arguments
        codename : {"Swiftest", "Swifter", "Swift"}, default "Swiftest"
            Name of the n-body code that will be used.
            Parameter input file equivalent is None
        integrator : {"symba","rmvs","whm","helio"}, default "symba"
            Name of the n-body integrator that will be used when executing a run.
            Parameter input file equivalent is None
        t0 : float, default 0.0
            The reference time for the start of the simulation. Defaults is 0.0.
            Parameter input file equivalent is `T0`
        tstart : float, default 0.0
            The start time for a restarted simulation. For a new simulation, tstart will be set to t0 automatically.
            Parameter input file equivalent is `TSTART`
        tstop : float, optional
            The stopping time for a simulation. `tstop` must be greater than `tstart`.
            Parameter input file equivalent is `TSTOP`
        dt : float, optional
            The step size of the simulation. `dt` must be less than or equal to `tstop-tstart`.
            Parameter input file equivalent is `DT`
        istep_out : int, optional
            The number of time steps between output saves to file. Only `istep_out` or `tstep_out` can be set.
            Parameter input file equivalent is `ISTEP_OUT`
        tstep_out : float, optional
            The approximate time between when outputs are written to file.
            Passing this computes::

                istep_out = floor(tstep_out / dt)

            Only `istep_out` or `tstep_out` can be set. `tstep_out` must be less than `tstop`.
            Parameter input file equivalent is None
        nstep_out : int, optional
            The total number of times that outputs are written to file. Passing this allows for a geometric progression of output
            steps::

                TSTART, f**0 * TSTEP_OUT, f**1 * TSTEP_OUT, f**2 * TSTEP_OUT, ..., f ** (nstep_out - 1) * TSTEP_OUT

            where `f` is a factor that can stretch (or shrink) the time between outputs.  Setting::

                nstep_out = int((tstart - tstop) / (tstep_out))

            is equivalent to the standard linear output (i.e. `f==1`) and is the same as not passing anything for this argument.
            Passing `nstep_out` requires passing either `istep_out` or `tstep_out` as well.
        dump_cadence : int, optional
            The number of output steps (given by `istep_out`) between when the saved data is dumped to a file. Setting it to 0
            is equivalent to only dumping data to file at the end of the simulation. Default value is 10.
            Parameter input file equivalent is `DUMP_CADENCE`
        init_cond_file_type : {"NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"}, default "NETCDF_DOUBLE"
            The file type containing initial conditions for the simulation.

            * "NETCDF_DOUBLE" A single initial conditions input file in NetCDF file format of type NETCDF_DOUBLE.
            * "NETCDF_FLOAT" A single initial conditions input file in NetCDF file format of type NETCDF_FLOAT.
            * "ASCII"  Three initial conditions files in ASCII format. The individual files define the central body,

            massive body, and test particle initial conditions.
            Parameter input file equivalent is `IN_TYPE`
        init_cond_file_name : str, path-like, or dict, optional
            Name of the input initial condition file or files. Whether to pass a single file name or a dictionary
            depends on the argument passed to `init_cond_file_type`. If `init_cond_file_type={"NETCDF_DOUBLE","NETCDF_FLOAT"}`,
            then this will be a single file name. If `init_cond_file_type="ASCII"` then this must be a dictionary where::

                init_cond_file_name = {
                    "CB" - [path to central body initial conditions file] (Swiftest only),
                    "PL" - [path to massive body initial conditions file],
                    "TP" - [path to test particle initial conditions file] }

            If no file name is provided then the following default file names will be used.

            * "NETCDF_DOUBLE", "NETCDF_FLOAT" `init_cond_file_name = "init_cond.nc"`
            * "ASCII" `init_cond_file_name = {"CB" : "cb.in", "PL" : "pl.in", "TP" : "tp.in"}`

            Parameter input file equivalent is `NC_IN`, `CB_IN`, `PL_IN`, `TP_IN`
        init_cond_format : {"EL", "XV"}, default "XV"
            Indicates whether the input initial conditions are given as orbital elements or cartesian position and
            velocity vectors.
            If `codename` is "Swift" or "Swifter", EL initial conditions are converted to XV.
            Parameter input file equivalent is `IN_FORM`
        output_file_type : {"NETCDF_DOUBLE", "NETCDF_FLOAT","REAL4","REAL8","XDR4","XDR8"}, default "NETCDF_DOUBLE"
            The file type for the outputs of the simulation. Compatible file types depend on the `codename` argument.

            * Swiftest - Only "NETCDF_DOUBLE" or "NETCDF_FLOAT" supported.
            * Swifter - Only "REAL4","REAL8","XDR4" or "XDR8"  supported.
            * Swift - Only "REAL4" supported.

            Parameter input file equivalent is `OUT_TYPE`
        output_file_name : str or path-like, optional
            Name of output file to generate. If not supplied, then one of the default file names are used, depending on
            the value passed to `output_file_type`. The default is "data.nc".
            Parameter input file equivalent is `BIN_OUT`
        output_format : {"XV","XVEL"}, default "XVEL"
            Specifies the format for the data saved to the output file. If "XV" then cartesian position and velocity
            vectors for all bodies are stored. If "XVEL" then the orbital elements are also stored.
            Parameter input file equivalent is `OUT_FORM`
        simdir : PathLike, optional
            Directory where simulation data will be stored, including the parameter file, initial conditions file, output file
        MU : str, optional
           The mass unit system to use. Case-insensitive valid options are

           * "Msun"   - Solar mass
           * "Mearth" - Earth mass
           * "kg"     - kilograms
           * "g","gm" - grams

        DU : str, optional
            The distance unit system to use. Case-insensitive valid options are

            * "AU"     - Astronomical Unit
            * "Rearth" - Earth radius
            * "km"     - kilometer
            * "m"      - meter
            * "cm"     - centimeter

        TU : str, optional
            The time unit system to use. Case-insensitive valid options are

            * "y","YR","year","years"      - Year
            * "d","day","days"             - Julian day
            * "s","sec","seconds","second" - second

            Parameter input file equivalent is None
        MU2KG : float, optional
            The conversion factor to multiply by the mass unit that would convert it to kilogram.
            Setting this overrides MU
            Parameter input file equivalent is `MU2KG`
        DU2M : float, optional
            The conversion factor to multiply by the distance unit that would convert it to meter.
            Setting this overrides DU
            Parameter input file equivalent is `DU2M`
        TU2S : float, optional
            The conversion factor to multiply by the time unit that would convert it to seconds.
            Setting this overrides TU
            Parameter input file equivalent is `TU2S`
        MU_name : str, optional
            The name of the mass unit. When setting one of the standard units via `MU` a name will be
            automatically set for the unit, so this argument will override the automatic name.
            Parameter input file equivalent is None
        DU_name : str, optional
            The name of the distance unit. When setting one of the standard units via `DU` a name will be
            automatically set for the unit, so this argument will override the automatic name.
            Parameter input file equivalent is None
        TU_name : str, optional
            The name of the time unit. When setting one of the standard units via `TU` a name will be
            automatically set for the unit, so this argument will override the automatic name.
            Parameter input file equivalent is None
        rmin : float, default value is the radius of the central body in the unit system defined by the unit input arguments.
            Minimum distance of the simulation
            Parameter input file equivalent are `CHK_QMIN`, `CHK_RMIN`, `CHK_QMIN_RANGE[0]`
        rmax : float, default value is 10000 AU in the unit system defined by the unit input arguments.
            Maximum distance of the simulation (CHK_RMAX, CHK_QMIN_RANGE[1])
            Parameter input file equivalent are `CHK_RMAX`, `CHK_QMIN_RANGE[1]`
        qmin_coord : str, {"HELIO", "BARY"}, default "HELIO"
            coordinate frame to use for checking the minimum periapsis distance
            Parameter input file equivalent is `QMIN_COORD`
        mtiny : float, optional
            The minimum mass of fully interacting bodies. Bodies below this mass interact with the larger bodies,
            but not each other (SyMBA only). Only mtiny or gmtiny is accepted, not both.
            Parameter input file equivalent is None
        gmtiny : float, optional
            The minimum G*mass of fully interacting bodies. Bodies below this mass interact with the larger bodies,
            but not each other (SyMBA only). Only mtiny or gmtiny is accepted, not both.
            Parameter input file equivalent is `GMTINY`
        close_encounter_check : bool, default True
            Check for close encounters between bodies. If set to True, then the radii of massive bodies must be included
            in initial conditions.
            Parameter input file equivalent is `CHK_CLOSE`
        encounter_save : {"NONE","TRAJECTORY","CLOSEST", "BOTH"}, default "NONE"
            Indicate if and how encounter data should be saved. If set to "TRAJECTORY", the position and velocity vectors
            of all bodies undergoing close encounters are saved at each intermediate step to the encounter files.
            If set to "CLOSEST", the position  and velocities at the point of closest approach between pairs of bodies are
            computed and stored to the encounter files. If set to "BOTH", then this stores the values that would be computed
            in "TRAJECTORY" and "CLOSEST". If set to "NONE" no trajectory information is saved.
            WARNING - Enabling this feature could lead to very large files.
        general_relativity : bool, default True
            Include the post-Newtonian correction in acceleration calculations.
            Parameter input file equivalent is "GR"
        collision_model : {"MERGE","BOUNCE","FRAGGLE"}, default "MERGE"
            This is used to set the collision/fragmentation model.
            This argument only applies to Swiftest-SyMBA simulations and will be ignored otherwise.
            Parameter input file equivalent is "COLLISION_MODEL"
        minimum_fragment_gmass : float, optional
            If fragmentation is turned on, this sets the mimimum G*mass of a collisional fragment that can be generated if a
            fragmentation model is enabled. Ignored otherwise.
            Only set one of minimum_fragment_gmass or minimum_fragment_mass
            Parameter input file equivalent is None
        minimum_fragment_mass : float, optional
            If fragmentation is turned on, this sets the mimimum mass of a collisional fragment that can be generated. if a
            fragmentation model is enabled. Ignored otherwise
            Only set one of minimum_fragment_gmass or minimum_fragment_mass
            Parameter input file equivalent is `MIN_GMFRAG`
        nfrag_reduction : float, optional
            If fragmentation is turned on, this is a reduction factor used to limit the number of fragments generated in a collision.
            For instance, if the SFD of the collision would generated 300 fragments above the `minimum_fragment_mass`, then a value
            of `nfrag_reduction = 30.0` would reduce it to 10.
            Currently only used by the Fraggle collision model.
        rotation : bool, default True
            If set to True, this turns on rotation tracking and radius, rotation vector, and moments of inertia values
            must be included in the initial conditions.
            This argument only applies to Swiftest-SyMBA simulations. It will be ignored otherwise.
            Parameter input file equivalent is `ROTATION`
        compute_conservation_values : bool, default True for SyMBA, otherwise False
            Turns on the computation of energy, angular momentum, and mass conservation and reports the values
            every output step of a running simulation.
            Parameter input file equivalent is `ENERGY`
        extra_force : bool, default False
            Turns on user-defined force function.
            Parameter input file equivalent is `EXTRA_FORCE`
        big_discard : bool, default False
            Includes big bodies when performing a discard (Swifter only)
            Parameter input file equivalent is `BIG_DISCARD`
        rhill_present : bool, default False
            Include the Hill's radius with the input files .
            Parameter input file equivalent is `RHILL_PRESENT`
        restart : bool, default False
            If true, will restart an old run. The file given by `output_file_name` must exist before the run can
            execute. If false, will start a new run. If the file given by `output_file_name` exists, it will be replaced
            when the run is executed.
            Parameter input file equivalent is `OUT_STAT`
        interaction_loops : {"TRIANGULAR","FLAT"}, default "TRIANGULAR"
            *Swiftest Experimental feature*
            Specifies which algorithm to use for the computation of body-body gravitational forces.

            * "TRIANGULAR" - Upper-triangular double-loops.
            * "FLAT" - Body-body interation pairs are flattened into a 1-D array.

            Parameter input file equivalent is `INTERACTION_LOOPS`
        encounter_check_loops : {"TRIANGULAR","SORTSWEEP"}, default "TRIANGULAR"
            *Swiftest Experimental feature*
            Specifies which algorithm to use for checking whether bodies are in a close encounter state or not.

            * "TRIANGULAR" - Upper-triangular double-loops.
            * "SORTSWEEP" - A Sort-Sweep algorithm is used to reduce the population of potential close encounter bodies.
              This algorithm is still in development, and does not necessarily speed up the encounter checking.
              Use with caution.

            Parameter input file equivalent is `ENCOUNTER_CHECK`
        dask : bool, default False
            Use Dask to lazily load data (useful for very large datasets)
        coarray : bool, default False
            If true, will employ Coarrays on test particle structures to run in single program/multiple data parallel mode.
            In order to use this capability, Swiftest must be compiled for Coarray support. Only certain integrators can use
            Coarrays. RMVS, WHM, Helio are all compatible, but SyMBA is not, due to the way tp-pl close encounters are handeled.
        verbose : bool, default False
            If set to True, then more information is printed by Simulation methods as they are executed. Setting to
            False suppresses most messages other than errors and some warnings.
        seed : ArrayLike of int, optional
            Seed for the random number generator. If not set, then a random seed is generated.
            Parameter input file equivalent is `SEED`

        Returns
        -------
        param : dict
            A dictionary of all Simulation parameters that changed
        """
        simdir = kwargs.pop("simdir", None)
        if simdir is not None:
            self.simdir = simdir
        param_file = kwargs.pop("param_file", None)
        if param_file is not None:
            self.param_file = param_file

        if "verbose" in kwargs:
            self.verbose = kwargs.pop("verbose")

        valid_args = list(self._param_to_argument.values()) + [
            "codename",
            "integrator",
            "coarray",
            "verbose",
            "dask",
            "MU",
            "DU",
            "TU",
            "MU_name",
            "DU_name",
            "TU_name",
            "tstep_out",
            "nstep_out",
            "mtiny",
            "minimum_fragment_mass",
            "ephemeris_date",
        ]

        unrecognized = [k for k, v in kwargs.items() if k not in valid_args]
        if len(unrecognized) > 0:
            for k in unrecognized:
                warnings.warn(f'Unrecognized argument "{k}"', stacklevel=2)

        # Setters returning parameter dictionary values
        param_dict = {}
        param_dict.update(self.set_unit_system(**kwargs))
        param_dict.update(self.set_integrator(**kwargs))
        param_dict.update(self.set_simulation_time(**kwargs))
        param_dict.update(self.set_init_cond_files(**kwargs))
        param_dict.update(self.set_output_files(**kwargs))
        param_dict.update(self.set_distance_range(**kwargs))
        param_dict.update(self.set_feature(**kwargs))

        # Non-returning setters
        if "ephemeris_date" in kwargs:
            self.ephemeris_date = kwargs["ephemeris_date"]

        return param_dict

    def get_parameter(self, **kwargs: Any) -> dict[str, Any]:
        """
        Setter for all possible parameters. Calls each of the specialized setters using keyword arguments.

        Parameters
        ----------
        **kwargs : Any
            Any of the arguments defined in Simulation. If none provided, it returns all arguments.

        Returns
        -------
        param : dict
            A dictionary of all Simulation parameters requested
        """
        # Getters returning parameter dictionary values
        param_dict = {}
        param_dict.update(self.get_unit_system(**kwargs))
        param_dict.update(self.get_integrator(**kwargs))
        param_dict.update(self.get_simulation_time(**kwargs))
        param_dict.update(self.get_init_cond_files(**kwargs))
        param_dict.update(self.get_output_files(**kwargs))
        param_dict.update(self.get_distance_range(**kwargs))
        param_dict.update(self.get_feature(**kwargs))

        self.get_ephemeris_date(**kwargs)

        return param_dict

    def set_integrator(
        self,
        codename: None | Literal["Swiftest", "Swifter", "Swift"] = "Swiftest",
        integrator: Literal["symba", "rmvs", "whm", "helio"] | None = None,
        mtiny: float | None = None,
        gmtiny: float | None = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """
        Sets the integrator to be used when running a simulation. Returns a dictionary of the parameters that were set.

        Parameters
        ----------
        codename : {"swiftest", "swifter", "swift"}, optional
            The name of the code to use. Case-insensitive valid options are swiftest, swifter, and swift. Currently only swiftest is
            is supported for excuting runs with the run() method.
        integrator : {"symba","rmvs","whm","helio"}, optional
            Name of the n-body integrator that will be used when executing a run.
        mtiny : float, optional
            The minimum mass of fully interacting bodies. Bodies below this mass interact with the larger bodies,
            but not each other (SyMBA only). Only mtiny or gmtiny is accepted, not both.
        gmtiny : float, optional
            The minimum G*mass of fully interacting bodies. Bodies below this mass interact with the larger bodies,
            but not each other (SyMBA only). Only mtiny or gmtiny is accepted, not both.
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        integrator_dict: dict
            A dictionary containing the subset of the parameter dictonary that was updated by this setter

        """
        verbose = kwargs.pop("verbose", self.verbose)
        update_list = []

        if codename is not None:
            self.codename = codename
            self.param["! VERSION"] = f"{self.codename} input file"
            update_list.append("codename")

        if self.codename == "Swifter":
            J2 = self.param.pop("J2", 0.0)
            J4 = self.param.pop("J4", 0.0)
            self.param = io.swiftest2swifter_param(self.param, J2, J4)

        if integrator is not None:
            self.integrator = integrator

        if mtiny is not None or gmtiny is not None:
            if verbose:
                if self.integrator != "symba":
                    warnings.warn("mtiny and gmtiny are only used by SyMBA.", stacklevel=2)
                if mtiny is not None and gmtiny is not None:
                    warnings.warn("Only set mtiny or gmtiny, not both. Using gmtiny", stacklevel=2)

            if gmtiny is not None:
                self.param["GMTINY"] = gmtiny
                update_list.append("gmtiny")
            elif mtiny is not None:
                self.param["GMTINY"] = self.GU * mtiny
                update_list.append("gmtiny")

        integrator_dict = self.get_integrator(update_list, verbose=verbose)

        return integrator_dict

    def get_integrator(self, arg_list: str | list[str] | None = None, **kwargs: Any) -> dict[str, Any]:
        """
        Returns a subset of the parameter dictionary containing the current values of the distance range parameters.

        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list: str | List[str], optional
            A single string or list of strings containing the names of the features to extract. Default is all of:
            ["integrator"]
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        integrator_dict : dict
            The subset of the dictionary containing the code name if codename is selected
        """
        verbose = kwargs.pop("verbose", self.verbose)

        valid_arg = ["gmtiny"]
        valid_var = self._create_valid_var(valid_arg)

        valid_instance_vars = {
            "codename": self.codename,
            "integrator": self.integrator,
            "param_file": str(self.param_file),
        }

        try:
            self.integrator
        except:
            if verbose:
                warnings.warn("integrator is not set", stacklevel=2)
            return {}

        try:
            self.codename
        except:
            if verbose:
                warnings.warn("codename is not set", stacklevel=2)
            return {}

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

    def set_feature(
        self,
        close_encounter_check: bool | None = None,
        general_relativity: bool | None = None,
        collision_model: Literal["MERGE", "BOUNCE", "FRAGGLE"] | None = None,
        minimum_fragment_gmass: float | None = None,
        minimum_fragment_mass: float | None = None,
        nfrag_reduction: float | None = None,
        rotation: bool | None = None,
        compute_conservation_values: bool | None = None,
        extra_force: bool | None = None,
        big_discard: bool | None = None,
        rhill_present: bool | None = None,
        tides: bool | None = None,
        interaction_loops: Literal["TRIANGULAR", "FLAT"] | None = None,
        encounter_check_loops: Literal["TRIANGULAR", "SORTSWEEP"] | None = None,
        encounter_save: Literal["NONE", "TRAJECTORY", "CLOSEST", "BOTH"] | None = None,
        coarray: bool | None = None,
        simdir: str | os.PathLike = None,
        seed: int | list[int] | npt.NDArray[np.int_] | None = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
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
            WARNING - Enabling this feature could lead to very large files.
        general_relativity : bool, optional
            Include the post-Newtonian correction in acceleration calculations.
        collision_model : {"MERGE","BOUNCE","FRAGGLE"}, default "MERGE"
            This is used to set the collision/fragmentation model.  This argument only applies to Swiftest-SyMBA simulations. It
            will be ignored otherwise. Parameter input file equivalent is `COLLISION_MODEL`
        minimum_fragment_gmass : float, optional
            If fragmentation is turned on, this sets the mimimum G*mass of a collisional fragment that can be generated if a
            fragmentation model is enabled. Ignored otherwise.
            Only set one of minimum_fragment_gmass or minimum_fragment_mass
            Parameter input file equivalent is None
        minimum_fragment_mass : float, optional
            If fragmentation is turned on, this sets the mimimum mass of a collisional fragment that can be generated. if a
            fragmentation model is enabled. Ignored otherwise
            Only set one of minimum_fragment_gmass or minimum_fragment_mass
            Parameter input file equivalent is `MIN_GMFRAG`
        nfrag_reduction : float, optional
            If fragmentation is turned on, this is a reduction factor used to limit the number of fragments generated in a collision.
            For instance, if the SFD of the collision would generated 300 fragments above the `minimum_fragment_mass`, then a value
            of `nfrag_reduction = 30.0` would reduce it to 10.
            Currently only used by the Fraggle collision model.
        rotation : bool, optional
            If set to True, this turns on rotation tracking and radius, rotation vector, and moments of inertia values
            must be included in the initial conditions.
            This argument only applies to Swiftest-SyMBA simulations. It will be ignored otherwise.
        compute_conservation_values : bool, optional
            Turns on the computation of energy, angular momentum, and mass conservation and reports the values
            every output step of a running simulation.
        extra_force : bool, optional
            Turns on user-defined force function.
        big_discard : bool, optional
            Includes big bodies when performing a discard (Swifter only)
        rhill_present : bool, optional
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
            In order to use this capability, Swiftest must be compiled for Coarray support. Only certain integrators
            can use Coarrays: RMVS, WHM, Helio are all compatible, but SyMBA is not, due to the way tp-pl close encounters
            are handeled.
        tides : bool, optional
            Turns on tidal model (IN DEVELOPMENT - IGNORED)
        Yarkovsky : bool, optional
            Turns on Yarkovsky model (IN DEVELOPMENT - IGNORED)
        YORP : bool, optional
            Turns on YORP model (IN DEVELOPMENT - IGNORED)
        simdir : PathLike, optional
            Directory where simulation data will be stored, including the parameter file, initial conditions file, output file,
            dump files, and log files.
        seed : ArrayLike of int, optional
            Seed for the random number generator. If not set, then a random seed is generated.
            Parameter input file equivalent is `SEED`
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        feature_dict : dict
            A dictionary containing the requested features.
        """
        verbose = kwargs.pop("verbose", self.verbose)

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
            self.simdir = simdir

        if self.codename == "Swiftest":
            if general_relativity is not None:
                self.param["GR"] = general_relativity
                update_list.append("general_relativity")

            fragmentation_models = ["FRAGGLE"]
            if collision_model is not None:
                collision_model = collision_model.upper()
                fragmentation = collision_model in fragmentation_models
                if self.codename != "Swiftest" and self.integrator != "symba" and fragmentation:
                    if verbose:
                        warnings.warn(
                            "Fragmentation is only available on Swiftest SyMBA.",
                            stacklevel=2,
                        )
                    self.param["COLLISION_MODEL"] = "MERGE"
                else:
                    self.param["COLLISION_MODEL"] = collision_model
                    update_list.append("collision_model")
                    if fragmentation:
                        if "MIN_GMFRAG" not in self.param and minimum_fragment_mass is None and minimum_fragment_gmass is None:
                            if verbose:
                                warnings.warn(
                                    "Minimum fragment mass is not set. Set it using minimum_fragment_gmass or minimum_fragment_mass",
                                    stacklevel=2,
                                )
                        else:
                            update_list.append("minimum_fragment_gmass")

            if verbose:
                if minimum_fragment_gmass is not None and minimum_fragment_mass is not None:
                    warnings.warn(
                        "Only set either minimum_fragment_mass or minimum_fragment_gmass, but not both!",
                        stacklevel=2,
                    )

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
                if self.integrator == "symba":
                    if not rotation:
                        if verbose:
                            warnings.warn(
                                "Rotation is on by default for SyMBA. This option is ignored",
                                stacklevel=2,
                            )
                    self.param["ROTATION"] = True
                else:
                    self.param["ROTATION"] = rotation
                update_list.append("rotation")

            if self.param["COLLISION_MODEL"] == "FRAGGLE" and not self.param["ROTATION"]:
                self.param["ROTATION"] = True
                update_list.append("rotation")

            if self.integrator == "symba":
                self.param["ENERGY"] = True
            else:
                self.param["ENERGY"] = False

            if compute_conservation_values is not None:
                if self.integrator == "symba":
                    if not compute_conservation_values:
                        if verbose:
                            warnings.warn(
                                "Energy, angular momentum, and mass conservation values are computed by default for SyMBA. This option is ignored",
                                stacklevel=2,
                            )
                    self.param["ENERGY"] = True
                else:
                    self.param["ENERGY"] = compute_conservation_values
                update_list.append("compute_conservation_values")

            if interaction_loops is not None:
                valid_vals = ["TRIANGULAR", "FLAT"]
                interaction_loops = interaction_loops.upper()
                if interaction_loops not in valid_vals:
                    if verbose:
                        msg = f"{interaction_loops} is not a valid option for interaction loops."
                        msg += f"\nMust be one of {valid_vals}"
                        warnings.warn(msg, stacklevel=2)
                    if "INTERACTION_LOOPS" not in self.param:
                        self.param["INTERACTION_LOOPS"] = valid_vals[0]
                else:
                    self.param["INTERACTION_LOOPS"] = interaction_loops
                    update_list.append("interaction_loops")

            if encounter_check_loops is not None:
                valid_vals = ["TRIANGULAR", "SORTSWEEP"]
                encounter_check_loops = encounter_check_loops.upper()
                if encounter_check_loops not in valid_vals:
                    if verbose:
                        msg = f"{encounter_check_loops} is not a valid option for interaction loops."
                        msg += f"\nMust be one of {valid_vals}"
                        warnings.warn(msg, stacklevel=2)
                    if "ENCOUNTER_CHECK" not in self.param:
                        self.param["ENCOUNTER_CHECK"] = valid_vals[0]
                else:
                    self.param["ENCOUNTER_CHECK"] = encounter_check_loops
                    update_list.append("encounter_check_loops")

            if encounter_save is not None:
                valid_vals = ["NONE", "TRAJECTORY", "CLOSEST", "BOTH"]
                encounter_save = encounter_save.upper()
                if encounter_save not in valid_vals:
                    if verbose:
                        msg = f"{encounter_save} is not a valid option for encounter_save."
                        msg += f"\nMust be one of {valid_vals}"
                        warnings.warn(msg, stacklevel=2)
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

            if seed is not None:
                if self.codename == "Swiftest":
                    if isinstance(seed, (int, np.integer)):
                        seed = [seed]
                    seed = np.array(seed, dtype=np.int64)
                    self.param["SEED"] = seed
                    update_list.append("seed")

        feature_dict = self.get_feature(update_list, verbose=verbose)
        return feature_dict

    def get_feature(self, arg_list: str | list[str] | None = None, **kwargs: Any) -> dict[str, Any]:
        """
        Returns a subset of the parameter dictionary containing the current value of the feature boolean values.

        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list : str | List[str], optional
            A single string or list of strings containing the names of the features to extract. Default is all of:
            ["close_encounter_check", "general_relativity", "collision_model", "rotation", "compute_conservation_values"]
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        feature_dict : dict
           A dictionary containing the requested features.
        """
        verbose = kwargs.pop("verbose", self.verbose)

        valid_arg = [
            "close_encounter_check",
            "collision_model",
            "encounter_save",
            "minimum_fragment_gmass",
            "nfrag_reduction",
            "rotation",
            "general_relativity",
            "compute_conservation_values",
            "rhill_present",
            "extra_force",
            "big_discard",
            "interaction_loops",
            "encounter_check_loops",
            "coarray",
            "seed",
        ]

        valid_var = self._create_valid_var(valid_arg)

        valid_arg, feature_dict = self._get_valid_arg_list(arg_list, valid_var)

        if verbose:
            for arg in valid_arg:
                key = valid_var[arg]
                if key in feature_dict:
                    if arg == "minimum_fragment_gmass":
                        print(f"{arg:<{self._getter_column_width}} {feature_dict[key]} {self.DU_name}^3 / {self.TU_name}^2 ")
                        print(
                            f"{'minimum_fragment_mass':<{self._getter_column_width}} {feature_dict[key] / self.GU} {self.MU_name}"
                        )
                    else:
                        print(f"{arg:<{self._getter_column_width}} {feature_dict[key]}")
                else:
                    print(f"{arg:<{self._getter_column_width}} NOT SET")

        return feature_dict

    def set_init_cond_files(
        self,
        init_cond_file_type: Literal["NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"] | None = None,
        init_cond_file_name: str | os.PathLike | dict[str, str] | dict[str, os.PathLike] | None = None,
        init_cond_format: Literal["EL", "XV"] | None = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """
        Sets the initial condition file parameters in the parameters dictionary.

        Parameters
        ----------
        init_cond_file_type : {"NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"}, optional
            The file type containing initial conditions for the simulation

            * NETCDF_DOUBLE: A single initial conditions input file in NetCDF file format of type NETCDF_DOUBLE
            * NETCDF_FLOAT: A single initial conditions input file in NetCDF file format of type NETCDF_FLOAT
            * ASCII : Three initial conditions files in ASCII format. The individual files define the central body, massive body, and test particle initial conditions.

        init_cond_file_name : str, path-like, or dict, optional
            Name of the input initial condition file or files. Whether to pass a single file name or a dictionary
            depends on the argument passed to `init_cond_file_type`: If `init_cond_file_type={"NETCDF_DOUBLE","NETCDF_FLOAT"}`,
            then this will be a single file name. If `init_cond_file_type="ASCII"` then this must be a dictionary where::

                init_cond_file_name = {
                    "CB" : *path to central body initial conditions file* (Swiftest only),
                    "PL" : *path to massive body initial conditions file*,
                    "TP" : *path to test particle initial conditions file*
                    }

            If no file name is provided then the following default file names will be used.

            * NETCDF_DOUBLE, NETCDF_FLOAT: `init_cond_file_name = "init_cond.nc"`
            * ASCII: `init_cond_file_name = {"CB" : "cb.in", "PL" : "pl.in", "TP" : "tp.in"}`

        init_cond_format : {"EL", "XV"}
            Indicates whether the input initial conditions are given as orbital elements or cartesian position and
            velocity vectors.
            If `codename` is "Swift" or "Swifter", EL initial conditions are converted to XV.
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        init_cond_file_dict : dict
           A dictionary containing the requested parameters
        """
        verbose = kwargs.pop("verbose", self.verbose)

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
            msg = "in set_init_cond_files: init_cond_file_name must be a dictionary of the form: "
            msg += "\n {"
            if codename == "Swiftest":
                msg += '\n"CB" : *path to central body initial conditions file*,'
            msg += '\n"PL" : *path to massive body initial conditions file*,'
            msg += '\n"TP" : *path to test particle initial conditions file*'
            msg += "\n}"
            warnings.warn(msg, stacklevel=2)
            return {}

        if init_cond_format is None:
            if "IN_FORM" in self.param:
                init_cond_format = self.param["IN_FORM"]
            else:
                init_cond_format = "XV"

        if init_cond_file_type is None:
            if "IN_TYPE" in self.param:
                init_cond_file_type = self.param["IN_TYPE"]
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
                if verbose:
                    warnings.warn(
                        f"{init_cond_file_type} is not supported by {self.codename}. Using ASCII instead",
                        stacklevel=2,
                    )
                init_cond_file_type = "ASCII"
            if init_cond_format != "XV":
                if verbose:
                    warnings.warn(
                        f"{init_cond_format} is not supported by {self.codename}. Using XV instead",
                        stacklevel=2,
                    )
                init_cond_format = "XV"

        valid_formats = {"EL", "XV"}
        if init_cond_format not in valid_formats:
            if verbose:
                warnings.warn(f"{init_cond_format} is not a valid input format", stacklevel=2)
        else:
            self.param["IN_FORM"] = init_cond_format

        valid_types = {"NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"}
        if init_cond_file_type not in valid_types:
            if verbose:
                warnings.warn(f"{init_cond_file_type} is not a valid input type", stackevel=2)
        else:
            self.param["IN_TYPE"] = init_cond_file_type

        if init_cond_file_type == "ASCII":
            if init_cond_file_name is None:
                # No file names passed, so we will just use the defaults.
                for key in init_cond_keys:
                    self.param[f"{key}_IN"] = f"{key.lower()}.in"
            elif type(init_cond_file_name) is not dict:
                # Oops, accidentally passed a single string or path-like instead of the expected dictionary for ASCII
                # input type.
                if verbose:
                    ascii_file_input_error_msg(self.codename)
            elif not all(key in init_cond_file_name for key in init_cond_keys):
                # This is the case where the dictionary doesn't have all the keys we expect. Print an error message.
                if verbose:
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
                if verbose:
                    warnings.warn(
                        "Only a single input file is used for NetCDF files",
                        stacklevel=2,
                    )
            else:
                self.param["NC_IN"] = init_cond_file_name

        init_cond_file_dict = self.get_init_cond_files(update_list)

        return init_cond_file_dict

    def get_init_cond_files(self, arg_list: str | list[str] | None = None, **kwargs: Any) -> dict[str, Any]:
        """
        Returns a subset of the parameter dictionary containing the current initial condition file parameters.

        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list : str | List[str], optional
            A single string or list of strings containing the names of the simulation time parameters to extract.
            Default is all of:
            ["init_cond_file_type", "init_cond_file_name", "init_cond_format"]
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        init_cond_file_dict : dict
           A dictionary containing the requested parameters
        """
        verbose = kwargs.pop("verbose", self.verbose)

        valid_arg = [
            "init_cond_file_type",
            "init_cond_format",
            "init_cond_file_name",
            "init_cond_file_name['CB']",
            "init_cond_file_name['PL']",
            "init_cond_file_name['TP']",
        ]
        valid_var = self._create_valid_var(valid_arg)

        three_file_args = [
            "init_cond_file_name['CB']",
            "init_cond_file_name['PL']",
            "init_cond_file_name['TP']",
        ]

        if self.codename == "Swifter":
            three_file_args.remove("init_cond_file_name['CB']")
            valid_var.pop("init_cond_file_name['CB']", None)

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

        if verbose:
            for arg in valid_arg:
                key = valid_var[arg]
                print(f"{arg:<{self._getter_column_width}} {init_cond_file_dict[key]}")

        return init_cond_file_dict

    def set_output_files(
        self,
        output_file_type: Literal["NETCDF_DOUBLE", "NETCDF_FLOAT", "REAL4", "REAL8", "XDR4", "XDR8"] | None = None,
        output_file_name: os.PathLike | str | None = None,
        output_format: Literal["XV", "XVEL"] | None = None,
        restart: bool | None = None,
        simdir: os.PathLike | str | None = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
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
        restart : bool, optional
            Indicates whether this is a restart of an old run or a new run.
        simdir : PathLike, optional
            Directory where simulation data will be stored, including the parameter file, initial conditions file, output file
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        output_file_dict : dict
           A dictionary containing the requested parameters
        """
        verbose = kwargs.pop("verbose", self.verbose)

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
        if simdir is not None:
            self.simdir = simdir
            update_list.append("simdir")
        if len(update_list) == 0:
            return {}

        if self.codename == "Swiftest":
            if output_file_type is None:
                output_file_type = self.param.pop("OUT_TYPE", None)
                if output_file_type is None:
                    output_file_type = "NETCDF_DOUBLE"
            elif output_file_type not in ["NETCDF_DOUBLE", "NETCDF_FLOAT"]:
                if verbose:
                    warnings.warn(
                        f"{output_file_type} is not compatible with Swiftest. Setting to NETCDF_DOUBLE",
                        stacklevel=2,
                    )
                output_file_type = "NETCDF_DOUBLE"
        elif self.codename == "Swifter":
            if output_file_type is None:
                output_file_type = self.param.pop("OUT_TYPE", None)
                if output_file_type is None:
                    output_file_type = "REAL8"
            elif output_file_type not in ["REAL4", "REAL8", "XDR4", "XDR8"]:
                if verbose:
                    warnings.warn(
                        f"{output_file_type} is not compatible with Swifter. Setting to REAL8",
                        stacklevel=2,
                    )
                output_file_type = "REAL8"
        elif self.codename == "Swift":
            if output_file_type is None:
                output_file_type = self.param.pop("OUT_TYPE", None)
                if output_file_type is None:
                    output_file_type = "REAL4"
            if output_file_type not in ["REAL4"]:
                if verbose:
                    warnings.warn(
                        f"{output_file_type} is not compatible with Swift. Setting to REAL4",
                        stacklevel=2,
                    )
                output_file_type = "REAL4"

        self.param["OUT_TYPE"] = output_file_type
        if output_file_name is None:
            if output_file_type in ["NETCDF_DOUBLE", "NETCDF_FLOAT"]:
                self.param["BIN_OUT"] = "data.nc"
            else:
                self.param["BIN_OUT"] = "bin.dat"
        else:
            self.param["BIN_OUT"] = output_file_name

        if output_format is not None:
            if output_format != "XV" and self.codename != "Swiftest":
                if verbose:
                    warnings.warn(
                        f"{output_format} is not compatible with {self.codename}. Setting to XV",
                        stacklevel=2,
                    )
                output_format = "XV"
            self.param["OUT_FORM"] = output_format

        if self.restart:
            self.param["OUT_STAT"] = "APPEND"
        else:
            self.param["OUT_STAT"] = "REPLACE"

        output_file_dict = self.get_output_files(update_list)

        return output_file_dict

    def get_output_files(self, arg_list: str | list[str] | None = None, **kwargs: Any) -> dict[str, Any]:
        """
        Returns a subset of the parameter dictionary containing the current output file parameters.

        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list : str | List[str], optional
            A single string or list of strings containing the names of the simulation time parameters to extract.
            Default is all of:
            ["output_file_type", "output_file_name", "output_format"]
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        output_file_dict : dict
           A dictionary containing the requested parameters

        """
        verbose = kwargs.pop("verbose", self.verbose)

        valid_arg = [
            "output_file_type",
            "output_file_name",
            "output_format",
            "restart",
        ]

        valid_var = self._create_valid_var(valid_arg)

        valid_instance_vars = {"simdir": self.simdir}

        if not bool(kwargs) and arg_list is None:
            arg_list = list(valid_instance_vars.keys())
            arg_list += [a for a in valid_var.keys() if a not in valid_instance_vars]

        valid_arg, output_file_dict = self._get_valid_arg_list(arg_list, valid_var)

        if verbose:
            for arg in valid_arg:
                key = valid_var[arg]
                print(f"{arg:<{self._getter_column_width}} {output_file_dict[key]}")

        return output_file_dict

    def set_unit_system(
        self,
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
        **kwargs: Any,
    ) -> dict[str, Any]:
        """
        Setter for setting the unit conversion between one of the standard sets.

        The units can be set one of two ways

        #. The user can supply string values to the arguments MU, DU, and TU to select between common systems
        #. The user can supply float values to the arguments MU2KG, DU2M, and TU2S to manually set the conversion
           factor between the desired unit and the SI unit (kg-m-s).

        The two sets of arguments are mutually exclusive. Any values passed to MU2KG, DU2M, or TU2S will override any
        specified in MU, DU, or TU, respectively. The default system is Msun-AU-YR. MU, DU, and TU are case-insenstive

        Parameters
        ----------
        MU : str, optional
           The mass unit system to use. Case-insensitive valid options are

           * "Msun"   - Solar mass
           * "Mearth" - Earth mass
           * "kg"     - kilograms
           * "g","gm" - grams

        DU : str, optional
            The distance unit system to use. Case-insensitive valid options are

            * "AU"     - Astronomical Unit
            * "Rearth" - Earth radius
            * "km"     - kilometer
            * "m"      - meter
            * "cm"     - centimeter

        TU : str, optional
            The time unit system to use. Case-insensitive valid options are

            * "y","YR","year","years"      - Year
            * "d","day","days"             - Julian day
            * "s","sec","seconds","second" - second

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
            This is a destructive operation, however if not executed then the values contained in the parameter
            file and input/output data files computed previously may not be consistent with the new unit conversion
            factors.
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        unit_dict : dict
           A dictionary containing the requested unit conversion parameters
        """
        verbose = kwargs.pop("verbose", self.verbose)

        MU2KG_old = None
        DU2M_old = None
        TU2S_old = None

        update_list = []
        if MU is not None or MU2KG is not None:
            update_list.append("MU")
        if DU is not None or DU2M is not None:
            update_list.append("DU")
        if TU is not None or TU2S is not None:
            update_list.append("TU")

        if MU2KG is not None or MU is not None:
            MU2KG_old = self.param.pop("MU2KG", None)
            if MU2KG is not None:
                self.param["MU2KG"] = MU2KG
                self.MU_name = None
            else:
                if MU.upper() == "MSUN":
                    self.param["MU2KG"] = constants.MSun
                    self.MU_name = "MSun"
                elif MU.upper() == "MEARTH":
                    self.param["MU2KG"] = constants.MEarth
                    self.MU_name = "MEarth"
                elif MU.upper() == "KG":
                    self.param["MU2KG"] = 1.0
                    self.MU_name = "kg"
                elif MU.upper() == "G" or MU.upper() == "GM":
                    self.param["MU2KG"] = 1e-3
                    self.MU_name = "g"
                else:
                    if verbose:
                        warnings.warn(
                            f"{MU} not a recognized unit system. Using MSun as a default.",
                            stacklevel=2,
                        )
                    self.param["MU2KG"] = constants.MSun
                    self.MU_name = "MSun"
            self.MU2KG = self.param["MU2KG"]
            self.KG2MU = 1.0 / self.MU2KG

        if DU2M is not None or DU is not None:
            DU2M_old = self.param.pop("DU2M", None)
            if DU2M is not None:
                self.param["DU2M"] = DU2M
                self.DU_name = None
            else:
                if DU.upper() == "AU":
                    self.param["DU2M"] = constants.AU2M
                    self.DU_name = "AU"
                elif DU.upper() == "REARTH":
                    self.param["DU2M"] = constants.REarth
                    self.DU_name = "REarth"
                elif DU.upper() == "KM":
                    self.param["DU2M"] = 1000.0
                    self.DU_name = "km"
                elif DU.upper() == "M":
                    self.param["DU2M"] = 1.0
                    self.DU_name = "m"
                elif DU.upper() == "CM":
                    self.param["DU2M"] = 1e-2
                    self.DU_name = "cm"
                else:
                    if verbose:
                        warnings.warn(
                            f"{DU} not a recognized unit system. Using AU as a default.",
                            stacklevel=2,
                        )
                    self.param["DU2M"] = constants.AU2M
                    self.DU_name = "AU"
            self.DU2M = self.param["DU2M"]
            self.M2DU = 1.0 / self.DU2M

        if TU2S is not None or TU is not None:
            TU2S_old = self.param.pop("TU2S", None)
            if TU2S is not None:
                self.param["TU2S"] = TU2S
                self.TU_name = None
            else:
                if TU.upper() == "YR" or TU.upper() == "Y" or TU.upper() == "YEAR" or TU.upper() == "YEARS":
                    self.param["TU2S"] = constants.YR2S
                    self.TU_name = "y"
                elif TU.upper() == "DAY" or TU.upper() == "D" or TU.upper() == "JD" or TU.upper() == "DAYS":
                    self.param["TU2S"] = constants.JD2S
                    self.TU_name = "d"
                elif TU.upper() == "S" or TU.upper() == "SECONDS" or TU.upper() == "SEC" or TU.upper() == "SECOND":
                    self.param["TU2S"] = 1.0
                    self.TU_name = "s"
                else:
                    warnings.warn(
                        f"{TU} not a recognized unit system. Using YR as a default.",
                        stacklevel=2,
                    )
                    self.param["TU2S"] = constants.YR2S
                    self.TU_name = "y"
            self.TU2S = self.param["TU2S"]
            self.S2TU = 1.0 / self.TU2S

        if MU_name is not None:
            self.MU_name = MU_name
        if DU_name is not None:
            self.DU_name = DU_name
        if TU_name is not None:
            self.TU_name = TU_name

        if "DU_name" in dir(self) and "TU_name" in dir(self):
            self.VU_name = f"{self.DU_name}/{self.TU_name}"
        if all(key in self.param for key in ["MU2KG", "DU2M", "TU2S"]):
            self.GU = constants.GC * self.param["TU2S"] ** 2 * self.param["MU2KG"] / self.param["DU2M"] ** 3

        if "MU2KG" in self.param and "DU2M" in self.param and "TU2S" in self.param:
            if (
                recompute_unit_values
                and MU2KG_old != self.param["MU2KG"]
                or DU2M_old != self.param["DU2M"]
                or TU2S_old != self.param["TU2S"]
            ):
                self._update_param_units(MU2KG_old, DU2M_old, TU2S_old)

        unit_dict = self.get_unit_system(update_list, verbose=verbose)

        return unit_dict

    def get_unit_system(self, arg_list: str | list[str] | None = None, **kwargs) -> dict[str, Any]:
        """
        Returns a subset of the parameter dictionary containing the current simulation unit system.

        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list : str | List[str], optional
            A single string or list of strings containing the names of the simulation unit system
            Default is all of ["MU", "DU", "TU"]
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        unit_dict : dict
           A dictionary containing the requested unit conversion parameters

        """
        verbose = kwargs.pop("verbose", self.verbose)

        valid_arg = ["MU", "DU", "TU"]
        valid_var = self._create_valid_var(valid_arg)

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

        units1 = {"MU": MU_name, "DU": DU_name, "TU": TU_name}
        units2 = {
            "MU": f"kg / {MU_name}",
            "DU": f"m / {DU_name}",
            "TU": f"s / {TU_name}",
        }

        valid_arg, unit_dict = self._get_valid_arg_list(arg_list, valid_var)

        if verbose:
            for arg in valid_arg:
                key = valid_var[arg]
                col_width = str(int(self._getter_column_width) - 4)
                print(f"{arg}: {units1[arg]:<{col_width}} {unit_dict[key]} {units2[arg]}")

        return unit_dict

    def _update_param_units(self, MU2KG_old, DU2M_old, TU2S_old):
        """
        Updates the values of parameters that have units when the units have changed.

        Parameters
        ----------
        MU2KG_old : float
            Old value of the mass unit conversion factor
        DU2M_old : float
            Old value of the distance unit conversion factor
        TU2S_old : float
            Old value of the time unit conversion factor

        Returns
        -------
        Updates the set of param dictionary values for the new unit system
        """
        mass_keys = ["GMTINY", "MIN_GMFRAG"]
        distance_keys = ["CHK_QMIN", "CHK_RMIN", "CHK_RMAX", "CHK_EJECT"]
        time_keys = ["T0", "TSTOP", "DT"]

        if MU2KG_old is not None:
            for k in mass_keys:
                if k in self.param:
                    self.param[k] *= MU2KG_old / self.param["MU2KG"]

        if DU2M_old is not None:
            for k in distance_keys:
                if k in self.param:
                    self.param[k] *= DU2M_old / self.param["DU2M"]

            CHK_QMIN_RANGE = self.param.pop("CHK_QMIN_RANGE", None)
            if CHK_QMIN_RANGE is not None:
                CHK_QMIN_RANGE = CHK_QMIN_RANGE.split(" ")
                for i, v in enumerate(CHK_QMIN_RANGE):
                    if float(v) > 0.0:
                        CHK_QMIN_RANGE[i] = float(v) * DU2M_old / self.param["DU2M"]
                self.param["CHK_QMIN_RANGE"] = f"{CHK_QMIN_RANGE[0]} {CHK_QMIN_RANGE[1]}"

        if TU2S_old is not None:
            for k in time_keys:
                if k in self.param:
                    self.param[k] *= TU2S_old / self.param["TU2S"]

        return

    def set_distance_range(
        self,
        rmin: float | None = None,
        rmax: float | None = None,
        qmin_coord: Literal["HELIO", "BARY"] | None = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
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
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        range_dict : dict
           A dictionary containing the requested parameters.
        """
        verbose = kwargs.pop("verbose", self.verbose)

        if rmax is None and rmin is None and qmin_coord is None:
            return {}

        update_list = []
        CHK_QMIN_RANGE = self.param.pop("CHK_QMIN_RANGE", None)
        if CHK_QMIN_RANGE is None:
            CHK_QMIN_RANGE = [-1, -1]
        else:
            CHK_QMIN_RANGE = CHK_QMIN_RANGE.split(" ")
        if rmin is not None:
            self.param["CHK_QMIN"] = rmin
            self.param["CHK_RMIN"] = rmin
            CHK_QMIN_RANGE[0] = rmin
            update_list.append("rmin")
        if rmax is not None:
            self.param["CHK_RMAX"] = rmax
            self.param["CHK_EJECT"] = rmax
            CHK_QMIN_RANGE[1] = rmax
            update_list.append("rmax")
        if qmin_coord is not None:
            valid_qmin_coord = ["HELIO", "BARY"]
            if qmin_coord.upper() not in valid_qmin_coord:
                if verbose:
                    warnings.warn(
                        f"qmin_coord = {qmin_coord} is not a valid option.  Must be one of",
                        ",".join(valid_qmin_coord),
                        stacklevel=2,
                    )
                self.param["CHK_QMIN_COORD"] = valid_qmin_coord[0]
            else:
                self.param["CHK_QMIN_COORD"] = qmin_coord.upper()
            update_list.append("qmin_coord")

        self.param["CHK_QMIN_RANGE"] = f"{CHK_QMIN_RANGE[0]} {CHK_QMIN_RANGE[1]}"

        range_dict = self.get_distance_range(update_list, verbose=verbose)

        return range_dict

    def get_distance_range(self, arg_list: str | list[str] | None = None, **kwargs: Any) -> dict[str, Any]:
        """
        Returns a subset of the parameter dictionary containing the current values of the distance range parameters.

        If the verbose option is set in the Simulation object, then it will also print the values.

        Parameters
        ----------
        arg_list : str | List[str], optional
            A single string or list of strings containing the names of the features to extract. Default is all of ["rmin", "rmax"]
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        range_dict : dict
           A dictionary containing the requested parameters.
        """
        verbose = kwargs.pop("verbose", self.verbose)

        valid_var = {
            "rmin": "CHK_RMIN",
            "rmax": "CHK_RMAX",
            "qmin_coord": "CHK_QMIN_COORD",
            "qmin": "CHK_QMIN",
            "qminR": "CHK_QMIN_RANGE",
        }

        units = {
            "rmin": self.DU_name,
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

    def add_solar_system_body(
        self,
        name: str | list[str] | None = None,
        ephemeris_id: int | list[int] | None = None,
        date: str | None = None,
        align_to_central_body_rotation: bool = False,
        **kwargs: Any,
    ) -> None:
        """
        Adds a solar system body to an existing simulation Dataset from the JPL Horizons ephemeris service.

        The JPL Horizons service will be searched for a body matching the string passed by `name`, or alternatively `ephemeris_id` if passed. Bodies will be named in the Swiftest initial conditions Dataset using `name`. Use `ephemeris_id` to have finer control over which body is searched in Horizons while using a custom name.

        If `name` is not passed, then the target name property is used as the name. You must pass either `name` and/or `ephemeris_id`.

        When passing `name` == "Earth" or `name` == "Pluto", it a body is generated that has initial conditions matching the system barycenter and mass equal to the sum of Earth+Moon or Pluto+Charon. To obtain initial conditions for either Earth or Pluto alone, pass `ephemeris_id` == "399" for Earth or `ephemeris_id` == "999" for Pluto.

        Parameters
        ----------
        name : str, optional | List[str]
            Add solar system body by name. This will be the name used in the Swiftest initial conditions Dataset unless not supplied
        ephemeris_id : int | List[int], optional but must be the same length as `name` if passed.
            In that case, the id is passed to the ephemeris service and the name is used.
        date : str, optional
            ISO-formatted date sto use when obtaining the ephemerides in the format YYYY-MM-DD. Defaults to value
            set by `ephemeris_date`.
        align_to_central_body_rotation : bool, default False
            If True, the cartesian coordinates will be aligned to the rotation pole of the central body. Otherwise, the This is only valid for when
            rotation is enabled.
        verbose : bool, default True
            If True, then warnings will be printed if the name is already in use in the Dataset.
        **kwargs : Any
            Additional keyword arguments to pass to the query method (i.e. astroquery.Horizons)

        Returns
        -------
        None
            initial conditions data stored as a SwiftestDataset in the init_cond instance variable
        """
        from .constants import CB_TYPE_NAME

        verbose = kwargs.pop("verbose", self.verbose)

        if name is None and ephemeris_id is None:
            if verbose:
                warnings.warn("Either `name` and/or `ephemeris_id` must be supplied to add_solar_system_body", stacklevel=2)
            return None
        if name is not None and (type(name) is str or type(name) is int):
            name = [name]

        if ephemeris_id is not None:
            if type(ephemeris_id) is int or type(ephemeris_id) is str:
                ephemeris_id = [ephemeris_id]
            if name is None:
                name = [None] * len(ephemeris_id)
            elif len(ephemeris_id) != len(name):
                raise ValueError(
                    f"The length of ephemeris_id ({len(ephemeris_id)}) does not match the length of name ({len(name)})"
                )
        else:
            ephemeris_id = [None] * len(name)

        if self.ephemeris_date is None:
            self.ephemeris_date = "MCBL"

        if date is None:
            date = self.ephemeris_date
        try:
            datetime.datetime.fromisoformat(date)
        except Exception:
            if verbose:
                warnings.warn(
                    f"{date} is not a valid date format. Must be 'YYYY-MM-DD'. Setting to {self.ephemeris_date}",
                    stacklevel=2,
                )
            date = self.ephemeris_date

        # Sun is the default central body
        if "particle_type" in self.data.variables and CB_TYPE_NAME in self.data.particle_type:
            cbname = self.data["name"].where(self.data.isel(time=0).particle_type == CB_TYPE_NAME, drop=True).values[0]
        else:
            cbname = "Sun"

        # Check to make sure we don't already have a body with the same name
        if "name" in self.data:
            bad_names = [n for n in name if n in self.data.name.values]
            if len(bad_names) > 0:
                name = [n for n in name if n not in bad_names]
                bad_names = ", ".join(bad_names)
                if verbose:
                    warnings.warn(
                        f"The following names are already in use and will not be added: {bad_names}",
                        stacklevel=2,
                    )

        if verbose:
            if None not in name:
                name_str = ", ".join(name)
                print(f"Adding solar system bodies by name: {name_str}")
            else:
                eph_str = ", ".join(ephemeris_id)
                print(f"Adding solar system bodies by id: {eph_str}")
        body_list = []
        for i, n in enumerate(name):
            body = init_cond.get_solar_system_body(
                name=n,
                ephemerides_start_date=date,
                ephemeris_id=ephemeris_id[i],
                central_body_name=cbname,
                verbose=verbose,
                **kwargs,
            )
            if body is None:
                if verbose:
                    warnings.warn(
                        f"Body with name {n} not found in the ephemerides. Skipping",
                        stacklevel=2,
                    )
                continue
            if "name" in self.data and body["name"] in self.data.name.values:
                if verbose:
                    warnings.warn(
                        f"Body with name {body['name']} already exists in the dataset. Skipping",
                        stacklevel=2,
                    )
                continue
            body_list.append(body)

        # Convert the list receieved from the get_solar_system_body output and turn it into arguments to vec2xr
        if len(body_list) == 0:
            if verbose:
                print("No valid bodies found")
            return
        else:
            vec2xr_kwargs = {}
            for d in body_list:
                if d is None:
                    continue
                for key, value in d.items():
                    if key not in vec2xr_kwargs:
                        vec2xr_kwargs[key] = [value]
                    else:
                        vec2xr_kwargs[key].append(value)

        scalar_floats = ["Gmass", "mass", "radius", "rhill", "j2rp2", "j4rp4"]
        vector_floats = ["rh", "vh", "Ip", "rot"]
        scalar_ints = ["id"]

        # Unit conversion factors
        for k, v in vec2xr_kwargs.items():
            if k in scalar_ints:
                v[v == None] = -1
                vec2xr_kwargs[k] = np.array(v, dtype=int)
            elif k in scalar_floats:
                vec2xr_kwargs[k] = np.array(v, dtype=np.float64)
                if all(np.isnan(vec2xr_kwargs[k])):
                    vec2xr_kwargs[k] = None
            elif k in vector_floats:
                vec2xr_kwargs[k] = np.vstack(v)
                vec2xr_kwargs[k] = vec2xr_kwargs[k].astype(np.float64)
                if np.all(np.isnan(vec2xr_kwargs[k])):
                    vec2xr_kwargs[k] = None
            else:
                vec2xr_kwargs[k] = np.array(v)
            if vec2xr_kwargs[k] is not None:
                if k == "Gmass":
                    vec2xr_kwargs[k] *= self.M2DU**3 / self.S2TU**2
                if k == "mass":
                    vec2xr_kwargs[k] *= self.KG2MU
                if k == "rh" or k == "radius" or k == "rhill":
                    vec2xr_kwargs[k] *= self.M2DU
                if k == "vh":
                    vec2xr_kwargs[k] *= self.M2DU / self.S2TU
                if k == "rot":
                    vec2xr_kwargs[k] /= self.S2TU
                if k == "j2rp2":
                    vec2xr_kwargs[k] *= self.M2DU**2
                if k == "j4rp4":
                    vec2xr_kwargs[k] *= self.M2DU**4

        vec2xr_kwargs["time"] = np.array([self.param["TSTART"]])

        # Create a Dataset containing the new bodies
        dsnew = self._vec2xr(**vec2xr_kwargs)
        dsnew = self._set_id_number(dsnew)
        dsnew = self._set_particle_type(dsnew)
        if "particle_type" in self.data:
            if CB_TYPE_NAME in self.data["particle_type"] and CB_TYPE_NAME in dsnew["particle_type"]:
                self.data = self._set_particle_type(
                    self.data
                )  # Make sure we update the original dataset if there is going to be a central body change

        if CB_TYPE_NAME in dsnew["particle_type"]:
            cbname = dsnew["name"].where(dsnew["particle_type"] == CB_TYPE_NAME, drop=True).values[0]
            GMcb = dsnew["Gmass"].sel(name=cbname)
        elif CB_TYPE_NAME in self.data.particle_type:
            cbname = self.data["name"].where(self.data["particle_type"] == CB_TYPE_NAME, drop=True).values[0]
            GMcb = self.data["Gmass"].sel(name=cbname)
        else:
            raise ValueError("No central body found in either the new dataset or the existing dataset")

        if (
            align_to_central_body_rotation and cbname not in dsnew["name"]
        ):  # If a new central body is being added, then the rotation occurs after the two datasets are merged
            rot = init_cond.get_solar_system_body_mass_rotation(cbname, ephemerides_start_date=date, verbose=verbose)["rot"]
            if rot is not None:
                rot *= self.param["TU2S"]
                dsnew = dsnew.rotate(pole=rot)

        dsnew = dsnew.xv2el(GMcb)

        dsnew = self._combine_and_fix_dsnew(dsnew, align_to_central_body_rotation, verbose=verbose, **kwargs)
        self.save(verbose=verbose)

        return

    def get_ephemeris_date(self, arg_list: str | list[str] | None = None, **kwargs: Any) -> str:
        """
        Prints the current value of the ephemeris date.

        Parameters
        ----------
        arg_list : str | List[str], optional
            A single string or list of strings containing the names of the features to extract. Default is all of:
            ["integrator"]
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        ephemeris_date : str
            The ISO-formatted date string for the ephemeris computation
        """
        verbose = kwargs.pop("verbose", self.verbose)

        try:
            self.ephemeris_date
        except:
            if verbose:
                warnings.warn("ephemeris_date is not set", stacklevel=2)
            return

        valid_arg = {"ephemeris_date": self.ephemeris_date}

        ephemeris_date = self._get_instance_var(arg_list, valid_arg, **kwargs)

        return ephemeris_date

    def _get_instance_var(
        self,
        arg_list: str | list[str] | None = None,
        valid_arg: dict | None = None,
        **kwargs: Any,
    ) -> tuple[Any, ...]:
        """
        Prints the current value of an instance variable.

        Parameters
        ----------
        arg_list : str | List[str]
            A single string or list of strings containing the names of the the instance variable to get.
        valid_arg : dict
            A dictionary where the key is the parameter argument and the value is the equivalent instance variable value.
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            get_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.

        Returns
        -------
        Tuple[Any, ...]
            Instance variable values given by the arg_list
        """
        verbose = kwargs.pop("verbose", self.verbose)

        arg_vals = []
        if verbose:
            if arg_list is None:
                arg_list = list(valid_arg.keys())
            for arg in arg_list:
                if arg in valid_arg:
                    print(f"{arg:<{self._getter_column_width}} {valid_arg[arg]}")
                    arg_vals.append(valid_arg[arg])

        return tuple(arg_vals)

    def _validate_body_arguments(
        self,
        name: str | list[str] | npt.NDArray[np.str_] | None = None,
        id: int | list[int] | npt.NDArray[np.int_] | None = None,
        a: float | list[float] | npt.NDArray[np.float_] | None = None,
        e: float | list[float] | npt.NDArray[np.float_] | None = None,
        inc: float | list[float] | npt.NDArray[np.float_] | None = None,
        capom: float | list[float] | npt.NDArray[np.float_] | None = None,
        omega: float | list[float] | npt.NDArray[np.float_] | None = None,
        capm: float | list[float] | npt.NDArray[np.float_] | None = None,
        rh: list[float] | list[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
        vh: list[float] | list[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
        mass: float | list[float] | npt.NDArray[np.float_] | None = None,
        Gmass: float | list[float] | npt.NDArray[np.float_] | None = None,
        radius: float | list[float] | npt.NDArray[np.float_] | None = None,
        rhill: float | list[float] | npt.NDArray[np.float_] | None = None,
        rot: list[float] | list[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
        Ip: list[float] | npt.NDArray[np.float_] | None = None,
        rotphase: float | list[float] | npt.NDArray[np.float_] | None = None,
        j2rp2: float | list[float] | npt.NDArray[np.float_] | None = None,
        j4rp4: float | list[float] | npt.NDArray[np.float_] | None = None,
        c_lm: list[float] | list[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Validates and formats the input the add_body and modify_body methods.

        Parameters
        ----------
        name : str or array-like of str, optional
            Name or names of Bodies.
        id : int or array-like of int, optional
            ID or IDs of Bodies.
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
        rot : (3) or (n,3) array-like of float, optional
            Rotation rate vectors if these are massive bodies with rotation enabled.
        Ip : (3) or (n,3) array-like of float, optional
            Principal axes moments of inertia vectors if these are massive bodies with rotation enabled.
        rotphase : float, optional
            rotation phase angle in degreesif these are massive bodies with rotation enabled
        j2rp2 : float, optional
            Non-normalized J2 values (e.g. J2*R**2, where R is the body radius) if this is a central body (only one of J2 or c_lm can be passed)
        j4rp4 : float, optional
            Non-normalized J4 values (e.g. J4*R**4, where R is the body radius) if this is a central body (only one of J4 or c_lm can be passed)
        c_lm : (2,l_max+1,l_max+1) array-like of float, optional
            Spherical harmonics coefficients if this is a central body (only one of J2/J4 or c_lm can be passed)
        align_to_central_body_rotation : bool, default False
            If True, the cartesian coordinates will be aligned to the rotation pole of the central body. This is only valid for when
            rotation is enabled.

        Returns
        -------
        None
            Sets the data and init_cond instance variables each with a SwiftestDataset containing the body or bodies that were added
        """

        # convert all inputs to numpy arrays
        def input_to_array(val, t, n=None):
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
                val = np.array([val], dtype=t)
            else:
                try:
                    val = np.array(val, dtype=t)
                except:
                    raise ValueError(f"{val} cannot be converted to a numpy array")

            if n is None:
                return val, len(val)
            else:
                if n != len(val):
                    raise ValueError(f"Mismatched array lengths in add_body. Got {len(val)} when expecting {n}")
                return val, n

        def input_to_array_3d(val, n=None):
            if val is None:
                return None, n
            else:
                try:
                    val = np.array(val, dtype=np.float64)
                except:
                    raise ValueError(f"{val} cannot be converted to a numpy array")
                if n is None:
                    ndims = len(val.shape)
                    if ndims > 2 or ndims == 0:
                        raise ValueError(f"Argument must be an (n,3) or (3,) array. This one is {val.shape}")
                    else:
                        if val.shape[-1] != 3:
                            raise ValueError(f"Argument must be a 3-dimensional vector. This one has {val.shape[0]}!")
                        if val.size == 3:
                            n = 1
                        else:
                            n = val.shape[0]
                if n == 1:
                    if val.shape != (1, 3) and val.shape != (3,):
                        raise ValueError(f"Argument is an incorrect shape. Expected {(n, 3)} or {(3, 1)}. Got {val.shape} instead")
                    elif val.shape == (3,):
                        val = np.expand_dims(val, axis=0)
                elif val.shape != (n, 3) and val.shape != (3, n):
                    raise ValueError(f"Argument is an incorrect shape. Expected {(n, 3)} or {(3, n)}. Got {val.shape} instead")
                elif val.shape == (3, n):
                    val = val.T

            return val, n

        def input_to_clm_array(val, n):
            # Create function to convert c_lm array to numpy array
            if val is None:
                return None, n
            elif isinstance(val, np.ndarray):
                pass
            else:
                try:
                    val = np.array(val, dtype=np.float64)
                except:
                    raise ValueError(f"{val} cannot be converted to a numpy array")
                ndims = len(val.shape)
                if ndims != 3 or val.shape[0] != 2 or val.shape[1] != val.shape[2]:
                    raise ValueError(f"C_lm is an incorrect shape. Expected (2, l_max + 1, l_max + 1). got {val.shape} instead.")
            return val, n

        nbodies = None
        name, nbodies = input_to_array(name, "s", nbodies)
        id, nbodies = input_to_array(id, "i", nbodies)
        a, nbodies = input_to_array(a, "f", nbodies)
        e, nbodies = input_to_array(e, "f", nbodies)
        inc, nbodies = input_to_array(inc, "f", nbodies)
        capom, nbodies = input_to_array(capom, "f", nbodies)
        omega, nbodies = input_to_array(omega, "f", nbodies)
        capm, nbodies = input_to_array(capm, "f", nbodies)
        mass, nbodies = input_to_array(mass, "f", nbodies)
        Gmass, nbodies = input_to_array(Gmass, "f", nbodies)
        rhill, nbodies = input_to_array(rhill, "f", nbodies)
        radius, nbodies = input_to_array(radius, "f", nbodies)

        rh, nbodies = input_to_array_3d(rh, nbodies)
        vh, nbodies = input_to_array_3d(vh, nbodies)
        rot, nbodies = input_to_array_3d(rot, nbodies)
        Ip, nbodies = input_to_array_3d(Ip, nbodies)
        rotphase, nbodies = input_to_array(rotphase, "f", nbodies)

        j2rp2, nbodies = input_to_array(j2rp2, "f", nbodies)
        j4rp4, nbodies = input_to_array(j4rp4, "f", nbodies)
        c_lm, nbodies = input_to_clm_array(c_lm, nbodies)

        if mass is not None and Gmass is not None:
            raise ValueError("Cannot use mass and Gmass inputs simultaneously!")

        if rh is not None or vh is not None:
            if a is not None or e is not None or inc is not None or capom is not None or omega is not None or capm is not None:
                raise ValueError("Only cartesian values or orbital elements may be passed, but not both.")

        if j2rp2 is not None or j4rp4 is not None:
            if c_lm is not None:
                raise ValueError("Cannot use J2/J4 and c_lm inputs simultaneously!")
        if a is not None:
            a = np.abs(a)  # Ensure that the semimajor axis is positive, per Swiftest convention for hyperbolic orbits

        validated_arguments = locals().copy()
        validated_arguments.pop("self")

        return validated_arguments

    def add_body(
        self,
        name: str | list[str] | npt.NDArray[np.str_] | None = None,
        a: float | list[float] | npt.NDArray[np.float_] | None = None,
        e: float | list[float] | npt.NDArray[np.float_] | None = None,
        inc: float | list[float] | npt.NDArray[np.float_] | None = None,
        capom: float | list[float] | npt.NDArray[np.float_] | None = None,
        omega: float | list[float] | npt.NDArray[np.float_] | None = None,
        capm: float | list[float] | npt.NDArray[np.float_] | None = None,
        rh: list[float] | list[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
        vh: list[float] | list[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
        mass: float | list[float] | npt.NDArray[np.float_] | None = None,
        Gmass: float | list[float] | npt.NDArray[np.float_] | None = None,
        radius: float | list[float] | npt.NDArray[np.float_] | None = None,
        rhill: float | list[float] | npt.NDArray[np.float_] | None = None,
        rot: list[float] | list[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
        Ip: list[float] | npt.NDArray[np.float_] | None = None,
        rotphase: float | list[float] | npt.NDArray[np.float_] | None = None,
        j2rp2: float | list[float] | npt.NDArray[np.float_] | None = None,
        j4rp4: float | list[float] | npt.NDArray[np.float_] | None = None,
        c_lm: list[float] | list[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
        align_to_central_body_rotation: bool = False,
        **kwargs: Any,
    ) -> None:
        """
        Adds a body (test particle or massive body) to the internal Dataset given a set of either orbital elements or cartesian state vectors.

        If orbital elements are passed, cartesian state vectors are computed and vice versa, using the currently-assigned central body, so cannot both be passed. Input all angles in degrees and dimensional quantities in the unit system defined in the current Simulation instance.

        This method will update the data attribute with the new body or bodies added to the existing Dataset.

        Parameters
        ----------
        name : str or array-like of str, optional
            Name or names of Bodies. If none passed, name will be "Body{id}"
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
        rot : (3) or (n,3) array-like of float, optional
            Rotation rate vectors if these are massive bodies with rotation enabled.
        Ip : (3) or (n,3) array-like of float, optional
            Principal axes moments of inertia vectors if these are massive bodies with rotation enabled.
        rotphase : float, optional
            rotation phase angle in degreesif these are massive bodies with rotation enabled
        j2rp2 : float, optional
            Non-normalized J2 values (e.g. J2*R**2, where R is the body radius) if this is a central body (only one of J2 or c_lm can be passed)
        j4rp4 : float, optional
            Non-normalized J4 values (e.g. J4*R**4, where R is the body radius) if this is a central body (only one of J4 or c_lm can be passed)
        c_lm : (2,l_max+1,l_max+1) array-like of float, optional
            Spherical harmonics coefficients if this is a central body (only one of J2/J4 or c_lm can be passed)
        align_to_central_body_rotation : bool, default False
            If True, the cartesian coordinates will be aligned to the rotation pole of the central body. This is only valid for when
            rotation is enabled.

        Returns
        -------
        None
            Sets the data and init_cond instance variables each with a SwiftestDataset containing the body or bodies that were added
        """
        from .constants import CB_TYPE_NAME

        verbose = kwargs.pop("verbose", self.verbose)

        # This allows us to re-use the same validation function for both add_body and modify_body
        arguments = locals().copy()
        arguments.pop("self")
        arguments = self._validate_body_arguments(**arguments)
        name = arguments["name"]
        a = arguments["a"]
        e = arguments["e"]
        inc = arguments["inc"]
        capom = arguments["capom"]
        omega = arguments["omega"]
        capm = arguments["capm"]
        rh = arguments["rh"]
        vh = arguments["vh"]
        mass = arguments["mass"]
        Gmass = arguments["Gmass"]
        radius = arguments["radius"]
        rhill = arguments["rhill"]
        rot = arguments["rot"]
        Ip = arguments["Ip"]
        rotphase = arguments["rotphase"]
        j2rp2 = arguments["j2rp2"]
        j4rp4 = arguments["j4rp4"]
        c_lm = arguments["c_lm"]
        nbodies = arguments["nbodies"]

        # Adding new bodies imposes additional constraints on arguments that are not present when modifying existing bodies
        if rh is not None and vh is None:
            raise ValueError("If rh is passed, vh must also be passed")
        if vh is not None and rh is None:
            raise ValueError("If vh is passed, rh must also be passed")
        if rh is None:
            if a is None:
                raise ValueError("Orbital element input requires at least a value for a (semimajor axis)")
            if e is None:
                e = np.zeros_like(a)
            if inc is None:
                inc = np.zeros_like(a)
            if capom is None:
                capom = np.zeros_like(a)
            if omega is None:
                omega = np.zeros_like(a)
            if capm is None:
                capm = np.zeros_like(a)

        if Gmass is not None:
            mass = Gmass / self.GU
        if mass is not None:
            Gmass = mass * self.GU

        if name is None or len(name) == 0:
            if len(self.data) == 0:
                maxid = -1
            else:
                maxid = self.data.id.max().values[()]
            id = np.arange(start=maxid + 1, stop=maxid + 1 + nbodies, dtype=int)
            name = np.char.mod("Body%d", id)

        time = [self.param["TSTART"]]

        if "name" in self.data:
            bad_names = [n for n in name if n in self.data.name.values]
            if len(bad_names) > 0:
                name = [n for n in name if n not in bad_names]
                bad_names = ", ".join(bad_names)
                if verbose:
                    warnings.warn(
                        f"The following names are already in use and will not be added: {bad_names}",
                        stacklevel=2,
                    )
        if len(name) == 0:
            if verbose:
                print("No valid names found")
            return
        if verbose:
            name_str = ", ".join(name)
            print(f"Adding bodies: {name_str}")
        dsnew = self._vec2xr(
            name=name,
            a=a,
            e=e,
            inc=inc,
            capom=capom,
            omega=omega,
            capm=capm,
            Gmass=Gmass,
            mass=mass,
            radius=radius,
            rhill=rhill,
            Ip=Ip,
            rh=rh,
            vh=vh,
            rot=rot,
            j2rp2=j2rp2,
            j4rp4=j4rp4,
            c_lm=c_lm,
            rotphase=rotphase,
            time=time,
        )

        dsnew = self._set_id_number(dsnew)
        dsnew = self._set_particle_type(dsnew)
        if "particle_type" in self.data:
            if CB_TYPE_NAME in self.data["particle_type"] and CB_TYPE_NAME in dsnew["particle_type"]:
                self.data = self._set_particle_type(
                    self.data
                )  # Make sure we update the original dataset if there is going to be a central body change

        if CB_TYPE_NAME in dsnew["particle_type"]:
            cbname = dsnew["name"].where(dsnew["particle_type"] == CB_TYPE_NAME, drop=True).values[0]
            GMcb = dsnew["Gmass"].sel(name=cbname)
        elif CB_TYPE_NAME in self.data.particle_type:
            cbname = self.data["name"].where(self.data["particle_type"] == CB_TYPE_NAME, drop=True).values[0]
            GMcb = self.data["Gmass"].sel(name=cbname)
        else:
            raise ValueError("No central body found in either the old or new Dataset")

        if rh is None or vh is None:
            dsnew = dsnew.el2xv(GMcb)
        if a is None:
            dsnew = dsnew.xv2el(GMcb)

        dsnew = self._combine_and_fix_dsnew(dsnew, align_to_central_body_rotation, verbose=verbose, **kwargs)
        self.save(verbose=verbose)

        return

    def _vec2xr(
        self,
        name: str | npt.ArrayLike[str],
        id: int | npt.ArrayLike[int] | None = None,
        a: float | npt.ArrayLike[float] | None = None,
        e: float | npt.ArrayLike[float] | None = None,
        inc: float | npt.ArrayLike[float] | None = None,
        capom: float | npt.ArrayLike[float] | None = None,
        omega: float | npt.ArrayLike[float] | None = None,
        capm: float | npt.ArrayLike[float] | None = None,
        rh: npt.ArrayLike[float] | None = None,
        vh: npt.ArrayLike[float] | None = None,
        Gmass: float | npt.ArrayLike[float] | None = None,
        mass: float | npt.ArrayLike[float] | None = None,
        radius: float | npt.ArrayLike[float] | None = None,
        rhill: float | npt.ArrayLike[float] | None = None,
        rot: npt.ArrayLike[float] | None = None,
        rotphase: float | None = None,
        Ip: npt.ArrayLike[float] | None = None,
        j2rp2: float | npt.ArrayLike[float] | None = None,
        j4rp4: float | npt.ArrayLike[float] | None = None,
        c_lm: npt.ArrayLike[float] | None = None,
        time: npt.ArrayLike[float] | None = None,
        **kwargs: Any,
    ) -> SwiftestDataset:
        """
        Converts and stores the variables of all bodies in an xarray dataset.

        Parameters
        ----------
        param : dict
            Swiftest simulation parameters.
        name : str or array-like of str
            Name or names of bodies. Bodies are indexed by name, so these must be unique
        id : int or array-like of int, optional
            Unique id values.
        a : float or array-like of float, optional
            semimajor axis for param['IN_FORM'] == "EL"
        e : float or array-like of float, optional
            eccentricity  for param['IN_FORM'] == "EL"
        inc : float or array-like of float, optional
            inclination for param['IN_FORM'] == "EL"
        capom : float or array-like of float, optional
            longitude of periapsis for param['IN_FORM'] == "EL"
        omega : float or array-like of float, optional
            argument of periapsis for param['IN_FORM'] == "EL"
        capm : float or array-like of float, optional
            mean anomaly for param['IN_FORM'] == "EL"
        rh : (n,3) array-like of float, optional
            Position vector array.
        vh : (n,3) array-like of float, optional
            Velocity vector array.
        Gmass : float or array-like of float, optional
            G*mass values if these are massive bodies. If mass is passed, but Gmass is not, then Gmass is calculated from mass.
        mass : float or array-like of float, optional
            Mass values if these are massive bodies. If Gmass is passed, then mass is calculated from Gmass, regardless of whether mass is passed.
        radius : float or array-like of float, optional
            Radius values if these are massive bodies
        rhill : float or array-like of float, optional
            Hill's radius values if these are massive bodies
        rot:  (n,3) array-like of float, optional
            Rotation rate vectors if these are massive bodies with rotation enabled in deg/TU
        rotphase : float
            rotational phase angle of the central body in degrees
        Ip: (n,3) array-like of flaot, optional
            Principal axes moments of inertia vectors if these are massive bodies with rotation enabled. This can be used
            instead of passing Ip1, Ip2, and Ip3 separately
        j2rp2 : float or array-like of float, optional
            J_2R^2 value for the body
        j4rp4 : float or array-like of float, optional
            J_4R^4 value for the body
        c_lm : (2, lmax + 1, lmax + 1) array of floats, optional
            Spherical Harmonics coefficients; lmax = max spherical harmonics order
        time : array of floats
            Time at start of simulation
        **kwargs : Any
            Additional keyword arguments. These are ignored.

        Returns
        -------
        ds : SwiftestDataset
            Dataset containing the variables of all bodies passed in kwargs
        """
        # Validate the inputs
        if name is None or len(name) == 0:
            raise ValueError("Name must be passed")

        if isinstance(name, str):
            nbody = 1
        else:
            nbody = len(name)

        scalar_dims = ["name"]
        vector_dims = ["name", "space"]
        space_coords = np.array(["x", "y", "z"])

        vector_vars = ["rh", "vh", "Ip", "rot"]
        scalar_vars = [
            "id",
            "a",
            "e",
            "inc",
            "capom",
            "omega",
            "capm",
            "mass",
            "Gmass",
            "radius",
            "rhill",
            "j2rp2",
            "j4rp4",
            "rotphase",
        ]
        sph_vars = ["c_lm"]
        time_vars = [
            "status",
            "rh",
            "vh",
            "Ip",
            "rot",
            "a",
            "e",
            "inc",
            "capom",
            "omega",
            "capm",
            "mass",
            "Gmass",
            "radius",
            "rhill",
            "j2rp2",
            "j4rp4",
            "rotphase",
            "c_lm",
        ]

        if "ROTATION" in self.param and self.param["ROTATION"] == True:
            if rot is None and Gmass is not None:
                rot = np.zeros((nbody, 3))
            if Ip is None and Gmass is not None:
                Ip = np.full((nbody, 3), 0.4)
            if rotphase is None and Gmass is not None:
                rotphase = np.zeros(nbody)

        if time is None:
            if "time" in self.data:
                time = self.data["time"].values[-1:]
            else:
                time = np.array([0.0])

        if self.param["CHK_CLOSE"]:
            if Gmass is not None and radius is None:
                raise ValueError("If Gmass is passed, then radius must also be passed when CHK_CLOSE is True")

        if Gmass is not None:
            Gmass[np.isnan(Gmass)] = 0.0
            mass = Gmass / self.GU
        elif mass is not None:
            mass[np.isnan(mass)] = 0.0
            Gmass = mass * self.GU
        else:
            mass = np.zeros(nbody)
            Gmass = np.zeros(nbody)

        valid_vars = vector_vars + scalar_vars + sph_vars + ["time", "id"]

        input_vars = {k: v for k, v in locals().items() if k in valid_vars and v is not None}

        data_vars = {k: (scalar_dims, v) for k, v in input_vars.items() if k in scalar_vars}
        data_vars.update({k: (vector_dims, v) for k, v in input_vars.items() if k in vector_vars})
        ds = xr.Dataset(
            data_vars=data_vars,
            coords={
                "name": (["name"], name),
                "space": (["space"], space_coords),
            },
        )
        # create a C_lm Dataset and combine
        if c_lm is not None:
            clm_xr = xr.DataArray(
                data=c_lm,
                coords={
                    "sign": (["sign"], [1, -1]),
                    "l": (["l"], range(0, c_lm.shape[1])),
                    "m": (["m"], range(0, c_lm.shape[2])),
                },
            ).to_dataset(name="c_lm")
            clm_xr = clm_xr.expand_dims(dim={"name": nbody}, axis=0)

            ds = xr.combine_by_coords([ds, clm_xr])

        time_vars = [v for v in time_vars if v in ds]
        for v in time_vars:
            ds[v] = ds[v].expand_dims(dim={"time": 1}, axis=0).assign_coords({"time": time})

        return SwiftestDataset(ds)

    def _set_id_number(self, ds: SwiftestDataset) -> SwiftestDataset:
        """
        Sets the id numbers for new bodies to be added to the Dataset.

        It will set the most massive body of both the old and new Dataset to have id=0 to indicate that it is to be considered the central body.

        Parameters
        ----------
        ds : SwiftestDataset
            Dataset to evaluate

        Returns
        -------
        ds : SwiftestDataset
            Dataset with updated id values.

        Notes
        -----
        If a body from the new Dataset is found to be more massive than one from the existing Dataset, then the id of the old central
        body will be modified to no longer be id==0.
        """
        # Make sure neither the old nor the new Dataset is indexed by id, as these will shift
        if "id" in ds.dims:
            ds = ds.swap_dims({"id": "name"})
        if "id" in self.data.dims:
            self.data = self.data.swap_dims({"id": "name"})

        nnew = ds.name.size
        # Increment id numbers
        if len(self.data) == 0 or self.data.id.sizes["name"] == 0:
            maxid = 0
        else:
            maxid = self.data.id.max().values[()]

        # Find out which body will be the central body (i.e. the most massive)
        if "Gmass" in ds:
            new_bigidx = ds["Gmass"].argmax(dim="name")
            new_bigname = ds["name"].isel(name=new_bigidx).values[()]
            new_Gmass = ds["Gmass"].isel(name=new_bigidx).values[0]
        else:
            new_bigname = None

        if "Gmass" in self.data and self.data.Gmass.sizes["name"] > 0:
            old_bigidx = self.data["Gmass"].argmax(dim="name")
            old_bigname = self.data["name"].isel(name=old_bigidx).values[()]
            old_Gmass = self.data["Gmass"].isel(name=old_bigidx).values[0]
        else:
            old_bigname = None

        if new_bigname is None and old_bigname is None:
            raise ValueError("No central body found in either the old or new Dataset")

        # Establish a new id variable for the new Dataset
        ds["id"] = xr.DataArray(np.arange(start=maxid + 1, stop=maxid + 1 + nnew, dtype=int), dims="name")
        if old_bigname is None:
            oldid = ds["id"].sel(name=new_bigname).values[()]
            ds["id"].loc[{"name": new_bigname}] = 0
            # Ensure we don't have a gap:
            ds["id"] = xr.where(ds["id"] > oldid, ds["id"] - 1, ds["id"])
        elif new_bigname is not None:
            if new_Gmass > old_Gmass:
                oldid = ds["id"].sel(name=new_bigname).values[()]
                ds["id"].loc[{"name": new_bigname}] = 0
                self.data["id"].loc[{"name": old_bigname}] = oldid

        return SwiftestDataset(ds)

    def _set_particle_type(self, ds: SwiftestDataset) -> SwiftestDataset:
        """
        Sets the particle type based on the values of Gmass.

        Parameters
        ----------
        ds : SwiftestDataset
            Dataset to evaluate

        Returns
        -------
        ds : SwiftestDataset
            Dataset with updated particle type values.
        """
        from .constants import (
            CB_TYPE_NAME,
            PL_TINY_TYPE_NAME,
            PL_TYPE_NAME,
            TP_TYPE_NAME,
        )

        if "Gmass" in ds:
            Gmass = ds.isel(time=0).Gmass
            ds["particle_type"] = xr.where(
                ds["id"] == 0,
                CB_TYPE_NAME,
                xr.where(np.isnan(Gmass) | (Gmass == 0.0), TP_TYPE_NAME, PL_TYPE_NAME),
            )
            if self.integrator == "symba" and "GMTINY" in self.param and self.param["GMTINY"] is not None:
                ds["particle_type"] = xr.where(
                    (ds["particle_type"] == PL_TYPE_NAME) & (Gmass < self.param["GMTINY"]),
                    PL_TINY_TYPE_NAME,
                    ds["particle_type"],
                )
        else:
            ds["particle_type"] = xr.full_like(ds["name"], TP_TYPE_NAME, dtype="<U32")
        return ds

    def _get_nvals(self, ds: SwiftestDataset) -> SwiftestDataset:
        """
        Computes the values of ntp, npl, and nplm.

        Parameters
        ----------
        ds : SwiftestDataset
            Dataset to evaluate

        Returns
        -------
        ds : SwiftestDataset
            Dataset with updated values of ntp, npl, and nplm values.
        """
        if "name" in ds.dims:
            count_dim = "name"
        elif "id" in ds.dims:
            count_dim = "id"

        if "particle_type" not in ds:
            ds = self._set_particle_type(ds)

        if "Gmass" in ds:
            Gmass = ds.Gmass
            ds["ntp"] = Gmass.where((ds.id != 0) & (~np.isnan(Gmass) & (Gmass == 0.0))).count(dim=[count_dim])
            ds["npl"] = Gmass.where((ds.id != 0) & (~np.isnan(Gmass) & (Gmass > 0.0))).count(dim=[count_dim])
            if self.integrator == "symba" and "GMTINY" in self.param and self.param["GMTINY"] is not None:
                ds["nplm"] = Gmass.where((ds.id != 0) & (~np.isnan(Gmass) & (Gmass > self.param["GMTINY"]))).count(dim=[count_dim])
        else:
            ds["ntp"] = ds.id.where(ds.id != 0).count(dim=[count_dim])
            ds["npl"] = xr.zeros_like(ds["ntp"])
            if self.integrator == "symba" and "GMTINY" in self.param and self.param["GMTINY"] is not None:
                ds["nplm"] = xr.zeros_like(ds["ntp"])
        return ds

    def _get_valid_body_list(
        self,
        name: str | list[str] | npt.NDArray[np.str_] | None = None,
        id: int | list[int] | npt.NDArray[np.int_] | None = None,
    ):
        """
        Returns a list of valid body names in the dataset given the input name or id.

        Parameters
        ----------
        name : str or array-like of str, optional
            Name or names of bodies to select
        id : int or array-like of int, optional
            Id value or values of bodies to select

        Returns
        -------
        name : list of str
            List of valid body names
        """
        if name is None and id is None:
            raise ValueError("You must pass either name or id as arguments")
        if name is not None and id is not None:
            raise ValueError("You must pass only name or id as arguments, but not both")

        if name is not None:
            if type(name) is str or type(name) is int:
                name = [name]
            invalid_names = ", ".join([n for n in name if n not in self.data.name.values])
            if len(invalid_names) > 0:
                warnings.warn(f"{invalid_names} not found in the Dataset. remove_body is ignoring these names.")
        else:
            if type(id) is int or type(id) is name:
                id = [id]
            invalid_ids = ", ".join([f"{i}" for i in id if i not in self.data.id.values])
            if len(invalid_ids) > 0:
                warnings.warn(f"id number(s) {invalid_ids} not found in the Dataset. remove_body is ignoring these ids.")
            name = self.data.name.where(self.data.id == id, drop=True).values.tolist()

        return name

    def remove_body(
        self,
        name: str | list[str] | npt.NDArray[np.str_] | None = None,
        id: int | list[int] | npt.NDArray[np.int_] | None = None,
    ):
        """
        Removes a body (test particle or massive body) from the internal Dataset.

        This method will update `data` and `init_cond` attributes with the body or bodies removed.

        Only `name` or `id` may be passed, but not both.

        Parameters
        ----------
        name : str or array-like of str, optional
            Name or names of bodies to remove.
        id : int or array-like of int, optional
            Id value or values of bodies to removed.

        Returns
        -------
        None
        """
        names = self._get_valid_body_list(name=name, id=id)
        keepnames = [n for n in self.data.name.values if n not in names]
        if len(keepnames) == 0:
            warnings.warn("No bodies left in the Dataset after remove_body")
        if len(keepnames) == len(self.data.name):
            warnings.warn("No bodies found that can be removed from the Dataset.")
            return

        self.data = self.data.sel(name=keepnames)

        self.data = self._get_nvals(self.data)

        return

    def modify_body(
        self,
        name: str | list[str] | npt.NDArray[np.str_] | None = None,
        id: int | list[int] | npt.NDArray[np.int_] | None = None,
        a: float | list[float] | npt.NDArray[np.float_] | None = None,
        e: float | list[float] | npt.NDArray[np.float_] | None = None,
        inc: float | list[float] | npt.NDArray[np.float_] | None = None,
        capom: float | list[float] | npt.NDArray[np.float_] | None = None,
        omega: float | list[float] | npt.NDArray[np.float_] | None = None,
        capm: float | list[float] | npt.NDArray[np.float_] | None = None,
        rh: list[float] | list[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
        vh: list[float] | list[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
        mass: float | list[float] | npt.NDArray[np.float_] | None = None,
        Gmass: float | list[float] | npt.NDArray[np.float_] | None = None,
        radius: float | list[float] | npt.NDArray[np.float_] | None = None,
        rhill: float | list[float] | npt.NDArray[np.float_] | None = None,
        rot: list[float] | list[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
        Ip: list[float] | npt.NDArray[np.float_] | None = None,
        rotphase: float | list[float] | npt.NDArray[np.float_] | None = None,
        j2rp2: float | list[float] | npt.NDArray[np.float_] | None = None,
        j4rp4: float | list[float] | npt.NDArray[np.float_] | None = None,
        c_lm: list[float] | list[npt.NDArray[np.float_]] | npt.NDArray[np.float_] | None = None,
        align_to_central_body_rotation: bool = False,
        framenum: int = -1,
        **kwargs: Any,
    ) -> None:
        """
        Modifies an existing body in the internal Dataset given a new value of either the orbital elements or cartesian state vectors, or the physical property of the body (mass, radius, etc).

        Input all angles in degrees and dimensions in the units defined in the current Simulation instance. Currently, this will only modify the last entry of the body in the time dimension.  This method will update the data attribute with the modified body or bodies added to the existing Dataset.

        Parameters
        ----------
        name : str or array-like of str, optional
            Name or names of bodies to modify.
        id : int or array-like of int, optional
            Id value or values of bodies to modify. Only one of name or id may be passed, but not both.
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
        rot : (3) or (n,3) array-like of float, optional
            Rotation rate vectors if these are massive bodies with rotation enabled.
        Ip : (3) or (n,3) array-like of float, optional
            Principal axes moments of inertia vectors if these are massive bodies with rotation enabled.
        rotphase : float, optional
            rotation phase angle in degreesif these are massive bodies with rotation enabled
        j2rp2 : float, optional
            Non-normalized J2 values (e.g. J2*R**2, where R is the body radius) if this is a central body (only one of J2 or c_lm can be passed)
        j4rp4 : float, optional
            Non-normalized J4 values (e.g. J4*R**4, where R is the body radius) if this is a central body (only one of J4 or c_lm can be passed)
        c_lm : (2,l_max+1,l_max+1) array-like of float, optional
            Spherical harmonics coefficients if this is a central body (only one of J2/J4 or c_lm can be passed)
        align_to_central_body_rotation : bool, default False
            If True, the cartesian coordinates will be aligned to the rotation pole of the central body. This is only valid for when
            rotation is enabled.
        framenum : int, default -1
            Frame number to modify. If -1, the last frame is modified.

        Returns
        -------
        None
            Sets the data and init_cond instance variables each with a SwiftestDataset containing the body or bodies that were added
        """
        from .constants import CB_TYPE_NAME

        verbose = kwargs.pop("verbose", self.verbose)

        # This allows us to re-use the same validation function for both add_body and modify_body
        arguments = locals().copy()
        arguments.pop("self")
        arguments = self._validate_body_arguments(**arguments)
        name = arguments["name"]
        id = arguments["id"]
        modnames = self._get_valid_body_list(name=name, id=id)
        if modnames is None or len(modnames) == 0:
            return
        name_str = ", ".join(modnames)
        if verbose:
            print(f"Modifying bodies: {name_str}")
        dsnew = self.data.sel(name=modnames).isel(time=[framenum])

        if arguments["c_lm"] is not None:
            if "j2rp2" in dsnew:
                dsnew["j2rp2"] = xr.full_like(dsnew["j2rp2"], np.nan)
            if "j4rp4" in dsnew:
                dsnew["j4rp4"] = xr.full_like(dsnew["j4rp4"], np.nan)
        if arguments["j2rp2"] is not None or arguments["j4rp4"] is not None and "c_lm" in dsnew:
            dsnew["c_lm"] = xr.full_like(dsnew["c_lm"], np.nan)

        dsmod = self._vec2xr(**arguments)
        if "Gmass" in dsnew and arguments["Gmass"] is None and arguments["mass"] is None:
            dsmod["Gmass"] = dsnew["Gmass"]
            dsmod["mass"] = dsnew["mass"]
        dsnew.update(dsmod)
        if arguments["mass"] is not None or arguments["Gmass"] is not None:
            dsnew = self._set_particle_type(dsnew)
            if "particle_type" in self.data:
                if CB_TYPE_NAME in self.data["particle_type"] and CB_TYPE_NAME in dsnew["particle_type"]:
                    self.data = self._set_particle_type(
                        self.data
                    )  # Make sure we update the original dataset if there is going to be a central body change

            if CB_TYPE_NAME in dsnew["particle_type"]:
                cbname = dsnew["name"].where(dsnew["particle_type"] == CB_TYPE_NAME, drop=True).values[0]
                GMcb = dsnew["Gmass"].sel(name=cbname).isel(time=framenum)
            elif CB_TYPE_NAME in self.data.particle_type:
                cbname = self.data["name"].where(self.data["particle_type"] == CB_TYPE_NAME, drop=True).values[0]
                GMcb = self.data["Gmass"].sel(name=cbname).isel(time=framenum)
            else:
                raise ValueError("No central body found in either the old or new Dataset")
        else:
            cbname = self.data["name"].where(self.data["particle_type"] == CB_TYPE_NAME, drop=True).values[0]
            GMcb = self.data["Gmass"].sel(name=cbname).isel(time=framenum)

        if any(arguments[var] is not None for var in ["rh", "vh"]):
            dsnew = dsnew.xv2el(GMcb)
        if any(arguments[var] is not None for var in ["a", "e", "inc", "capom", "omega", "capm"]):
            dsnew = dsnew.el2xv(GMcb)
        if any(arguments[var] is not None for var in ["mass", "Gmass"]):
            dsnew = dsnew.xv2el(GMcb)

        dsnew = self._combine_and_fix_dsnew(dsnew, align_to_central_body_rotation, verbose=verbose, **kwargs)
        self.save(verbose=verbose)

        return

    def _combine_and_fix_dsnew(
        self,
        dsnew: SwiftestDataset,
        align_to_central_body_rotation: bool = False,
        **kwargs: Any,
    ) -> SwiftestDataset:
        """
        Combines the new Dataset with the old one. Also computes the values of ntp and npl and sets the proper types.

        Parameters
        ----------
        dsnew : SwiftestDataset
            Dataset with new bodies
        align_to_central_body_rotation : bool, default False
            If True, the cartesian coordinates will be aligned to the rotation pole of the central body. This is only valid for when
            rotation is enabled.

        Returns
        -------
        dsnew : SwiftestDataset
            Updated Dataset with ntp, npl values and types fixed.
        """
        verbose = kwargs.pop("verbose", self.verbose)

        if not isinstance(dsnew, SwiftestDataset):
            dsnew = SwiftestDataset(dsnew)

        if verbose and "id" not in self.data.dims and len(np.unique(dsnew["name"])) != len(dsnew["name"]):
            msg = "Non-unique names detected for bodies. The Dataset will be dimensioned by integer id instead of name."
            msg += "\nConsider using unique names instead."
            print(msg)
        dsnew["status"] = xr.zeros_like(dsnew["id"]).expand_dims(dim={"time": dsnew.time.values}, axis=0)

        if "name" in self.data.coords and "name" in dsnew.coords:
            names_new = dsnew.coords["name"].to_numpy()
            names_old = self.data.coords["name"].to_numpy() if "name" in self.data.coords else np.array([])
            overlap = np.intersect1d(names_new, names_old)

            if overlap.size == 0:  # No overlapping names, so we can concat
                self.data = xr.concat([self.data, dsnew], dim="name")
            else:  # Modify the values corresponding to the overlapping names.
                for name in overlap:
                    for time in dsnew.coords["time"].values:
                        # Create a selector for the combination of name and time
                        selector_name_time = {"name": name, "time": time}
                        selector_name = {"name": name}
                        # Check if this combination exists in the old dataset
                        if ((self.data.coords["name"] == name) & (self.data.coords["time"] == time)).any():
                            # If it exists, update the corresponding values
                            for var in dsnew.data_vars:
                                if "name" in dsnew[var].dims:
                                    if var in self.data:
                                        if "name" not in self.data[var].dims:
                                            self.data[var] = (
                                                self.data[var]
                                                .expand_dims(
                                                    dim={"name": self.data.name.values},
                                                    axis=1,
                                                )
                                                .copy(deep=True)
                                            )
                                        # Check if the variable depends on both name and time
                                        if "time" in self.data[var].dims:
                                            self.data[var].loc[selector_name_time] = dsnew[var].loc[selector_name_time].values
                                        else:
                                            # Update based only on name if time is not a dimension
                                            self.data[var].loc[selector_name] = dsnew[var].loc[selector_name].values
                                    else:
                                        self.data[var] = dsnew[var]
        else:
            # If `name` is not a coord in either dataset, fall back to a safe outer-merge.
            self.data = xr.merge([self.data, dsnew], compat="override", join="outer")

        if not isinstance(self.data, SwiftestDataset):
            self.data = SwiftestDataset(self.data)

        self._set_central_body(align_to_central_body_rotation, verbose=verbose)
        dsnew = self._get_nvals(dsnew)
        self.data = self._get_nvals(self.data)
        self.data = self.data.sortby("id")
        self.data = io.reorder_dims(self.data)

        if self.param["OUT_TYPE"] == "NETCDF_DOUBLE":
            dsnew = io.fix_types(dsnew, ftype=np.float64)
            self.data = io.fix_types(self.data, ftype=np.float64)
        elif self.param["OUT_TYPE"] == "NETCDF_FLOAT":
            dsnew = io.fix_types(dsnew, ftype=np.float32)
            self.data = io.fix_types(self.data, ftype=np.float32)

        return dsnew

    def read_param(
        self,
        param_file: os.PathLike | str | None = None,
        codename: Literal["Swiftest", "Swifter", "Swift"] | None = None,
        **kwargs: Any,
    ) -> bool:
        """
        Reads in an input parameter file and stores the values in the param dictionary.

        Parameters
        ----------
        param_file : str or path-like, default is the value of the Simulation object's internal `param_file`.
            File name of the input parameter file
        codename : {"Swiftest", "Swifter", "Swift"}, default is the value of the Simulation object's internal `codename` instance variable
            Type of parameter file, either "Swift", "Swifter", or "Swiftest"

        Returns
        -------
        None
        """
        verbose = kwargs.pop("verbose", self.verbose)

        if param_file is None:
            param_file = self.simdir / self.param_file
        if codename is None:
            codename = self.codename

        if not os.path.exists(param_file):
            raise FileNotFoundError(f"Parameter file {param_file} not found.")

        if codename == "Swiftest":
            self.param = io.read_swiftest_param(param_file, self.param, verbose=verbose)
        elif codename == "Swifter":
            self.param = io.read_swifter_param(param_file, verbose=verbose)
        elif codename == "Swift":
            self.param = io.read_swift_param(param_file, verbose=verbose)
        else:
            if verbose:
                warnings.warn(
                    f'{codename} is not a recognized code name. Valid options are "Swiftest", "Swifter", or "Swift".',
                    stacklevel=2,
                )
            return

        # Set all the Simulation parameters based on values retrieved from the parameter file
        kwargs = {self._param_to_argument[key]: value for key, value in self.param.items() if key in self._param_to_argument}
        self.set_parameter(**kwargs)

        return

    def read_init_cond(
        self,
        init_cond_file_name: os.PathLike | str | None = None,
        init_cond_file_type: Literal["NETCDF", "ASCII"] | None = None,
        init_cond_format: Literal["EL", "XV"] | None = None,
        dask: bool = False,
        **kwargs: Any,
    ) -> bool:
        """
        Reads in a set of initial conditions from file and stores them as the init_cond and data instance variables.

        Parameters
        ----------
        init_cond_file_name : str, path-like, or dict, optional
            Name of the input initial condition file or files. Whether to pass a single file name or a dictionary
            depends on the argument passed to `init_cond_file_type`. If `init_cond_file_type={"NETCDF_DOUBLE","NETCDF_FLOAT"}`,
            then this will be a single file name. If `init_cond_file_type="ASCII"` then this must be a dictionary where::

                init_cond_file_name = {
                    "CB" - [path to central body initial conditions file] (Swiftest only),
                    "PL" - [path to massive body initial conditions file],
                    "TP" - [path to test particle initial conditions file] }

            If no file name is provided then the following default file names will be used.

            * "NETCDF_DOUBLE", "NETCDF_FLOAT" `init_cond_file_name = "init_cond.nc"`
            * "ASCII" `init_cond_file_name = {"CB" : "cb.in", "PL" : "pl.in", "TP" : "tp.in"}`

            Parameter input file equivalent is `NC_IN`, `CB_IN`, `PL_IN`, `TP_IN`
        init_cond_file_type : {"NETCDF_DOUBLE", "NETCDF_FLOAT", "ASCII"}, default "NETCDF_DOUBLE"
            The file type containing initial conditions for the simulation.

            * "NETCDF_DOUBLE" A single initial conditions input file in NetCDF file format of type NETCDF_DOUBLE.
            * "NETCDF_FLOAT" A single initial conditions input file in NetCDF file format of type NETCDF_FLOAT.
            * "ASCII"  Three initial conditions files in ASCII format. The individual files define the central body,

            massive body, and test particle initial conditions.
            Parameter input file equivalent is `IN_TYPE`

        init_cond_format : {"EL", "XV"}, default "XV"
            Indicates whether the input initial conditions are given as orbital elements or cartesian position and
            velocity vectors.
            If `codename` is "Swift" or "Swifter", EL initial conditions are converted to XV.
            Parameter input file equivalent is `IN_FORM`
        dask : bool, default False
            Use Dask to lazily load data (useful for very large datasets)

        Returns
        -------
        None
        """
        verbose = kwargs.pop("verbose", self.verbose)

        if init_cond_file_name is None:
            init_cond_file_name = self.param["NC_IN"]
        init_cond_file = self.simdir / self.param["NC_IN"]
        if not os.path.exists(init_cond_file):
            raise FileNotFoundError(f"Initial conditions file {init_cond_file} not found.")
        if verbose:
            print("Reading initial conditions file as .init_cond")
        if init_cond_file_type is None:
            init_cond_file_type = self.param["IN_TYPE"]
        if init_cond_format is None:
            init_cond_format = self.param["IN_FORM"]
        if "NETCDF" in init_cond_file_type:
            param_tmp = self.param.copy()
            param_tmp["BIN_OUT"] = init_cond_file
            param_tmp["IN_TYPE"] = init_cond_file_type
            param_tmp["IN_FORM"] = init_cond_format
            self.init_cond = io.swiftest2xr(param_tmp, verbose=verbose, dask=dask)
        elif verbose:
            warnings.warn(
                "Reading in ASCII initial conditions files in Python is not yet supported",
                stacklevel=2,
            )

        return

    def write_param(
        self,
        codename: Literal["Swiftest", "Swifter", "Swift"] | None = None,
        param_file: str | os.PathLike | None = None,
        param: dict | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Writes to a param.in file and determines whether the output format needs to be converted between Swift/Swifter/Swiftest.

        Parameters
        ----------
        codename : {"Swiftest", "Swifter", "Swift"}, optional
            Alternative name of the n-body code that the parameter file will be formatted for. Defaults to current instance
            variable codename
        param_file : str or path-like, optional
            Alternative file name of the input parameter file. Defaults to current instance variable self.param_file
        param : Dict, optional
            An alternative parameter dictionary to write out. Defaults to the current instance variable self.param
        **kwargs : Any
            A dictionary of additional keyword argument. These are ignored.

        Returns
        -------
        None
        """
        verbose = kwargs.pop("verbose", self.verbose)

        if codename is None:
            codename = self.codename
        if param_file is None:
            param_file = self.param_file
        if param is None:
            param = self.param

        if verbose:
            print(f"Writing parameter inputs to file {str(self.simdir / param_file)}")
        param["! VERSION"] = f"{codename} input file"

        self.simdir.mkdir(parents=True, exist_ok=True)
        # Check to see if the parameter type matches the output type. If not, we need to convert
        if codename == "Swifter" or codename == "Swiftest":
            if param["IN_TYPE"] == "ASCII":
                param.pop("NC_IN", None)
            else:
                param.pop("CB_IN", None)
                param.pop("PL_IN", None)
                param.pop("TP_IN", None)
            io.write_labeled_param(param, self.simdir / param_file)
        elif codename == "Swift":
            io.write_swift_param(param, self.simdir / param_file)
        elif verbose:
            warnings.warn(
                'Cannot process unknown code type. Call the read_param method with a valid code name. Valid options are "Swiftest", "Swifter", or "Swift".',
                stacklevel=2,
            )

        return

    def convert(
        self,
        param_file: str,
        newcodename: str = "Swiftest",
        plname: str = "pl.swiftest.in",
        tpname: str = "tp.swiftest.in",
        cbname: str = "cb.swiftest.in",
        conversion_questions: dict = {},
        dask: bool = False,
        **kwargs: Any,
    ) -> SwiftestDataset:
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
        oldparam : Dict
            The old parameter configuration.
        """
        oldparam = self.param
        if self.codename == newcodename:
            warnings.warn(
                f"This parameter configuration is already in {newcodename} format",
                stacklevel=2,
            )
            return oldparam
        if newcodename != "Swift" and newcodename != "Swifter" and newcodename != "Swiftest":
            warnings.warn(
                f'{newcodename} is an invalid code type. Valid options are "Swiftest", "Swifter", or "Swift".',
                stacklevel=2,
            )
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
            warnings.warn(
                f"Conversion from {self.codename} to {newcodename} is not supported.",
                stacklevel=2,
            )
        return oldparam

    def read_output_file(
        self,
        read_init_cond: bool = True,
        read_encounters: bool = True,
        read_collisions: bool = True,
        dask: bool = False,
        **kwargs: Any,
    ) -> None:
        """
        Reads in simulation data from an output file and stores it as a SwiftestDataset in the `data` instance variable.

        Parameters
        ----------
        read_init_cond : bool, default True
            Read in an initial conditions file along with the output file. Default is True
        dask : bool, default False
            Use Dask to lazily load data (useful for very large datasets)

        Returns
        -------
        None
            Sets the data instance variable xarray dataset
        """
        verbose = kwargs.pop("verbose", self.verbose)
        # Make a temporary copy of the parameter dictionary so we can supply the absolute path of the binary file
        # This is done to handle cases where the method is called from a different working directory than the simulation
        # results

        param_tmp = self.param.copy()
        param_tmp["BIN_OUT"] = self.simdir / self.param["BIN_OUT"]
        datafilefound = param_tmp["BIN_OUT"].exists()
        if self.codename == "Swiftest":
            if datafilefound:
                self.data = io.swiftest2xr(param_tmp, verbose=verbose, dask=dask)
                if verbose:
                    print("Swiftest simulation data stored as xarray Dataset .data")
            elif verbose:
                print(f"Output file {param_tmp['BIN_OUT']} not found.")
            if read_init_cond:
                try:
                    self.read_init_cond(dask=dask, verbose=verbose)
                    if not datafilefound:
                        self.data = self.init_cond.copy(deep=True)
                except:
                    if verbose:
                        print("No initial conditions file found. Generating from the first time step of the output file.")
                    self.init_cond = self.data.isel(time=[0]).copy(deep=True)
                    self._scrub_init_cond()

            if read_encounters:
                self.read_encounter_file(dask=dask, verbose=verbose)
            if read_collisions:
                self.read_collision_file(dask=dask, verbose=verbose)
            if verbose:
                print("Finished reading Swiftest dataset files.")

        elif self.codename == "Swifter":
            if datafilefound:
                self.data = io.swifter2xr(param_tmp, verbose=verbose)
                if verbose:
                    print("Swifter simulation data stored as SwiftestDataset .data")
            elif verbose:
                print(f"Output file {param_tmp['BIN_OUT']} not found.")
        elif self.codename == "Swift":
            if verbose:
                warnings.warn("Reading Swift simulation data is not implemented yet", stacklevel=2)
        else:
            if verbose:
                warnings.warn(
                    'Cannot process unknown code type. Call the read_param method with a valid code name. Valid options are "Swiftest", "Swifter", or "Swift".',
                    stacklevel=2,
                )
        return

    def read_encounter_file(self, dask: bool = False, **kwargs: Any) -> None:
        """
        Reads in an encounter history file and stores it as a SwiftestDataset in the `encounters` instance variable.

        Parameters
        ----------
        dask : bool, default False
            Use Dask to lazily load data (useful for very large datasets)

        Returns
        -------
        None
            Sets the encounters instance variable SwiftestDataset
        """
        verbose = kwargs.pop("verbose", self.verbose)

        enc_file = self.simdir / "encounters.nc"
        if not enc_file.exists():
            return

        if verbose:
            print("Reading encounter history file as .encounters")

        if dask:
            ds = xr.open_mfdataset(enc_file, engine="h5netcdf", mask_and_scale=False)
        else:
            with xr.open_dataset(enc_file, mask_and_scale=False) as ds:
                ds.load()
        self.encounters = io.process_netcdf_input(ds, self.param)
        ds.close()

        # Remove any overlapping time values
        tgood, tid = np.unique(self.encounters.time, return_index=True)
        self.encounters = self.encounters.isel(time=tid)
        # Remove any NaN values
        tgood = self.encounters.time.where(~np.isnan(self.encounters.time), drop=True).values
        self.encounters = self.encounters.sel(time=tgood)

        return

    def read_collision_file(self, dask: bool = False, **kwargs: Any) -> None:
        """
        Reads in a collision history file and stores it as an Xarray Dataset in the `collisions` instance variable.

        Parameters
        ----------
        dask : bool, default False
            Use Dask to lazily load data (useful for very large datasets)
        **kwargs : Any
            A dictionary of additional keyword arguments.

        Returns
        -------
        None
            Sets the collisions instance variable xarray dataset
        """
        verbose = kwargs.pop("verbose", self.verbose)

        col_file = self.simdir / "collisions.nc"
        if not col_file.exists():
            return

        if verbose:
            print("Reading collisions history file as .collisions")

        if dask:
            ds = xr.open_mfdataset(col_file, engine="h5netcdf", mask_and_scale=False, combine="nested")
        else:
            with xr.open_dataset(col_file, mask_and_scale=False) as ds:
                ds.load()

        self.collisions = io.process_netcdf_input(ds, self.param)
        ds.close()

        return

    def follow(self, codestyle: str = "Swifter", dask: bool = False, **kwargs: Any) -> SwiftestDataset:
        """
        An implementation of the Swift tool_follow algorithm. Under development. Currently only for Swift simulations.

        Parameters
        ----------
        codestyle : str, default "Swifter"
            Name of the desired format (Swift/Swifter/Swiftest)
        dask : bool, default False
            Use Dask to lazily load data (useful for very large datasets)
        **kwargs : Any
            A dictionary of additional keyword arguments.

        Returns
        -------
        xarray dataset
            Dataset containing the variables retrieved from the follow algorithm
        """
        verbose = kwargs.pop("verbose", self.verbose)

        if self.data is None:
            self.read_output_file(dask=dask, verbose=verbose)
        if codestyle == "Swift":
            try:
                with open("follow.in") as f:
                    line = f.readline()  # Parameter file (ignored because read_output_file already takes care of it
                    line = f.readline()  # PL file (ignored)
                    line = f.readline()  # TP file (ignored)
                    line = f.readline()  # ifol
                    i_list = [i for i in line.split(" ") if i.strip()]
                    ifol = int(i_list[0])
                    line = f.readline()  # nskp
                    i_list = [i for i in line.split(" ") if i.strip()]
                    nskp = int(i_list[0])
            except OSError:
                if verbose:
                    warnings.warn("No follow.in file found", stacklevel=2)
                ifol = None
                nskp = None
            fol = tool.follow_swift(self.data, ifol=ifol, nskp=nskp)
        else:
            fol = None

        if verbose:
            print("follow.out written")
        return fol

    def _scrub_init_cond(self):
        """
        Drops the oblateness terms from all bodies except the central body.
        """
        ds = self.init_cond

        # Drop any variables that may have been copied over from old runs. This is a list of all potential init_cond.nc variables
        ic_vars = [
            "rh",
            "vh",
            "gr_pseudo_vh",
            "a",
            "e",
            "inc",
            "capom",
            "omega",
            "capm",
            "varpi",
            "lam",
            "f",
            "cape",
            "capf",
            "Gmass",
            "mass",
            "radius",
            "rhill",
            "j2rp2",
            "j4rp4",
            "rot",
            "rotphase",
            "Ip",
            "id",
            "particle_type",
            "status",
            "c_lm",
            "ntp",
            "npl",
            "nplm",
        ]

        vars = [k for k in ic_vars if k in ds]
        ds = ds[vars]
        iscb = ds["particle_type"].compute() == constants.CB_TYPE_NAME
        cbname = ds["name"].where(iscb, drop=True).values[0]
        if "j2rp2" in ds:
            if "name" in ds.j2rp2.dims:
                j2rp2 = ds.j2rp2.sel(name=cbname)
            elif "id" in ds.j2rp2.dims:
                j2rp2 = ds.j2rp2.sel(id=0)
            else:
                j2rp2 = ds.j2rp2
            if np.isnan(j2rp2.values[0]):
                ds = ds.drop_vars("j2rp2")
            else:
                ds["j2rp2"] = j2rp2
        if "j4rp4" in ds:
            if "name" in ds.j4rp4.dims:
                j4rp4 = ds.j4rp4.sel(name=cbname)
            elif "id" in ds.j4rp4.dims:
                j4rp4 = ds.j4rp4.sel(id=0)
            else:
                j4rp4 = ds.j4rp4
            if np.isnan(j4rp4.values[0]):
                ds = ds.drop_vars("j4rp4")
            else:
                ds["j4rp4"] = j4rp4
        if "c_lm" in ds:
            if "name" in ds.c_lm.dims:
                c_lm = ds.c_lm.sel(name=cbname)
            elif "id" in ds.c_lm.dims:
                c_lm = ds.c_lm.sel(id=0)
            else:
                c_lm = ds.c_lm
            if np.any(np.isnan(c_lm.values)):
                ds = ds.drop_vars(["c_lm", "sign", "l", "m"], errors="ignore")
            else:
                ds["c_lm"] = c_lm
        if "rotphase" in ds:
            if "name" in ds.rotphase.dims:
                rotphase = ds.rotphase.sel(name=cbname)
            elif "id" in ds.rotphase.dims:
                rotphase = ds.rotphase.sel(id=0)
            else:
                rotphase = ds.rotphase
            if np.isnan(rotphase.values[0]):
                ds = ds.drop_vars("rotphase")
            else:
                ds["rotphase"] = rotphase

        # Drop any variables that may have been copied over from old runs.
        if self.param["IN_FORM"] == "EL":
            ds = ds.drop_vars(["rh", "vh"], errors="ignore")
        elif self.param["IN_FORM"] == "XV":
            ds = ds.drop_vars(
                [
                    "a",
                    "e",
                    "inc",
                    "capom",
                    "omega",
                    "capm",
                    "cape",
                    "capf",
                    "varpi",
                    "lam",
                    "f",
                ],
                errors="ignore",
            )
        self.init_cond = ds
        return

    def save(
        self,
        codename: Literal["Swiftest", "Swifter", "Swift"] | None = None,
        param_file: str | os.PathLike | None = None,
        param: dict | None = None,
        framenum: int = -1,
        **kwargs: Any,
    ) -> None:
        """
        Saves an xarray dataset to a set of input files.

        Parameters
        ----------
        codename : {"Swiftest", "Swifter", "Swift"}, optional
            Alternative name of the n-body code that the parameter file will be formatted for. Defaults to current instance
            variable self.codename
        param_file : str or path-like, optional
            Alternative file name of the input parameter file. Defaults to current instance variable self.param_file
        param : Dict, optional
            An alternative parameter dictionary to write out. Defaults to the current instance variable self.param
        framenum : int Default=-1
            Time frame to use to generate the initial conditions. If this argument is not passed, the default is to use the last frame in the dataset.
        **kwargs : Any
            A dictionary of additional keyword argument.

        Returns
        -------
        None
        """
        verbose = kwargs.pop("verbose", self.verbose)

        if codename is None:
            codename = self.codename
        if param_file is None:
            param_file = self.param_file
        if param is None:
            param = self.param

        if not self.simdir.exists():
            self.simdir.mkdir(parents=True, exist_ok=True)

        if not self.restart:
            if len(self.data) == 0:
                self.data = self.init_cond.copy(deep=True)
            else:
                self.init_cond = self.data.isel(time=[framenum]).copy(deep=True)
            self._scrub_init_cond()

        if codename == "Swiftest":
            infile_name = Path(self.simdir) / param["NC_IN"]
            io.swiftest_xr2infile(
                ds=self.init_cond,
                param=param,
                in_type=self.param["IN_TYPE"],
                infile_name=infile_name,
                verbose=verbose,
            )
            self.write_param(param_file=param_file, verbose=verbose, **kwargs)
        elif codename == "Swifter":
            swifter_param = io.swiftest2swifter_param(param)
            if "rhill" in self.data:
                swifter_param["RHILL_PRESENT"] = "YES"
            io.swifter_xr2infile(self.init_cond, swifter_param, self.simdir)
            self.write_param(
                codename=codename,
                param_file=param_file,
                param=swifter_param,
                verbose=verbose,
                **kwargs,
            )
        elif verbose:
            warnings.warn(f"Saving to {codename} not supported", stacklevel=2)

        return

    def clean(self, deep: bool = False, **kwargs) -> None:
        """
        Cleans up simulation directory by deleting old files (dump, logs, and output files).

        Parameters
        ----------
        deep : bool, default False
            If True, also delete the entire simulation directory.
        **kwargs : Any
            A dictionary of additional keyword arguments.

        Returns
        -------
        None
        """
        verbose = kwargs.pop("verbose", self.verbose)
        self.encounters = SwiftestDataset()
        self.collisions = SwiftestDataset()
        if deep:
            if verbose and self.simdir.exists():
                print(f"Removing simulation directory {self.simdir}")
            shutil.rmtree(self.simdir, ignore_errors=True)
            self.init_cond = SwiftestDataset()
            self.data = SwiftestDataset()
            return

        if verbose:
            print(f"Cleaning up simulation directory {self.simdir} of old output files.")

        old_files = [
            self.simdir / self.param["BIN_OUT"],
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

        # Clean out data structure and reset it to initial conditions
        if "time" in self.data:
            self.data = self.data.isel(time=[0])
        else:
            self.data = self.init_cond.copy(deep=True)
        return

    def _set_central_body(self, align_to_central_body_rotation: bool = False, **kwargs: Any):
        """
        Sets the central body to be the most massive body in the dataset. Cartesian position and velocity Cartesian coordinates are rotated  If align_to_central_body_rotation is True, the rotation pole is set to the z-axis.

        Parameters
        ----------
        align_to_central_body_rotation : bool, default False
            If True, the rotation pole is set to the z-axis.

        Returns
        -------
        None

        """
        verbose = kwargs.pop("verbose", self.verbose)

        if "Gmass" not in self.data:
            if verbose:
                warnings.warn(
                    "No bodies with Gmass values found in dataset. Cannot set central body.",
                    stacklevel=2,
                )
            return

        if "name" in self.data.dims:
            cbidx = self.data.Gmass.argmax(dim="name").values[0]
            cbid = self.data.id.isel(name=cbidx).values[()]
            cbname = self.data.name.isel(name=cbidx).values[()]
        elif "id" in self.data.dims:
            cbidx = self.data.Gmass.argmax(dim="id").values[0]
            cbid = self.data.id.isel(id=cbidx).values[()]
            cbname = self.data.name.isel(id=cbidx).values[()]
        else:
            raise ValueError("No 'name' or 'id' dimensions found in dataset.")

        recompute_el = False
        if cbid != 0:
            recompute_el = True
            if "name" in self.data.dims:
                if 0 in self.data.id.values:
                    name_0 = self.data.name.where(self.data.id == 0, drop=True).values[()]
                    self.data["id"].loc[dict(name=name_0)] = cbid
                self.data["id"].loc[dict(name=cbname)] = 0
                self.data["particle_type"].loc[dict(name=cbname)] = constants.CB_TYPE_NAME
            else:
                if 0 in self.data.id.values:
                    self.data["id"].loc[dict(id=0)] = cbid
                self.data["id"].loc[dict(id=cbid)] = 0
                self.data["particle_type"].loc[dict(id=cbid)] = constants.CB_TYPE_NAME

        # Ensure that the central body is at the origin
        if "name" in self.data.dims:
            cbda = self.data.sel(name=cbname)
        else:
            cbda = self.data.sel(id=cbid)

        pos_skip = ["space", "Ip", "rot"]
        for var in self.data.variables:
            if "space" in self.data[var].dims and var not in pos_skip and not np.isnan(self.data[var].values).any():
                if np.any(cbda[var].values != 0.0):
                    recompute_el = True
                    self.data[var] = self.data[var] - cbda[var]

        # If the central body origin has changed and we expect the system to be aligned with its rotation frame, then rotate the system
        # before computing the orbital elements
        if align_to_central_body_rotation and "rot" in cbda and recompute_el:
            if "rot" in cbda and not np.isnan(cbda.rot.isel(time=0).values).any():
                self.data = self.data.rotate(pole=cbda.rot.isel(time=0).values[()])

        if recompute_el:
            self.data = self.data.xv2el()

        self.set_distance_range(rmin=cbda.radius.values.item(), verbose=verbose)

        return

    @property
    def param(self) -> dict:
        """
        Dict: A dictionary of simulation parameters. These are stored in the param.in file.
        """
        return self._param

    @param.setter
    def param(self, value: dict) -> None:
        if not isinstance(value, dict):
            raise TypeError("Parameter value must be a dictionary")
        self._param = value
        return

    @property
    def data(self) -> SwiftestDataset:
        """
        SwiftestDataset: A dataset containing the simulation data.
        """
        return self._data

    @data.setter
    def data(self, value: SwiftestDataset) -> None:
        if not isinstance(value, SwiftestDataset):
            if isinstance(value, xr.Dataset):
                value = SwiftestDataset(value)
            else:
                raise TypeError("Data value must be a SwiftestDataset")
        self._data = value
        return

    @property
    def init_cond(self) -> SwiftestDataset:
        """
        SwiftestDataset: A dataset containing the initial conditions.
        """
        return self._init_cond

    @init_cond.setter
    def init_cond(self, value: SwiftestDataset) -> None:
        if not isinstance(value, SwiftestDataset):
            if isinstance(value, xr.Dataset):
                value = SwiftestDataset(value)
            else:
                raise TypeError("Initial conditions value must be a SwiftestDataset")
        self._init_cond = value
        return

    @property
    def encounters(self) -> SwiftestDataset:
        """
        SwiftestDataset: A dataset containing the encounter history
        """
        return self._encounters

    @encounters.setter
    def encounters(self, value: SwiftestDataset) -> None:
        if not isinstance(value, SwiftestDataset):
            if isinstance(value, xr.Dataset):
                value = SwiftestDataset(value)
            else:
                raise TypeError("Encounters value must be a SwiftestDataset")
        self._encounters = value
        return

    @property
    def collisions(self) -> SwiftestDataset:
        """
        SwiftestDataset: A dataset containing the collision history.
        """
        return self._collisions

    @collisions.setter
    def collisions(self, value: SwiftestDataset) -> None:
        if not isinstance(value, SwiftestDataset):
            if isinstance(value, xr.Dataset):
                value = SwiftestDataset(value)
            else:
                raise TypeError("Collisions value must be a SwiftestDataset")
        self._collisions = value
        return

    @property
    def simdir(self) -> Path:
        """
        Path: The simulation directory.
        """
        return self._simdir

    @simdir.setter
    def simdir(self, value: os.PathLike | Path) -> None:
        try:
            value = Path(value).resolve()
        except:
            raise TypeError("Simulation directory value must be a a valid path")

        if value.exists():
            if not value.is_dir():
                msg = f"Cannot create the {self.simdir} directory: File exists."
                msg += "\nDelete the file or change the location of param_file"
                raise NotADirectoryError(msg)
        self._simdir = value
        return

    @property
    def param_file(self) -> Path:
        """
        Path: The parameter file.
        """
        return self._param_file

    @param_file.setter
    def param_file(self, value: os.PathLike | Path) -> None:
        if not os.path.exists(self.simdir / value):
            self.write_param(param_file=value)
        param_path = self.simdir / value
        if not param_path.exists():
            raise FileNotFoundError(f"Parameter file {param_path} not found.")
        self.read_param(param_file=param_path)
        self._param_file = value
        return

    @property
    def MU_name(self) -> str:
        """
        str: The name of the mass unit currently defined in the simulation.
        """
        return self._MU_name

    @MU_name.setter
    def MU_name(self, value: str) -> None:
        if value is None:
            value = "MU"
        if not isinstance(value, str):
            raise TypeError("Mass unit name value must be a string")
        self._MU_name = value
        return

    @property
    def DU_name(self) -> str:
        """
        str: The name of the distance unit currently defined in the simulation.
        """
        return self._DU_name

    @DU_name.setter
    def DU_name(self, value: str) -> None:
        if value is None:
            value = "DU"
        if not isinstance(value, str):
            raise TypeError("Distance unit name value must be a string")
        self._DU_name = value
        return

    @property
    def TU_name(self) -> str:
        """
        str: The name of the time unit currently defined in the simulation.
        """
        return self._TU_name

    @TU_name.setter
    def TU_name(self, value: str) -> None:
        if value is None:
            value = "TU"
        if not isinstance(value, str):
            raise TypeError("Time unit name value must be a string")
        self._TU_name = value
        return

    @property
    def MU2KG(self) -> FloatLike:
        """
        np.float64: The conversion factor from mass unit to kilograms.
        """
        return self._MU2KG

    @MU2KG.setter
    def MU2KG(self, value: FloatLike) -> None:
        if not isinstance(value, (float, int, np.number)):
            raise TypeError("Mass unit to kilogram conversion value must be a float")
        self._MU2KG = np.longdouble(value)
        self._KG2MU = 1.0 / self._MU2KG
        return

    @property
    def KG2MU(self) -> FloatLike:
        """
        np.float64: The conversion factor from kilograms to mass unit.
        """
        return self._KG2MU

    @KG2MU.setter
    def KG2MU(self, value: FloatLike) -> None:
        if not isinstance(value, (float, int, np.number)):
            raise TypeError("Kilogram to mass unit conversion value must be a float")
        self._KG2MU = np.longdouble(value)
        self._MU2KG = 1.0 / self._KG2MU
        return

    @property
    def DU2M(self) -> FloatLike:
        """
        np.float64: The conversion factor from distance unit to meters.
        """
        return self._DU2M

    @DU2M.setter
    def DU2M(self, value: FloatLike) -> None:
        if not isinstance(value, (float, int, np.number)):
            raise TypeError("Distance unit to meter conversion value must be a float")
        self._DU2M = np.longdouble(value)
        self._M2DU = 1.0 / self._DU2M
        return

    @property
    def M2DU(self) -> FloatLike:
        """
        np.float64: The conversion factor from meters to distance unit.
        """
        return self._M2DU

    @M2DU.setter
    def M2DU(self, value: FloatLike) -> None:
        if not isinstance(value, (float, int, np.number)):
            raise TypeError("Meter to distance unit conversion value must be a float")
        self._M2DU = np.longdouble(value)
        self._DU2M = 1.0 / self._M2DU
        return

    @property
    def TU2S(self) -> FloatLike:
        """
        np.float64: The conversion factor from time unit to seconds.
        """
        return self._TU2S

    @TU2S.setter
    def TU2S(self, value: FloatLike) -> None:
        if not isinstance(value, (float, int, np.number)):
            raise TypeError("Time unit to second conversion value must be a float")
        self._TU2S = np.longdouble(value)
        self._S2TU = 1.0 / self._TU2S
        return

    @property
    def S2TU(self) -> FloatLike:
        """
        np.float64: The conversion factor from seconds to time unit.
        """
        return self._S2TU

    @S2TU.setter
    def S2TU(self, value: FloatLike) -> None:
        if not isinstance(value, (float, int, np.number)):
            raise TypeError("Second to time unit conversion value must be a float")
        self._S2TU = np.longdouble(value)
        self._TU2S = 1.0 / self._S2TU
        return

    @property
    def GU(self) -> FloatLike:
        """
        np.float64: The gravitational constant in the simulation units.
        """
        return self._GU

    @GU.setter
    def GU(self, value: FloatLike) -> None:
        if not isinstance(value, (float, int, np.number)):
            raise TypeError("Gravitational constant value must be a float")
        self._GU = np.float64(value)
        return

    @property
    def integrator(self) -> str:
        """
        Literal["symba","rmvs","whm","helio"]: The name of the integrator used in the simulation. Currently supports "symba", "rmvs", "whm", and "helio".
        """
        return self._integrator

    @integrator.setter
    def integrator(self, value: str) -> None:
        if not isinstance(value, str):
            raise TypeError("Integrator value must be a string")

        valid_integrator = ["symba", "rmvs", "whm", "helio"]
        if value.lower() not in valid_integrator:
            raise ValueError(
                f"{value} is not a valid integrator. Valid options are ",
                ",".join(valid_integrator),
            )

        self._integrator = value.lower()
        return

    @property
    def codename(self) -> str:
        """
        Literal["Swiftest","Swifter","Swift"]: The name of the n-body code used in the simulation. Currently supports "Swiftest", "Swifter", and "Swift".
        """
        return self._codename

    @codename.setter
    def codename(self, value: str) -> None:
        if not isinstance(value, str):
            raise TypeError("Code name value must be a string")

        valid_codename = ["Swiftest", "Swifter", "Swift"]
        if value.title() not in valid_codename:
            raise ValueError(
                f"{value} is not a valid code name. Valid options are ",
                ",".join(valid_codename),
            )

        self._codename = value.title()
        return

    @property
    def verbose(self) -> bool:
        """
        bool: A boolean value indicating whether verbose output is enabled.
        """
        return self._verbose

    @verbose.setter
    def verbose(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("verbose value must be a boolean")
        self._verbose = value
        return

    @property
    def restart(self) -> bool:
        """
        bool: A boolean value indicating whether the run should be restarted from a previous output.
        """
        return self._restart

    @restart.setter
    def restart(self, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError("restart value must be a boolean")
        self._restart = value
        return

    @property
    def ephemeris_date(self) -> str:
        """
        ISO-formatted date sto use when obtaining the ephemerides in the format YYYY-MM-DD.
        """
        return self._ephemeris_date

    @ephemeris_date.setter
    def ephemeris_date(self, value: str) -> None:
        """
        Sets the date to use when obtaining the ephemerides.

        Parameters
        ----------
        ephemeris_date : str, optional
            Date to use when obtaining the ephemerides.
            Valid options are "today", "MBCL", or date in the format YYYY-MM-DD.
        **kwargs : Any
            A dictionary of additional keyword argument. This allows this method to be called by the more general
            set_parameter method, which takes all possible Simulation parameters as arguments, so these are ignored.
        if not isinstance(str, bool):
            raise TypeError("ephemeris_date value must be a string")
        """
        if value.upper() == "MBCL":
            value = constants.MINTON_BCL
        elif value.upper() == "TODAY":
            value = datetime.date.today().isoformat()
        else:
            try:
                datetime.datetime.fromisoformat(value)
            except:
                valid_date_args = ['"MBCL"', '"TODAY"', '"YYYY-MM-DD"']
                msg = (
                    f"{value} is not a valid format. Valid options include:",
                    ", ".join(valid_date_args),
                )
                msg += "\nUsing MBCL for date."
                warnings.warn(msg, stacklevel=2)
                value = constants.MINTON_BCL

        self._ephemeris_date = value

        return
