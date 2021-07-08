from swiftest import io
from swiftest import init_cond
from swiftest import tool
from swiftest import constants
from datetime import date
import xarray as xr

class Simulation:
    """
    This is a class that define the basic Swift/Swifter/Swiftest simulation object
    """
    def __init__(self, codename="Swiftest", param_file=""):
        self.ds = xr.Dataset()
        self.param = {
            '! VERSION': f"Swiftest parameter input",
            'T0': "0.0",
            'TSTOP': "0.0",
            'DT': "0.0",
            'PL_IN': "pl.in",
            'TP_IN': "tp.in",
            'CB_IN': "cb.in",
            'IN_TYPE': "ASCII",
            'ISTEP_OUT': "1",
            'ISTEP_DUMP': "1",
            'BIN_OUT': "bin.dat",
            'OUT_TYPE': 'REAL8',
            'OUT_FORM': "EL",
            'OUT_STAT': "REPLACE",
            'CHK_RMAX': "1000.0",
            'CHK_EJECT': "1000.0",
            'CHK_RMIN': constants.RSun / constants.AU2M,
            'CHK_QMIN': constants.RSun / constants.AU2M,
            'CHK_QMIN_COORD': "HELIO",
            'CHK_QMIN_RANGE': f"{constants.RSun / constants.AU2M} 1000.0",
            'ENC_OUT': "enc.dat",
            'MU2KG': constants.MSun,
            'TU2S': constants.JD2S,
            'DU2M': constants.AU2M,
            'EXTRA_FORCE': "NO",
            'BIG_DISCARD': "NO",
            'CHK_CLOSE': "YES",
            'FRAGMENTATION': "NO",
            'ROTATION': "NO",
            'TIDES': "NO",
            'ENERGY': "NO",
            'GR': "NO",
            'YARKOVSKY': "NO",
            'YORP': "NO",
            'MTINY' : "0.0"
        }
        self.codename = codename
        if param_file != "" :
            self.read_param(param_file, codename)
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
        self.ds = init_cond.solar_system_horizons(plname, idval, self.param, date, self.ds)
        return
    
    def read_param(self, param_file, codename="Swiftest"):
        if codename == "Swiftest":
            self.param = io.read_swiftest_param(param_file, self.param)
            self.codename = "Swiftest"
        elif codename == "Swifter":
            self.param = io.read_swifter_param(param_file)
            self.codename = "Swifter"
        elif codename == "Swift":
            self.param = io.read_swift_param(param_file)
            self.codename = "Swift"
        else:
            print(f'{codename} is not a recognized code name. Valid options are "Swiftest", "Swifter", or "Swift".')
            self.codename = "Unknown"
        return
    
    def write_param(self, param_file, param=None):
        if param is None:
            param = self.param
        # Check to see if the parameter type matches the output type. If not, we need to convert
        codename = param['! VERSION'].split()[0]
        if codename == "Swifter" or codename == "Swiftest":
            io.write_labeled_param(param, param_file)
        elif codename == "Swift":
            io.write_swift_param(param, param_file)
        else:
            print('Cannot process unknown code type. Call the read_param method with a valid code name. Valid options are "Swiftest", "Swifter", or "Swift".')
        return
    
    def convert(self, param_file, newcodename="Swiftest", plname="pl.swiftest.in", tpname="tp.swiftest.in", cbname="cb.swiftest.in", conversion_questions={}):
        """
        Converts simulation input files from one code type to another (Swift, Swifter, or Swiftest). Returns the old parameter configuration.
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
        if self.codename == "Swiftest":
            self.ds = io.swiftest2xr(self.param)
            print('Swiftest simulation data stored as xarray DataSet .ds')
        elif self.codename == "Swifter":
            self.ds = io.swifter2xr(self.param)
            print('Swifter simulation data stored as xarray DataSet .ds')
        elif self.codename == "Swift":
            print("Reading Swift simulation data is not implemented yet")
        else:
            print('Cannot process unknown code type. Call the read_param method with a valid code name. Valid options are "Swiftest", "Swifter", or "Swift".')
        return
    
    def follow(self, codestyle="Swifter"):
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
        
        print('follow.out written')
        return fol
    
    def save(self, param_file, framenum=-1, codename="Swiftest"):
        if codename == "Swiftest":
            io.swiftest_xr2infile(self.ds, self.param, framenum)
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