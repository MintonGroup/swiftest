from swiftest import swiftestio


class Simulation:
    """
    This is a class that define the basic Swift/Swifter/Swiftest simulation object
    """
    def __init__(self, codename="Swiftest", source=""):
        if source == "":
            self.param = {
                '! VERSION': f"Swiftest parameter input",
                'T0': "0.0",
                'TSTOP': "0.0",
                'DT': "0.0",
                'PL_IN': "",
                'TP_IN': "",
                'CB_IN': "",
                'IN_TYPE': "REAL8",
                'ISTEP_OUT': "1",
                'ISTEP_DUMP': "1",
                'BIN_OUT': "bin.dat",
                'OUT_TYPE': 'REAL8',
                'OUT_FORM': "XV",
                'OUT_STAT': "REPLACE",
                'J2': "0.0",
                'J4': "0.0",
                'CHK_RMIN': "-1.0",
                'CHK_RMAX': "-1.0",
                'CHK_EJECT': "-1.0",
                'CHK_QMIN': "-1.0",
                'CHK_QMIN_COORD': "HELIO",
                'CHK_QMIN_RANGE': "",
                'ENC_OUT': "",
                'MTINY': "-1.0",
                'MU2KG': "-1.0",
                'TU2S': "-1.0",
                'DU2M': "-1.0",
                'GU': "-1.0",
                'EXTRA_FORCE': "NO",
                'BIG_DISCARD': "NO",
                'CHK_CLOSE': "NO",
                'FRAGMENTATION': "NO",
                'MTINY_SET': "NO",
                'ROTATION': "NO",
                'TIDES': "NO",
                'ENERGY': "NO",
                'GR': "NO",
                'YARKOVSKY': "NO",
                'YORP': "NO",
            }
        else:
            self.read_param(source, codename)
        return
    
    def read_param(self, param_file_name, codename="Swiftest"):
        if codename == "Swiftest":
            self.param = swiftestio.read_swiftest_param(param_file_name)
            self.codename = "Swiftest"
        elif codename == "Swifter":
            self.param = swiftestio.read_swifter_param(param_file_name)
            self.codename = "Swifter"
        elif codename == "Swift":
            self.param = swiftestio.read_swift_param(param_file_name)
            self.codename = "Swift"
        else:
            print(f'{codename} is not a recognized code name. Valid options are "Swiftest", "Swifter", or "Swift".')
            self.codename = "Unknown"
        return
    
    def write_param(self, param_file_name):
        # Check to see if the parameter type matches the output type. If not, we need to convert
        codename = self.param['! VERSION'].split()[0]
        if codename == "Swifter" or codename == "Swiftest":
            swiftestio.write_labeled_param(self.param, param_file_name)
        elif codename == "Swift":
            swiftestio.write_swift_param(self.param, param_file_name)
        else:
            print('Cannot process unknown code type. Call the read_param method with a valid code name. Valid options are "Swiftest", "Swifter", or "Swift".')
        return
    
    def convert(self, param_file_name, newcodename="Swiftest", plname="pl.swiftest.in", tpname="tp.swiftest.in", cbname="cb.swiftest.in"):
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
                self.param = swiftestio.swifter2swiftest(self.param, plname, tpname, cbname)
            else:
                goodconversion = False
        elif self.codename == "Swift":
            if newcodename == "Swifter":
                self.param = swiftestio.swift2swifter(self.param, plname, tpname)
            elif newcodename == "Swiftest":
                self.param = swiftestio.swift2swiftest(self.param, plname, tpname, cbname)
            else:
                goodconversion = False

        if goodconversion:
            self.write_param(param_file_name)
        else:
            print(f"Conversion from {self.codename} to {newcodename} is not supported.")
        return oldparam
    
    def bin2xr(self):
        if self.codename == "Swiftest":
            self.ds = swiftestio.swiftest2xr(self.param)
        elif self.codename == "Swifter":
            self.ds = swiftestio.swifter2xr(self.param)
        elif self.codename == "Swift":
            print("Reading Swift simulation data is not implemented yet")
        else:
            print('Cannot process unknown code type. Call the read_param method with a valid code name. Valid options are "Swiftest", "Swifter", or "Swift".')
        return
    

