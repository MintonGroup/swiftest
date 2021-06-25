from swiftest import swiftestio
class Simulation:
    """
    This is a class that define the basic Swift/Swifter/Swiftest simulation object
    """
    def __init__(self, param_file_name, codename="Swiftest"):
        self.read_param(param_file_name, codename)
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
        if self.codename == "Swifter":
            if newcodename == "Swiftest":
                self.param = swiftestio.swifter2swiftest(self.param, plname, tpname, cbname)
        self.write_param(param_file_name)
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
    

