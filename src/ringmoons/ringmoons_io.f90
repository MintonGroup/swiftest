! Copyright 2026 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(ringmoons) s_ringmoons_io
    use swiftest
contains        
    module subroutine ringmoons_io_netcdf_dump(self, param)
        implicit none
        class(ringmoons_storage), intent(inout)        :: self   
            !! ringmoons storage object
        class(base_parameters),   intent(inout)        :: param  
            !! Current run configuration parameters 

        return
    end subroutine ringmoons_io_netcdf_dump

    module subroutine ringmoons_io_netcdf_initialize_output(self, param)
        implicit none
        class(ringmoons_netcdf_parameters), intent(inout) :: self    
            !! Parameters used to identify a particular NetCDF dataset
        class(base_parameters),             intent(in)    :: param   

        return
    end subroutine ringmoons_io_netcdf_initialize_output

    module subroutine ringmoons_io_netcdf_open(self, param, readonly)
        implicit none
        class(ringmoons_netcdf_parameters), intent(inout) :: self     
            !! Parameters used to identify a particular NetCDF dataset
        class(base_parameters),             intent(in)    :: param    
            !! Current run configuration parameters
        logical, optional,                  intent(in)    :: readonly 
            !! Logical flag indicating that this should be open read only

        return
    end subroutine ringmoons_io_netcdf_open

end submodule s_ringmoons_io
