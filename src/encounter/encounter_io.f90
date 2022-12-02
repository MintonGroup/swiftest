!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (encounter_classes) s_encounter_io
   use swiftest
contains

   module subroutine encounter_io_dump_storage_list(self, param)
      !! author: David A. Minton
      !!
      !! Dumps the time history of an encounter to file.
      implicit none
      ! Arguments
      class(encounter_storage(*)), intent(inout) :: self   !! Encounter storage object
      class(swiftest_parameters),  intent(inout) :: param  !! Current run configuration parameters 
   end subroutine encounter_io_dump_storage_list

   module subroutine encounter_io_initialize_output(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a NetCDF encounter file system and defines all variables.
      implicit none
      ! Arguments
      class(encounter_io_parameters), intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),     intent(in)    :: param   !! Current run configuration parameters

      return
   end subroutine encounter_io_initialize_output

   module subroutine encounter_io_open_file(self, param, readonly)
      !! author: David A. Minton
      !!
      !! Opens a NetCDF encounter file and does the variable inquiries to activate variable ids
      implicit none
      ! Arguments
      class(encounter_io_parameters), intent(inout) :: self     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),     intent(in)    :: param    !! Current run configuration parameters
      logical, optional,              intent(in)    :: readonly !! Logical flag indicating that this should be open read only

      return
   end subroutine encounter_io_open_file

end submodule s_encounter_io