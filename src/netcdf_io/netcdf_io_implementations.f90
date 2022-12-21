!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (netcdf_io) s_netcdf_io_implementations
   use netcdf
contains

   module subroutine netcdf_io_check(status, call_identifier)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Checks the status of all NetCDF operations to catch errors
      use netcdf
      implicit none
      ! Arguments
      integer, intent (in) :: status !! The status code returned by a NetCDF function
      character(len=*), intent(in), optional :: call_identifier !! String that indicates which calling function caused the error for diagnostic purposes

      if(status /= nf90_noerr) then
         if (present(call_identifier)) write(*,*) "NetCDF error in ",trim(call_identifier)
         write(*,*) trim(nf90_strerror(status))
         call util_exit(FAILURE)
      end if

      return
   end subroutine netcdf_io_check


   module subroutine netcdf_io_close(self)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Closes a NetCDF file
      use netcdf
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self   !! Parameters used to identify a particular NetCDF dataset

      if (self%lfile_is_open) then
         call netcdf_io_check( nf90_close(self%id), "netcdf_io_close" )
         self%lfile_is_open = .false.
      end if

      return
   end subroutine netcdf_io_close


   module subroutine netcdf_io_sync(self)
      !! author: David A. Minton
      !!
      !! Syncrhonize the disk and memory buffer of the NetCDF file (e.g. commit the frame files stored in memory to disk) 
      !!    
      use netcdf
      implicit none
      ! Arguments
      class(netcdf_parameters), intent(inout) :: self !! Parameters used to identify a particular NetCDF dataset

      call netcdf_io_check( nf90_sync(self%id), "netcdf_io_sync nf90_sync"  )

      return
   end subroutine netcdf_io_sync



end submodule s_netcdf_io_implementations
