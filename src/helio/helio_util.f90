!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest.
!! If not, see: https://www.gnu.org/licenses.

submodule(helio) s_helio_util
   use swiftest
contains

      module subroutine helio_util_setup_initialize_system(self, system_history, param)
      !! author: David A. Minton
      !!
      !! Initialize a Helio nbody system from files, converting all heliocentric quantities to barycentric.
      !!
      implicit none
      ! Arguments
      class(helio_nbody_system),               intent(inout) :: self   !! Helio nbody system object
      class(swiftest_storage),    allocatable, intent(inout) :: system_history !! Stores the system history between output dumps
      class(swiftest_parameters),              intent(inout) :: param          !! Current run configuration parameters 

      call swiftest_util_setup_initialize_system(self, system_history, param)
      call self%pl%sort("mass", ascending=.false.)
      call self%pl%vh2vb(self%cb)
      call self%tp%h2b(self%cb)

      ! Make sure that the discard list gets allocated initially
      call self%tp_discards%setup(0, param)
      call self%pl%set_mu(self%cb)
      call self%tp%set_mu(self%cb)

      if (param%lgr .and. param%in_type == "ASCII") then !! pseudovelocity conversion for NetCDF input files is handled by NetCDF routines
         call self%pl%v2pv(param)
         call self%tp%v2pv(param)
      end if

      call self%pl%flatten(param)

      return
   end subroutine helio_util_setup_initialize_system

end submodule s_helio_util
