!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(swiftest_classes) s_util_snapshot
   use swiftest
contains

   module subroutine util_snapshot_system(self, param, system, t, arg)
      !! author: David A. Minton
      !!
      !! Takes a snapshot of the system for later file storage
      implicit none
      ! Arguments
      class(swiftest_storage(*)),   intent(inout) :: self   !! Swiftest storage object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object to store
      real(DP),                     intent(in), optional :: t      !! Time of snapshot if different from system time
      character(*),                 intent(in), optional :: arg    !! Optional argument (needed for extended storage type used in collision snapshots)

      self%iframe = self%iframe + 1
      self%nt = self%iframe
      self%frame(self%iframe) = system ! Store a snapshot of the system for posterity
      self%nid = self%nid + 1 ! Central body
      if (allocated(system%pl)) self%nid = self%nid + system%pl%nbody
      if (allocated(system%tp)) self%nid = self%nid + system%tp%nbody

      return
   end subroutine util_snapshot_system

end submodule s_util_snapshot