!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(helio_classes) s_helio_util
   use swiftest
contains

   module subroutine helio_util_final_pl(self)
      !! author: David A. Minton
      !!
      !! Finalize the Helio massive body object - deallocates all allocatables
      implicit none
      ! Arguments
      type(helio_pl),  intent(inout) :: self !! Helio massive body object

      call self%dealloc()

      return
   end subroutine helio_util_final_pl


   module subroutine helio_util_final_system(self)
      !! author: David A. Minton
      !!
      !! Finalize the Helio nbody system object - deallocates all allocatables
      implicit none
      ! Arguments
      type(helio_nbody_system),  intent(inout) :: self !! Helio nbody system object

      call whm_util_final_system(self%whm_nbody_system)

      return
   end subroutine helio_util_final_system


   module subroutine helio_util_final_tp(self)
      !! author: David A. Minton
      !!
      !! Finalize the Helio test particle object - deallocates all allocatables
      implicit none
      ! Arguments
      type(helio_tp),  intent(inout) :: self !! Helio test particle object

      call self%dealloc()

      return
   end subroutine helio_util_final_tp

end submodule s_helio_util