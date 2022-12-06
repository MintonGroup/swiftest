!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_util_reset
   use swiftest
contains

   module subroutine util_reset_storage(self)
      !! author: David A. Minton
      !!
      !! Resets a storage object by deallocating all items and resetting the frame counter to 0
      implicit none
      ! Arguments
      class(swiftest_storage(*)), intent(inout) :: self !! Swiftest storage object
      ! Internals
      integer(I4B) :: i

      do i = 1, self%nframes
         if (allocated(self%frame(i)%item)) deallocate(self%frame(i)%item)
      end do
      self%iframe = 0

      return
   end subroutine util_reset_storage

end submodule s_util_reset