!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_util_final
   use swiftest
contains

   module subroutine util_final_storage(self)
      !! author: David A. Minton
      !!
      !! Finalizer for the storage data type
      implicit none
      ! Arguments
      type(swiftest_storage(*)) :: self
      ! Internals
      integer(I4B) :: i

      do i = 1, self%nframes
         if (allocated(self%frame(i)%item)) deallocate(self%frame(i)%item)
      end do

      return
   end subroutine util_final_storage

   module subroutine util_final_storage_frame(self)
      !! author: David A. Minton
      !!
      !! Finalizer for the storage frame data type
      implicit none
      type(swiftest_storage_frame) :: self

      if (allocated(self%item)) deallocate(self%item)

      return
   end subroutine util_final_storage_frame


   module subroutine util_final_system(self)
      !! author: David A. Minton
      !!
      !! Finalize the swiftest nbody system object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_nbody_system),  intent(inout) :: self !! Swiftest nbody system object

      if (allocated(self%cb)) deallocate(self%cb)
      if (allocated(self%pl)) deallocate(self%pl)
      if (allocated(self%tp)) deallocate(self%tp)
      if (allocated(self%tp_discards)) deallocate(self%tp_discards)
      if (allocated(self%pl_discards)) deallocate(self%pl_discards)

      return
   end subroutine util_final_system


end submodule s_util_final
