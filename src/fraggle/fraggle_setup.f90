!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (fraggle) s_fraggle_setup
   use swiftest
   use symba
contains

   module subroutine fraggle_setup_fragments_system(self, nfrag)
      !! author: David A. Minton
      !!
      !! Initializer for the fragments of the collision system. 
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: self  !! Encounter collision system object
      integer(I4B),          intent(in)    :: nfrag !! Number of fragments to create

      if (allocated(self%fragments)) deallocate(self%fragments)
      allocate(fraggle_fragments(nbody=nfrag) :: self%fragments)
      self%fragments%nbody = nfrag

      return
   end subroutine fraggle_setup_fragments_system

end submodule s_fraggle_setup