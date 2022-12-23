!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (collision) s_collision_setup
   use swiftest
contains

   module subroutine collision_setup_system(self, nbody_system)
      !! author: David A. Minton
      !!
      !! Initializer for the encounter collision system. Sets up impactors and the before/after snapshots,
      !! but not fragments. Those are setup later when the number of fragments is known.
      implicit none
      ! Arguments
      class(collision_merge),  intent(inout) :: self         !! Encounter collision system object
      class(base_nbody_system), intent(in)    :: nbody_system !! Current nbody system. Used as a mold for the before/after snapshots

      call self%setup_impactors()
      if (allocated(self%before)) deallocate(self%before)
      if (allocated(self%after)) deallocate(self%after)

      allocate(self%before, mold=nbody_system)
      allocate(self%after,  mold=nbody_system)

      return
   end subroutine collision_setup_system


   module subroutine collision_setup_impactors_system(self)
      !! author: David A. Minton
      !!
      !! Initializer for the impactors for the encounter collision system. Deallocates old impactors before creating new ones
      implicit none
      ! Arguments
      class(collision_merge), intent(inout) :: self   !! Encounter collision system object

      if (allocated(self%impactors)) deallocate(self%impactors)
      allocate(collision_impactors :: self%impactors)

      return
   end subroutine collision_setup_impactors_system


   module subroutine collision_setup_fragments_system(self, nfrag)
      !! author: David A. Minton
      !!
      !! Initializer for the fragments of the collision system. 
      implicit none
      ! Arguments
      class(collision_merge), intent(inout) :: self  !! Encounter collision system object
      integer(I4B),            intent(in)    :: nfrag !! Number of fragments to create

      if (allocated(self%fragments)) deallocate(self%fragments)
      allocate(collision_fragments(nfrag) :: self%fragments)
      self%fragments%nbody = nfrag

      return
   end subroutine collision_setup_fragments_system

end submodule s_collision_setup

 