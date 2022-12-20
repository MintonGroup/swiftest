!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (encounter) s_encounter_setup
   use swiftest
contains

   module subroutine encounter_setup_aabb(self, n, n_last)
      !! author: David A. Minton
      !!
      !! Sets up or modifies an axis-aligned bounding box structure.
      implicit none
      ! Arguments
      class(encounter_bounding_box), intent(inout) :: self   !! Swiftest encounter structure
      integer(I4B),                  intent(in)    :: n      !! Number of objects with bounding box extents
      integer(I4B),                  intent(in)    :: n_last !! Number of objects with bounding box extents the previous time this was called
      ! Internals
      integer(I4B) :: next, next_last, k, dim
      integer(I4B), dimension(:), allocatable :: itmp

      next = 2 * n
      next_last = 2 * n_last

      if (n > n_last) then ! The number of bodies has grown. Resize and append the new bodies
         do dim = 1, SWEEPDIM
            allocate(itmp(next))
            if (n_last > 0) itmp(1:next_last) = self%aabb(dim)%ind(1:next_last)
            call move_alloc(itmp, self%aabb(dim)%ind)
            self%aabb(dim)%ind(next_last+1:next) = [(k, k = next_last+1, next)]
         end do
      else ! The number of bodies has gone down. Resize and chop of the old indices
         do dim = 1, SWEEPDIM
            allocate(itmp(next))
            itmp(1:next) = pack(self%aabb(dim)%ind(1:next_last), self%aabb(dim)%ind(1:next_last) <= next)
            call move_alloc(itmp, self%aabb(dim)%ind)
         end do
      end if

      do dim = 1, SWEEPDIM
         if (allocated(self%aabb(dim)%ibeg)) deallocate(self%aabb(dim)%ibeg)
         allocate(self%aabb(dim)%ibeg(n))
         if (allocated(self%aabb(dim)%iend)) deallocate(self%aabb(dim)%iend)
         allocate(self%aabb(dim)%iend(n))
      end do

      return
   end subroutine encounter_setup_aabb


   module subroutine encounter_setup_list(self, n)
      !! author: David A. Minton
      !!
      !! A constructor that sets the number of encounters and allocates and initializes all arrays  
      !!
      implicit none
      ! Arguments
      class(encounter_list), intent(inout) :: self !! Swiftest encounter structure
      integer(I8B),          intent(in)    :: n    !! Number of encounters to allocate space for
      ! Internals
      integer(I8B) :: i

      if (n < 0) return
      call self%dealloc()

      self%nenc = n
      if (n == 0_I8B) return
      self%t = 0.0_DP

      allocate(self%tcollision(n))
      allocate(self%lvdotr(n))
      allocate(self%lclosest(n))
      allocate(self%status(n))
      allocate(self%index1(n))
      allocate(self%index2(n))
      allocate(self%id1(n))
      allocate(self%id2(n))
      allocate(self%r1(NDIM,n))
      allocate(self%r2(NDIM,n))
      allocate(self%v1(NDIM,n))
      allocate(self%v2(NDIM,n))

      self%tcollision(:) = 0.0_DP
      self%lvdotr(:) = .false.
      self%lclosest(:) = .false.
      self%status(:) = INACTIVE
      self%index1(:) = 0
      self%index2(:) = 0
      self%id1(:) = 0
      self%id2(:) = 0
      self%r1(:,:) = 0.0_DP
      self%r2(:,:) = 0.0_DP
      self%v1(:,:) = 0.0_DP
      self%v2(:,:) = 0.0_DP

      return
   end subroutine encounter_setup_list

end submodule s_encounter_setup

 