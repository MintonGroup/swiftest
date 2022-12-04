!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (encounter_classes) s_encounter_util
   use swiftest
contains

   module subroutine encounter_util_append_list(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(encounter_list), intent(inout) :: self         !! Swiftest encounter list object
      class(encounter_list), intent(in)    :: source       !! Source object to append
      logical, dimension(:), intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nold, nsrc

      nold = self%nenc
      nsrc = source%nenc
      call util_append(self%lvdotr, source%lvdotr, nold, nsrc, lsource_mask)
      call util_append(self%status, source%status, nold, nsrc, lsource_mask)
      call util_append(self%index1, source%index1, nold, nsrc, lsource_mask)
      call util_append(self%index2, source%index2, nold, nsrc, lsource_mask)
      call util_append(self%id1, source%id1, nold, nsrc, lsource_mask)
      call util_append(self%id2, source%id2, nold, nsrc, lsource_mask)
      call util_append(self%x1, source%x1, nold, nsrc, lsource_mask)
      call util_append(self%x2, source%x2, nold, nsrc, lsource_mask)
      call util_append(self%v1, source%v1, nold, nsrc, lsource_mask)
      call util_append(self%v2, source%v2, nold, nsrc, lsource_mask)
      self%nenc = nold + count(lsource_mask(1:nsrc))

      return
   end subroutine encounter_util_append_list


   module subroutine encounter_util_copy_list(self, source)
      !! author: David A. Minton
      !!
      !! Copies elements from the source encounter list into self.
      implicit none
      ! Arguments
      class(encounter_list), intent(inout) :: self   !! Encounter list 
      class(encounter_list), intent(in)    :: source !! Source object to copy into

      associate(n => source%nenc)
         self%nenc = n
         self%t = source%t
         self%lvdotr(1:n) = source%lvdotr(1:n) 
         self%status(1:n) = source%status(1:n) 
         self%index1(1:n) = source%index1(1:n)
         self%index2(1:n) = source%index2(1:n)
         self%id1(1:n) = source%id1(1:n)
         self%id2(1:n) = source%id2(1:n)
         self%x1(:,1:n) = source%x1(:,1:n)
         self%x2(:,1:n) = source%x2(:,1:n)
         self%v1(:,1:n) = source%v1(:,1:n)
         self%v2(:,1:n) = source%v2(:,1:n)
      end associate

      return
   end subroutine encounter_util_copy_list

   module subroutine encounter_util_dealloc_aabb(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatables
      implicit none
      ! Arguments
      class(encounter_bounding_box_1D), intent(inout) :: self

      if (allocated(self%ind)) deallocate(self%ind)
      if (allocated(self%ibeg)) deallocate(self%ibeg)
      if (allocated(self%iend)) deallocate(self%iend)

      return
   end subroutine encounter_util_dealloc_aabb


   module subroutine encounter_util_dealloc_list(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatables
      implicit none
      ! Arguments
      class(encounter_list), intent(inout) :: self

      if (allocated(self%lvdotr)) deallocate(self%lvdotr)
      if (allocated(self%status)) deallocate(self%status)
      if (allocated(self%index1)) deallocate(self%index1)
      if (allocated(self%index2)) deallocate(self%index2)
      if (allocated(self%id1)) deallocate(self%id1)
      if (allocated(self%id2)) deallocate(self%id2)
      if (allocated(self%x1)) deallocate(self%x1)
      if (allocated(self%x2)) deallocate(self%x2)
      if (allocated(self%v1)) deallocate(self%v1)
      if (allocated(self%v2)) deallocate(self%v2)

      return
   end subroutine encounter_util_dealloc_list


   module subroutine encounter_util_final_aabb(self)
      !! author: David A. Minton
      !!
      !! Finalize the axis aligned bounding box (1D) - deallocates all allocatables
      implicit none
      ! Arguments
      type(encounter_bounding_box_1D), intent(inout) :: self

      call self%dealloc()

      return
   end subroutine encounter_util_final_aabb


   module subroutine encounter_util_final_list(self)
      !! author: David A. Minton
      !!
      !! Finalize the encounter list - deallocates all allocatables
      implicit none
      ! Arguments
      type(encounter_list), intent(inout) :: self

      call self%dealloc()

      return
   end subroutine encounter_util_final_list


   module subroutine encounter_util_resize_list(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of the encounter list against the required size and extends it by a factor of 2 more than requested if it is too small.
      !! Note: The reason to extend it by a factor of 2 is for performance. When there are many enounters per step, resizing every time you want to add an 
      !! encounter takes significant computational effort. Resizing by a factor of 2 is a tradeoff between performance (fewer resize calls) and memory managment
      !! Memory usage grows by a factor of 2 each time it fills up, but no more. 
      implicit none
      ! Arguments
      class(encounter_list), intent(inout) :: self !! Swiftest encounter list 
      integer(I8B),          intent(in)    :: nnew !! New size of list needed
      ! Internals
      class(encounter_list), allocatable :: enc_temp
      integer(I8B)                       :: nold
      logical                            :: lmalloc

      lmalloc = allocated(self%status)
      if (lmalloc) then
         nold = size(self%status)
      else
         nold = 0_I8B
      end if
      if (nnew > nold) then
         if (lmalloc) allocate(enc_temp, source=self)
         call self%setup(2_I8B * nnew)
         if (lmalloc) then
            call self%copy(enc_temp)
            deallocate(enc_temp)
         end if
      else
         self%status(nnew+1_I8B:nold) = INACTIVE
      end if
      self%nenc = nnew

      return
   end subroutine encounter_util_resize_list


   module subroutine encounter_util_spill_list(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest encounter structure from active list to discard list
      implicit none
      ! Arguments
      class(encounter_list), intent(inout) :: self         !! Swiftest encounter list 
      class(encounter_list), intent(inout) :: discards     !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
      ! Internals
      integer(I8B) :: nenc_old
  
      associate(keeps => self)
         call util_spill(keeps%lvdotr, discards%lvdotr, lspill_list, ldestructive)
         call util_spill(keeps%status, discards%status, lspill_list, ldestructive)
         call util_spill(keeps%index1, discards%index1, lspill_list, ldestructive)
         call util_spill(keeps%index2, discards%index2, lspill_list, ldestructive)
         call util_spill(keeps%id1, discards%id1, lspill_list, ldestructive)
         call util_spill(keeps%id2, discards%id2, lspill_list, ldestructive)
         call util_spill(keeps%x1, discards%x1, lspill_list, ldestructive)
         call util_spill(keeps%x2, discards%x2, lspill_list, ldestructive)
         call util_spill(keeps%v1, discards%v1, lspill_list, ldestructive)
         call util_spill(keeps%v2, discards%v2, lspill_list, ldestructive)

         nenc_old = keeps%nenc

         ! This is the base class, so will be the last to be called in the cascade. 
         ! Therefore we need to set the nenc values for both the keeps and discareds
         discards%nenc = count(lspill_list(1:nenc_old))
         if (ldestructive) keeps%nenc = nenc_old - discards%nenc
      end associate
   
      return
   end subroutine encounter_util_spill_list

end submodule s_encounter_util