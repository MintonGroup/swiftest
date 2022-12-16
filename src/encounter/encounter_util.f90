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
      call util_append(self%lclosest, source%lclosest, nold, nsrc, lsource_mask)
      call util_append(self%status, source%status, nold, nsrc, lsource_mask)
      call util_append(self%index1, source%index1, nold, nsrc, lsource_mask)
      call util_append(self%index2, source%index2, nold, nsrc, lsource_mask)
      call util_append(self%id1, source%id1, nold, nsrc, lsource_mask)
      call util_append(self%id2, source%id2, nold, nsrc, lsource_mask)
      call util_append(self%r1, source%r1, nold, nsrc, lsource_mask)
      call util_append(self%r2, source%r2, nold, nsrc, lsource_mask)
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
         self%lclosest(1:n) = source%lclosest(1:n) 
         self%status(1:n) = source%status(1:n) 
         self%index1(1:n) = source%index1(1:n)
         self%index2(1:n) = source%index2(1:n)
         self%id1(1:n) = source%id1(1:n)
         self%id2(1:n) = source%id2(1:n)
         self%r1(:,1:n) = source%r1(:,1:n)
         self%r2(:,1:n) = source%r2(:,1:n)
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
      if (allocated(self%lclosest)) deallocate(self%lclosest)
      if (allocated(self%status)) deallocate(self%status)
      if (allocated(self%index1)) deallocate(self%index1)
      if (allocated(self%index2)) deallocate(self%index2)
      if (allocated(self%id1)) deallocate(self%id1)
      if (allocated(self%id2)) deallocate(self%id2)
      if (allocated(self%r1)) deallocate(self%r1)
      if (allocated(self%r2)) deallocate(self%r2)
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


   module subroutine encounter_util_final_snapshot(self)
      !! author: David A. Minton
      !!
      !! Deallocates allocatable arrays in an encounter snapshot
      implicit none
      ! Arguments
      type(encounter_snapshot),  intent(inout) :: self !! Encounter storage object

      if (allocated(self%pl)) deallocate(self%pl)
      if (allocated(self%tp)) deallocate(self%tp)
      self%t = 0.0_DP

      return
   end subroutine encounter_util_final_snapshot


   module subroutine encounter_util_final_storage(self)
      !! author: David A. Minton
      !!
      !! Deallocates allocatable arrays in an encounter snapshot
      implicit none
      ! Arguments
      type(encounter_storage(*)),  intent(inout) :: self !! Encounter storage object

      call util_final_storage(self%swiftest_storage)

      return
   end subroutine encounter_util_final_storage


   module subroutine encounter_util_get_idvalues_snapshot(self, idvals)
      !! author: David A. Minton
      !!
      !! Returns an array of all id values saved in this snapshot
      implicit none
      ! Arguments
      class(encounter_snapshot),               intent(in)  :: self   !! Encounter snapshot object
      integer(I4B), dimension(:), allocatable, intent(out) :: idvals !! Array of all id values saved in this snapshot
      ! Internals
      integer(I4B) :: npl, ntp

      if (allocated(self%pl)) then
         npl = self%pl%nbody
      else
         npl = 0
      end if 
      if (allocated(self%tp)) then
         ntp = self%tp%nbody
      else
         ntp = 0
      end if

      if (npl + ntp == 0) return
      allocate(idvals(npl+ntp))

      if (npl > 0) idvals(1:npl) = self%pl%id(:)
      if (ntp >0) idvals(npl+1:npl+ntp) = self%tp%id(:)

      return

   end subroutine encounter_util_get_idvalues_snapshot


   subroutine encounter_util_get_vals_storage(storage, idvals, tvals)
      !! author: David A. Minton
      !!
      !! Gets the id values in a storage object, regardless of whether it is encounter of collision
      ! Argument
      class(swiftest_storage(*)), intent(in)              :: storage !! Swiftest storage object
      integer(I4B), dimension(:), allocatable, intent(out) :: idvals  !! Array of all id values in all snapshots
      real(DP),     dimension(:), allocatable, intent(out) :: tvals   !! Array of all time values in all snapshots
      ! Internals
      integer(I4B) :: i, n, nlo, nhi, ntotal
      integer(I4B), dimension(:), allocatable :: itmp

      associate(nsnaps => storage%iframe)

         allocate(tvals(nsnaps))

         tvals(:) = 0.0_DP

         ! First pass to get total number of ids
         ntotal = 0
         do i = 1, nsnaps
            if (allocated(storage%frame(i)%item)) then
               select type(snapshot => storage%frame(i)%item)
               class is (encounter_snapshot)
                  tvals(i) = snapshot%t
                  call snapshot%get_idvals(itmp)
                  if (allocated(itmp)) then
                     n = size(itmp)
                     ntotal = ntotal + n
                  end if
               end select
            end if
         end do

         allocate(idvals(ntotal))
         nlo = 1
         ! Second pass to store all ids get all of the ids stored
         do i = 1, nsnaps
            if (allocated(storage%frame(i)%item)) then
               select type(snapshot => storage%frame(i)%item)
               class is (encounter_snapshot)
                  tvals(i) = snapshot%t
                  call snapshot%get_idvals(itmp)
                  if (allocated(itmp)) then
                     n = size(itmp)
                     nhi = nlo + n - 1
                     idvals(nlo:nhi) = itmp(1:n)
                     nlo = nhi + 1 
                  end if
               end select
            end if
         end do

      end associate 
      return
   end subroutine encounter_util_get_vals_storage


   module subroutine encounter_util_index_map(self)
      !! author: David A. Minton
      !!
      !! Maps body id values to storage index values so we don't have to use unlimited dimensions for id.
      !! Basically this will make a unique list of ids that exist in all of the saved snapshots
      implicit none
      ! Arguments
      class(encounter_storage(*)), intent(inout) :: self !! Swiftest storage object
      ! Internals
      integer(I4B), dimension(:), allocatable :: idvals
      real(DP), dimension(:), allocatable :: tvals

      call encounter_util_get_vals_storage(self, idvals, tvals)

      ! Consolidate ids to only unique values
      call util_unique(idvals,self%idvals,self%idmap)
      self%nid = size(self%idvals)

      ! Consolidate time values to only unique values
      call util_unique(tvals,self%tvals,self%tmap)
      self%nt = size(self%tvals)

      return
   end subroutine encounter_util_index_map


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
         call util_spill(keeps%lclosest, discards%lclosest, lspill_list, ldestructive)
         call util_spill(keeps%status, discards%status, lspill_list, ldestructive)
         call util_spill(keeps%index1, discards%index1, lspill_list, ldestructive)
         call util_spill(keeps%index2, discards%index2, lspill_list, ldestructive)
         call util_spill(keeps%id1, discards%id1, lspill_list, ldestructive)
         call util_spill(keeps%id2, discards%id2, lspill_list, ldestructive)
         call util_spill(keeps%r1, discards%r1, lspill_list, ldestructive)
         call util_spill(keeps%r2, discards%r2, lspill_list, ldestructive)
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


   subroutine encounter_util_save_snapshot(encounter_history, snapshot)
      !! author: David A. Minton
      !!
      !! Checks the current size of the encounter storage against the required size and extends it by a factor of 2 more than requested if it is too small.
      !! Note: The reason to extend it by a factor of 2 is for performance. When there are many enounters per step, resizing every time you want to add an 
      !! encounter takes significant computational effort. Resizing by a factor of 2 is a tradeoff between performance (fewer resize calls) and memory managment
      !! Memory usage grows by a factor of 2 each time it fills up, but no more. 
      implicit none
      ! Arguments
      type(encounter_storage(*)), allocatable, intent(inout) :: encounter_history !! SyMBA encounter storage object
      class(encounter_snapshot),               intent(in)    :: snapshot          !! Encounter snapshot object
      ! Internals
      type(encounter_storage(nframes=:)), allocatable :: tmp
      integer(I4B) :: i, nnew, nold, nbig

      ! Advance the snapshot frame counter
      encounter_history%iframe = encounter_history%iframe + 1

      ! Check to make sure the current encounter_history object is big enough. If not, grow it by a factor of 2
      nnew = encounter_history%iframe
      nold = encounter_history%nframes

      if (nnew > nold) then
         nbig = nold
         do while (nbig < nnew)
            nbig = nbig * 2
         end do
         allocate(encounter_storage(nbig) :: tmp) 
         tmp%iframe = encounter_history%iframe
         call move_alloc(encounter_history%nc, tmp%nc)

         do i = 1, nold
            if (allocated(encounter_history%frame(i)%item)) call move_alloc(encounter_history%frame(i)%item, tmp%frame(i)%item)
         end do
         deallocate(encounter_history)
         call move_alloc(tmp,encounter_history)
         nnew = nbig
      end if

      ! Find out which time slot this belongs in by searching for an existing slot
      ! with the same value of time or the first available one
      encounter_history%frame(nnew) = snapshot

      return
   end subroutine encounter_util_save_snapshot


end submodule s_encounter_util