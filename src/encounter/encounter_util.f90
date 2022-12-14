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


   module subroutine encounter_util_final_collision_storage(self)
      !! author: David A. Minton
      !!
      !! Finalizer will deallocate all allocatables
      implicit none
      ! Arguments
      type(collision_storage(*)),  intent(inout) :: self !! SyMBA nbody system object

      call util_final_storage(self%swiftest_storage)

      return
   end subroutine encounter_util_final_collision_storage


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


   module subroutine encounter_util_index_map_encounter(self)
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
   end subroutine encounter_util_index_map_encounter



   module subroutine encounter_util_index_map_collision(self)
      !! author: David A. Minton
      !!
      !! Maps body id values to storage index values so we don't have to use unlimited dimensions for id
      implicit none
      ! Arguments
      class(collision_storage(*)), intent(inout) :: self  !! Swiftest storage object
      ! Internals
      integer(I4B), dimension(:), allocatable :: idvals
      real(DP), dimension(:), allocatable :: tvals

      call encounter_util_get_vals_storage(self, idvals, tvals)

      ! Consolidate ids to only unique values
      call util_unique(idvals,self%idvals,self%idmap)
      self%nid = size(self%idvals)

      ! Don't consolidate time values (multiple collisions can happen in a single time step)
      self%nt = size(self%tvals)

      return
   end subroutine encounter_util_index_map_collision


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



   subroutine encounter_util_save_collision(collision_history, snapshot)
      !! author: David A. Minton
      !!
      !! Checks the current size of the encounter storage against the required size and extends it by a factor of 2 more than requested if it is too small.
      !! Note: The reason to extend it by a factor of 2 is for performance. When there are many enounters per step, resizing every time you want to add an 
      !! encounter takes significant computational effort. Resizing by a factor of 2 is a tradeoff between performance (fewer resize calls) and memory managment
      !! Memory usage grows by a factor of 2 each time it fills up, but no more. 
      implicit none
      ! Arguments
      type(collision_storage(*)), allocatable, intent(inout) :: collision_history  !! Collision history object
      class(encounter_snapshot),               intent(in)    :: snapshot           !! Encounter snapshot object
      ! Internals
      type(collision_storage(nframes=:)), allocatable :: tmp
      integer(I4B) :: i, nnew, nold, nbig

      ! Advance the snapshot frame counter
      collision_history%iframe = collision_history%iframe + 1

      ! Check to make sure the current encounter_history object is big enough. If not, grow it by a factor of 2
      nnew = collision_history%iframe
      nold = collision_history%nframes

      if (nnew > nold) then
         nbig = nold
         do while (nbig < nnew)
            nbig = nbig * 2
         end do
         allocate(collision_storage(nbig) :: tmp) 
         tmp%iframe = collision_history%iframe
         call move_alloc(collision_history%nc, tmp%nc)

         do i = 1, nold
            if (allocated(collision_history%frame(i)%item)) call move_alloc(collision_history%frame(i)%item, tmp%frame(i)%item)
         end do
         deallocate(collision_history)
         call move_alloc(tmp,collision_history)
         nnew = nbig
      end if

      collision_history%frame(nnew) = snapshot

      return
   end subroutine encounter_util_save_collision


   subroutine encounter_util_save_encounter(encounter_history, snapshot)
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
   end subroutine encounter_util_save_encounter


   module subroutine encounter_util_snapshot_collision(self, param, system, t, arg)
      !! author: David A. Minton
      !!
      !! Takes a minimal snapshot of the state of the system during an encounter so that the trajectories
      !! can be played back through the encounter
      implicit none
      ! Internals
      class(collision_storage(*)),  intent(inout)        :: self   !! Swiftest storage object
      class(swiftest_parameters),   intent(inout)        :: param  !! Current run configuration parameters
      class(swiftest_nbody_system), intent(inout)        :: system !! Swiftest nbody system object to store
      real(DP),                     intent(in), optional :: t      !! Time of snapshot if different from system time
      character(*),                 intent(in), optional :: arg    !! "before": takes a snapshot just before the collision. "after" takes the snapshot just after the collision.
      ! Arguments
      class(fraggle_snapshot), allocatable :: snapshot
      type(symba_pl)                                 :: pl
      character(len=:), allocatable :: stage

      if (present(arg)) then
         stage = arg
      else
         stage = ""
      end if 

      select type (system)
      class is (symba_nbody_system)

         select case(stage)
         case("before")
            ! Saves the states of the bodies involved in the collision before the collision is resolved
            associate (idx => system%colliders%idx, ncoll => system%colliders%ncoll)
               call pl%setup(ncoll, param)
               pl%id(:) = system%pl%id(idx(:))
               pl%Gmass(:) = system%pl%Gmass(idx(:))
               pl%radius(:) = system%pl%radius(idx(:))
               pl%rot(:,:) = system%pl%rot(:,idx(:))
               pl%Ip(:,:) = system%pl%Ip(:,idx(:))
               pl%rh(:,:) = system%pl%rh(:,idx(:))
               pl%vh(:,:) = system%pl%vh(:,idx(:))
               pl%info(:) = system%pl%info(idx(:))
               !end select
               allocate(system%colliders%pl, source=pl)
            end associate
         case("after")
            allocate(fraggle_snapshot :: snapshot)
            allocate(snapshot%colliders, source=system%colliders) 
            allocate(snapshot%fragments, source=system%fragments)
            snapshot%t = t
            select type(param)
            class is (symba_parameters)
               call encounter_util_save_collision(param%collision_history,snapshot)
            end select
         case default
            write(*,*) "encounter_util_snapshot_collision requies either 'before' or 'after' passed to 'arg'"
         end select

      end select

      return
   end subroutine encounter_util_snapshot_collision


   module subroutine encounter_util_snapshot_encounter(self, param, system, t, arg)
      !! author: David A. Minton
      !!
      !! Takes a minimal snapshot of the state of the system during an encounter so that the trajectories
      !! can be played back through the encounter
      implicit none
      ! Internals
      class(encounter_storage(*)),  intent(inout)        :: self   !! Swiftest storage object
      class(swiftest_parameters),   intent(inout)        :: param  !! Current run configuration parameters
      class(swiftest_nbody_system), intent(inout)        :: system !! Swiftest nbody system object to store
      real(DP),                     intent(in), optional :: t      !! Time of snapshot if different from system time
      character(*),                 intent(in), optional :: arg    !! Optional argument (needed for extended storage type used in collision snapshots)
      ! Arguments
      class(encounter_snapshot), allocatable :: snapshot
      integer(I4B) :: i, pi, pj, k, npl_snap, ntp_snap, iflag
      real(DP), dimension(NDIM) :: rrel, vrel, rcom, vcom
      real(DP) :: Gmtot, a, q, capm, tperi
      real(DP), dimension(NDIM,2) :: rb,vb

      if (.not.present(t)) then
         write(*,*) "encounter_util_snapshot_encounter requires `t` to be passed"
         return
      end if

      if (.not.present(arg)) then
         write(*,*) "encounter_util_snapshot_encounter requires `arg` to be passed"
         return
      end if

      select type(param)
      class is (symba_parameters)
         select type (system)
         class is (symba_nbody_system)
            select type(pl => system%pl)
            class is (symba_pl)
               select type (tp => system%tp)
               class is (symba_tp)
                  associate(npl => pl%nbody,  ntp => tp%nbody)
                     if (npl + ntp == 0) return
                     allocate(encounter_snapshot :: snapshot)
                     allocate(snapshot%pl, mold=pl)
                     allocate(snapshot%tp, mold=tp)
                     snapshot%iloop = param%iloop

                     select type(pl_snap => snapshot%pl)
                     class is (symba_pl)
                        select type(tp_snap => snapshot%tp)
                        class is (symba_tp)

                           select case(arg)
                           case("trajectory")
                              snapshot%t = t
         
                              npl_snap = npl
                              ntp_snap = ntp
   
                              if (npl > 0) then
                                 pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE .and. pl%levelg(1:npl) == system%irec
                                 npl_snap = count(pl%lmask(1:npl))
                              end if
                              if (ntp > 0) then
                                 tp%lmask(1:ntp) = tp%status(1:ntp) /= INACTIVE .and. tp%levelg(1:ntp) == system%irec
                                 ntp_snap = count(tp%lmask(1:ntp))
                              end if

                              if (npl_snap + ntp_snap == 0) return ! Nothing to snapshot

                              pl_snap%nbody = npl_snap

                              ! Take snapshot of the currently encountering massive bodies
                              if (npl_snap > 0) then
                                 call pl_snap%setup(npl_snap, param)
                                 pl_snap%levelg(:) = pack(pl%levelg(1:npl), pl%lmask(1:npl))
                                 pl_snap%id(:) = pack(pl%id(1:npl), pl%lmask(1:npl))
                                 pl_snap%info(:) = pack(pl%info(1:npl), pl%lmask(1:npl))
                                 pl_snap%Gmass(:) = pack(pl%Gmass(1:npl), pl%lmask(1:npl))
                                 do i = 1, NDIM
                                    pl_snap%rh(i,:) = pack(pl%rh(i,1:npl), pl%lmask(1:npl))
                                    pl_snap%vh(i,:) = pack(pl%vb(i,1:npl), pl%lmask(1:npl))
                                 end do
                                 if (param%lclose) then
                                    pl_snap%radius(:) = pack(pl%radius(1:npl), pl%lmask(1:npl))
                                 end if

                                 if (param%lrotation) then
                                    do i = 1, NDIM
                                       pl_snap%Ip(i,:) = pack(pl%Ip(i,1:npl), pl%lmask(1:npl))
                                       pl_snap%rot(i,:) = pack(pl%rot(i,1:npl), pl%lmask(1:npl))
                                    end do
                                 end if
                                 call pl_snap%sort("id", ascending=.true.)
                              end if

                              ! Take snapshot of the currently encountering test particles
                              tp_snap%nbody = ntp_snap
                              if (ntp_snap > 0) then
                                 call tp_snap%setup(ntp_snap, param)
                                 tp_snap%id(:) = pack(tp%id(1:ntp), tp%lmask(1:ntp))
                                 tp_snap%info(:) = pack(tp%info(1:ntp), tp%lmask(1:ntp))
                                 do i = 1, NDIM
                                    tp_snap%rh(i,:) = pack(tp%rh(i,1:ntp), tp%lmask(1:ntp))
                                    tp_snap%vh(i,:) = pack(tp%vh(i,1:ntp), tp%lmask(1:ntp))
                                 end do
                              end if

                              ! Save the snapshot
                              param%encounter_history%nid = param%encounter_history%nid + ntp_snap + npl_snap
                              call encounter_util_save_encounter(param%encounter_history,snapshot)
                           case("closest")
                              associate(plplenc_list => system%plplenc_list, pltpenc_list => system%pltpenc_list)
                                 if (any(plplenc_list%lclosest(:))) then
                                    call pl_snap%setup(2, param)
                                    do k = 1, plplenc_list%nenc
                                       if (plplenc_list%lclosest(k)) then
                                          pi = plplenc_list%index1(k)
                                          pj = plplenc_list%index2(k)
                                          pl_snap%levelg(:) = pl%levelg([pi,pj])
                                          pl_snap%id(:) = pl%id([pi,pj])
                                          pl_snap%info(:) = pl%info([pi,pj])
                                          pl_snap%Gmass(:) = pl%Gmass([pi,pj])
                                          Gmtot = sum(pl_snap%Gmass(:))
                                          if (param%lclose) pl_snap%radius(:) = pl%radius([pi,pj])
                                          if (param%lrotation) then
                                             do i = 1, NDIM
                                                pl_snap%Ip(i,:) = pl%Ip(i,[pi,pj])
                                                pl_snap%rot(i,:) = pl%rot(i,[pi,pj])
                                             end do
                                          end if

                                          ! Compute pericenter passage time to get the closest approach parameters
                                          rrel(:) = plplenc_list%r2(:,k) - plplenc_list%r1(:,k)
                                          vrel(:) = plplenc_list%v2(:,k) - plplenc_list%v1(:,k)
                                          call orbel_xv2aqt(Gmtot, rrel(1), rrel(2), rrel(3), vrel(1), vrel(2), vrel(3), a, q, capm, tperi)
                                          snapshot%t = t + tperi
                                          if ((snapshot%t < maxval(pl_snap%info(:)%origin_time)) .or. &
                                              (snapshot%t > minval(pl_snap%info(:)%discard_time))) cycle

                                          ! Computer the center mass of the pair
                                          rcom(:) = (plplenc_list%r1(:,k) * pl_snap%Gmass(1) + plplenc_list%r2(:,k) * pl_snap%Gmass(2)) / Gmtot
                                          vcom(:) = (plplenc_list%v1(:,k) * pl_snap%Gmass(1) + plplenc_list%v2(:,k) * pl_snap%Gmass(2)) / Gmtot
                                          rb(:,1) = plplenc_list%r1(:,k) - rcom(:)
                                          rb(:,2) = plplenc_list%r2(:,k) - rcom(:)
                                          vb(:,1) = plplenc_list%v1(:,k) - vcom(:)
                                          vb(:,2) = plplenc_list%v2(:,k) - vcom(:)

                                          ! Drift the relative orbit to get the new relative position and velocity
                                          call drift_one(Gmtot, rrel(1), rrel(2), rrel(3), vrel(1), vrel(2), vrel(3), tperi, iflag)
                                          if (iflag /= 0) write(*,*) "Danby error in encounter_util_snapshot_encounter. Closest approach positions and vectors may not be accurate."

                                          ! Get the new position and velocity vectors
                                          rb(:,1) = -(pl_snap%Gmass(2) / Gmtot) * rrel(:)
                                          rb(:,2) =  (pl_snap%Gmass(1)) / Gmtot * rrel(:)

                                          vb(:,1) = -(pl_snap%Gmass(2) / Gmtot) * vrel(:)
                                          vb(:,2) =  (pl_snap%Gmass(1)) / Gmtot * vrel(:)

                                          ! Move the CoM assuming constant velocity over the time it takes to reach periapsis
                                          rcom(:) = rcom(:) + vcom(:) * tperi

                                          ! Compute the heliocentric position and velocity vector at periapsis
                                          pl_snap%rh(:,1) = rb(:,1) + rcom(:)
                                          pl_snap%rh(:,2) = rb(:,2) + rcom(:)
                                          pl_snap%vh(:,1) = vb(:,1) + vcom(:)
                                          pl_snap%vh(:,2) = vb(:,2) + vcom(:)

                                          call pl_snap%sort("id", ascending=.true.)
                                          call encounter_util_save_encounter(param%encounter_history,snapshot)
                                       end if
                                    end do

                                    plplenc_list%lclosest(:) = .false.
                                 end if

                                 if (any(pltpenc_list%lclosest(:))) then
                                    do k = 1, pltpenc_list%nenc
                                    end do
                                    pltpenc_list%lclosest(:) = .false.
                                 end if
                              end associate
                           case default
                              write(*,*) "encounter_util_snapshot_encounter requires `arg` to be either `trajectory` or `closest`"
                           end select
                        end select
                     end select
                  end associate
               end select
            end select
         end select
      end select

      return
   end subroutine encounter_util_snapshot_encounter

end submodule s_encounter_util