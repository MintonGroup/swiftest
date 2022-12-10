!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(symba_classes) s_symba_util
   use swiftest
contains


   module subroutine symba_util_append_arr_kin(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of kinship type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      type(symba_kinship), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      type(symba_kinship), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                                   intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,             dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine symba_util_append_arr_kin


   module subroutine symba_util_append_encounter_list(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one encounter list (pl-pl or pl-tp) body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(symba_encounter), intent(inout) :: self         !! SyMBA encounter list object
      class(encounter_list),  intent(in)    :: source       !! Source object to append
      logical, dimension(:),  intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nold, nsrc

      nold = self%nenc
      nsrc = source%nenc
      select type(source)
      class is (symba_encounter)
         call util_append(self%level, source%level, nold, nsrc, lsource_mask)
      end select
      call encounter_util_append_list(self, source, lsource_mask) 

      return
   end subroutine symba_util_append_encounter_list


   module subroutine symba_util_append_pl(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one massive body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      !! Arguments
      class(symba_pl),                 intent(inout) :: self         !! SyMBA massive body object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (symba_pl)
         associate(nold => self%nbody, nsrc => source%nbody)
            call util_append(self%lcollision, source%lcollision, nold, nsrc, lsource_mask)
            call util_append(self%lencounter, source%lencounter, nold, nsrc, lsource_mask)
            call util_append(self%lmtiny, source%lmtiny, nold, nsrc, lsource_mask)
            call util_append(self%nplenc, source%nplenc, nold, nsrc, lsource_mask)
            call util_append(self%ntpenc, source%ntpenc, nold, nsrc, lsource_mask)
            call util_append(self%levelg, source%levelg, nold, nsrc, lsource_mask)
            call util_append(self%levelm, source%levelm, nold, nsrc, lsource_mask)
            call util_append(self%isperi, source%isperi, nold, nsrc, lsource_mask)
            call util_append(self%peri, source%peri, nold, nsrc, lsource_mask)
            call util_append(self%atp, source%atp, nold, nsrc, lsource_mask)
            call util_append(self%kin, source%kin, nold, nsrc, lsource_mask)

            call util_append_pl(self, source, lsource_mask) ! Note: helio_pl does not have its own append method, so we skip back to the base class
         end associate
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class symba_pl or its descendents!"
         call util_exit(FAILURE)
      end select

      return
   end subroutine symba_util_append_pl


   module subroutine symba_util_append_merger(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one massive body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(symba_merger),             intent(inout) :: self         !! SyMBA massive body object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B), dimension(:), allocatable        :: ncomp_tmp    !! Temporary placeholder for ncomp incase we are appending a symba_pl object to a symba_merger
      integer(I4B) :: nold, nsrc, nnew

      nold = self%nbody
      nsrc = source%nbody
      nnew = count(lsource_mask)

      select type(source)
      class is (symba_merger)
         call util_append(self%ncomp, source%ncomp, nold, nsrc, lsource_mask)
         call symba_util_append_pl(self, source, lsource_mask) 
      class is (symba_pl)
         allocate(ncomp_tmp, mold=source%id)
         ncomp_tmp(:) = 0
         call util_append(self%ncomp, ncomp_tmp, nold, nsrc, lsource_mask)
         call symba_util_append_pl(self, source, lsource_mask) 
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class symba_pl or its descendents!"
         call util_exit(FAILURE)
      end select

      ! Save the number of appended bodies 
      self%ncomp(nold+1:nold+nnew) = nnew

      return
   end subroutine symba_util_append_merger


   module subroutine symba_util_append_tp(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from test particle object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      !! Arguments
      class(symba_tp),                 intent(inout) :: self         !! SyMBA test particle object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (symba_tp)
         associate(nold => self%nbody, nsrc => source%nbody)
            call util_append(self%nplenc, source%nplenc, nold, nsrc, lsource_mask)
            call util_append(self%levelg, source%levelg, nold, nsrc, lsource_mask)
            call util_append(self%levelm, source%levelm, nold, nsrc, lsource_mask)

            call util_append_tp(self, source, lsource_mask) ! Note: helio_tp does not have its own append method, so we skip back to the base class
         end associate
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class symba_tp or its descendents!"
         call util_exit(FAILURE)
      end select

      return
   end subroutine symba_util_append_tp


   module subroutine symba_util_copy_encounter_list(self, source)
      !! author: David A. Minton
      !!
      !! Copies elements from the source encounter list into self.
      implicit none
      ! Arguments
      class(symba_encounter),    intent(inout) :: self   !! Encounter list 
      class(encounter_list), intent(in)    :: source !! Source object to copy into
  
      select type(source)
      class is (symba_encounter)
         associate(n => source%nenc)
            self%level(1:n) = source%level(1:n) 
            self%tcollision(1:n) = source%tcollision(1:n) 
         end associate
      end select

      call encounter_util_copy_list(self, source)
   
      return
   end subroutine symba_util_copy_encounter_list


   module subroutine symba_util_dealloc_encounter_list(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Argumentse
      class(symba_encounter),  intent(inout) :: self !! SyMBA encounter list

      if (allocated(self%level)) deallocate(self%level)
      if (allocated(self%tcollision)) deallocate(self%tcollision)

      return
   end subroutine symba_util_dealloc_encounter_list


   module subroutine symba_util_dealloc_kin(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Arguments
      class(symba_kinship),  intent(inout) :: self !! SyMBA kinship object

      if (allocated(self%child)) deallocate(self%child)

      return
   end subroutine symba_util_dealloc_kin


   module subroutine symba_util_dealloc_merger(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Arguments
      class(symba_merger),  intent(inout) :: self !! SyMBA body merger object

      if (allocated(self%ncomp)) deallocate(self%ncomp)

      call symba_util_dealloc_pl(self)

      return
   end subroutine symba_util_dealloc_merger


   module subroutine symba_util_dealloc_pl(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Arguments
      class(symba_pl),  intent(inout) :: self !! SyMBA massive body object
      ! Internals
      integer(I4B) :: i

      if (allocated(self%lcollision)) deallocate(self%lcollision)
      if (allocated(self%lencounter)) deallocate(self%lencounter)
      if (allocated(self%lmtiny)) deallocate(self%lmtiny)
      if (allocated(self%nplenc)) deallocate(self%nplenc)
      if (allocated(self%ntpenc)) deallocate(self%ntpenc)
      if (allocated(self%levelg)) deallocate(self%levelg)
      if (allocated(self%levelm)) deallocate(self%levelm)
      if (allocated(self%isperi)) deallocate(self%isperi)
      if (allocated(self%peri)) deallocate(self%peri)
      if (allocated(self%atp)) deallocate(self%atp)

      if (allocated(self%kin)) then
         do i = 1, self%nbody
            call self%kin(i)%dealloc()
         end do
         deallocate(self%kin)
      end if

      call util_dealloc_pl(self)

      return
   end subroutine symba_util_dealloc_pl


   module subroutine symba_util_dealloc_tp(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Arguments
      class(symba_tp),  intent(inout) :: self !! SyMBA test particle object

      if (allocated(self%nplenc)) deallocate(self%nplenc)
      if (allocated(self%levelg)) deallocate(self%levelg)
      if (allocated(self%levelm)) deallocate(self%levelm)

      call util_dealloc_tp(self)

      return
   end subroutine symba_util_dealloc_tp


   module subroutine symba_util_fill_arr_kin(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of particle kinship types
      !! This is the inverse of a spill operation   
      implicit none
      ! Arguments
      type(symba_kinship), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      type(symba_kinship), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,             dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
      keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))
   
      return
   end subroutine symba_util_fill_arr_kin


   module subroutine symba_util_fill_pl(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new SyMBA test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      implicit none
      ! Arguments
      class(symba_pl),       intent(inout) :: self       !! SyMBA masive body object
      class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)
         select type(inserts)
         class is (symba_pl)
            call util_fill(keeps%lcollision, inserts%lcollision, lfill_list)
            call util_fill(keeps%lencounter, inserts%lencounter, lfill_list)
            call util_fill(keeps%lmtiny, inserts%lmtiny, lfill_list)
            call util_fill(keeps%nplenc, inserts%nplenc, lfill_list)
            call util_fill(keeps%ntpenc, inserts%ntpenc, lfill_list)
            call util_fill(keeps%levelg, inserts%levelg, lfill_list)
            call util_fill(keeps%levelm, inserts%levelm, lfill_list)
            call util_fill(keeps%isperi, inserts%isperi, lfill_list)
            call util_fill(keeps%peri, inserts%peri, lfill_list)
            call util_fill(keeps%atp, inserts%atp, lfill_list)
            call util_fill(keeps%kin, inserts%kin, lfill_list)
            
            call util_fill_pl(keeps, inserts, lfill_list)  ! Note: helio_pl does not have its own fill method, so we skip back to the base class
         class default
            write(*,*) "Invalid object passed to the fill method. Source must be of class symba_pl or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate

      return
   end subroutine symba_util_fill_pl


   module subroutine symba_util_fill_tp(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new SyMBA test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      implicit none
      ! Arguments
      class(symba_tp),       intent(inout) :: self       !! SyMBA test particle object
      class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)
         select type(inserts)
         class is (symba_tp)
            call util_fill(keeps%nplenc, inserts%nplenc, lfill_list)
            call util_fill(keeps%levelg, inserts%levelg, lfill_list)
            call util_fill(keeps%levelm, inserts%levelm, lfill_list)
            
            call util_fill_tp(keeps, inserts, lfill_list) ! Note: helio_tp does not have its own fill method, so we skip back to the base class
         class default
            write(*,*) "Invalid object passed to the fill method. Source must be of class symba_tp or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate

      return
   end subroutine symba_util_fill_tp


   module subroutine symba_util_flatten_eucl_plpl(self, param)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix. This also sets the lmtiny flag and computes the
      !! number of interactions that excludes semi-interacting bodies with each other (Gmass < GMTINY).
      !! This method will also sort the bodies in descending order by Mass
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751
      implicit none
      ! Arguments
      class(symba_pl),            intent(inout) :: self  !! SyMBA massive body object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I8B) :: npl, nplm

      associate(pl => self, nplplm => self%nplplm)
         npl = int(self%nbody, kind=I8B)
         select type(param)
         class is (symba_parameters)
            pl%lmtiny(1:npl) = pl%Gmass(1:npl) < param%GMTINY 
         end select
         nplm = count(.not. pl%lmtiny(1:npl))
         pl%nplm = int(nplm, kind=I4B)
         nplplm = nplm * npl - nplm * (nplm + 1_I8B) / 2_I8B ! number of entries in a strict lower triangle, npl x npl, minus first column including only mutually interacting bodies

         call util_flatten_eucl_plpl(pl, param)
      end associate

      return
   end subroutine symba_util_flatten_eucl_plpl


   module subroutine symba_util_final_encounter_list(self)
      !! author: David A. Minton
      !!
      !! Finalize the SyMBA encounter list object - deallocates all allocatables
      implicit none
      ! Argument
      type(symba_encounter),  intent(inout) :: self !! SyMBA encounter list object

      call self%dealloc()

      return
   end subroutine symba_util_final_encounter_list

   module subroutine symba_util_final_kin(self)
      !! author: David A. Minton
      !!
      !! Finalize the SyMBA kinship object - deallocates all allocatables
      implicit none
      ! Argument
      type(symba_kinship),  intent(inout) :: self !! SyMBA kinship object

      call self%dealloc()

      return
   end subroutine symba_util_final_kin

   module subroutine symba_util_final_merger(self)
      !! author: David A. Minton
      !!
      !! Finalize the SyMBA merger object - deallocates all allocatables
      implicit none
      ! Argument
      type(symba_merger),  intent(inout) :: self !! SyMBA merger object

      call self%dealloc()

      return
   end subroutine symba_util_final_merger


   module subroutine symba_util_final_pl(self)
      !! author: David A. Minton
      !!
      !! Finalize the SyMBA massive body object - deallocates all allocatables
      implicit none
      ! Argument
      type(symba_pl),  intent(inout) :: self !! SyMBA massive body object

      call self%dealloc()

      return
   end subroutine symba_util_final_pl


   module subroutine symba_util_final_system(self)
      !! author: David A. Minton
      !!
      !! Finalize the SyMBA nbody system object - deallocates all allocatables
      implicit none
      ! Argument
      type(symba_nbody_system),  intent(inout) :: self !! SyMBA nbody system object

      if (allocated(self%pl_adds)) deallocate(self%pl_adds)
      if (allocated(self%pltpenc_list)) deallocate(self%pltpenc_list)
      if (allocated(self%plplenc_list)) deallocate(self%plplenc_list)
      if (allocated(self%plplcollision_list)) deallocate(self%plplcollision_list)

      call helio_util_final_system(self%helio_nbody_system)

      return
   end subroutine symba_util_final_system


   module subroutine symba_util_final_tp(self)
      !! author: David A. Minton
      !!
      !! Finalize the SyMBA test particleobject - deallocates all allocatables
      implicit none
      ! Argument
      type(symba_tp),  intent(inout) :: self !! SyMBA test particle object

      call self%dealloc()

      return
   end subroutine symba_util_final_tp


   module subroutine symba_util_peri_pl(self, system, param)
      !! author: David A. Minton
      !!
      !! Determine system pericenter passages for planets in SyMBA
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_peri.f90
      !! Adapted from Hal Levison's Swift routine util_mass_peri.f
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! SyMBA massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)       :: i
      real(DP)           :: vdotr, e

      associate(pl => self, npl => self%nbody)
         if (pl%lfirst) then
            if (param%qmin_coord == "HELIO") then
               do i = 1, npl
                  if (pl%status(i) == ACTIVE) then
                     vdotr = dot_product(pl%rh(:,i), pl%vh(:,i))
                     if (vdotr > 0.0_DP) then
                        pl%isperi(i) = 1
                     else
                        pl%isperi(i) = -1
                     end if
                  end if
               end do
            else
               do i = 1, npl
                  if (pl%status(i) == ACTIVE) then
                     vdotr = dot_product(pl%xb(:,i), pl%vb(:,i))
                     if (vdotr > 0.0_DP) then
                        pl%isperi(i) = 1
                     else
                        pl%isperi(i) = -1
                     end if
                  end if
               end do
            end if
         else
            if (param%qmin_coord == "HELIO") then
               do i = 1, npl
                  if (pl%status(i) == ACTIVE) then
                     vdotr = dot_product(pl%rh(:,i), pl%vh(:,i))
                     if (pl%isperi(i) == -1) then
                        if (vdotr >= 0.0_DP) then
                           pl%isperi(i) = 0
                           CALL orbel_xv2aeq(pl%mu(i), pl%rh(1,i), pl%rh(2,i), pl%rh(3,i), pl%vh(1,i), pl%vh(2,i), pl%vh(3,i), &
                                  pl%atp(i), e, pl%peri(i))
                        end if
                     else
                        if (vdotr > 0.0_DP) then
                           pl%isperi(i) = 1
                        else
                           pl%isperi(i) = -1
                        end if
                     end if
                  end if
               end do
            else
               do i = 1, npl
                  if (pl%status(i) == ACTIVE) then
                     vdotr = dot_product(pl%xb(:,i), pl%vb(:,i))
                     if (pl%isperi(i) == -1) then
                        if (vdotr >= 0.0_DP) then
                           pl%isperi(i) = 0
                           CALL orbel_xv2aeq(system%Gmtot, pl%xb(1,i), pl%xb(2,i), pl%xb(3,i), pl%vb(1,i), pl%vb(2,i), pl%vb(3,i),&
                                             pl%atp(i), e, pl%peri(i))
                        end if
                     else
                        if (vdotr > 0.0_DP) then
                           pl%isperi(i) = 1
                        else
                           pl%isperi(i) = -1
                        end if
                     end if
                  end if
               end do
            end if
         end if
      end associate
 
     return
   end subroutine symba_util_peri_pl


   module subroutine symba_util_rearray_pl(self, system, param)
      !! Author: the Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Clean up the massive body structures to remove discarded bodies and add new bodies
      implicit none
      ! Arguments
      class(symba_pl),           intent(inout) :: self   !! SyMBA massive body object
      class(symba_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(symba_parameters),   intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      class(symba_pl), allocatable :: tmp !! The discarded body list.
      integer(I4B) :: i, k, npl, nadd, nencmin, nenc_old, idnew1, idnew2, idold1, idold2
      logical, dimension(:), allocatable :: lmask, ldump_mask
      class(symba_plplenc), allocatable :: plplenc_old
      logical :: lencounter
      integer(I4B), dimension(:), allocatable :: levelg_orig_pl, levelm_orig_pl, levelg_orig_tp, levelm_orig_tp
      integer(I4B), dimension(:), allocatable :: nplenc_orig_pl, nplenc_orig_tp, ntpenc_orig_pl

      associate(pl => self, pl_adds => system%pl_adds)

         npl = pl%nbody
         nadd = pl_adds%nbody
         if (npl == 0) return
         ! Deallocate any temporary variables
         if (allocated(pl%xbeg)) deallocate(pl%xbeg)
         if (allocated(pl%xend)) deallocate(pl%xend)

         ! Remove the discards and destroy the list, as the system already tracks pl_discards elsewhere
         allocate(lmask(npl))
         lmask(1:npl) = pl%ldiscard(1:npl)
         if (count(lmask(:)) > 0) then
            allocate(tmp, mold=self)
            call pl%spill(tmp, lspill_list=lmask, ldestructive=.true.)
            npl = pl%nbody
            call tmp%setup(0,param)
            deallocate(tmp)
            deallocate(lmask)
         end if

         ! Store the original plplenc list so we don't remove any of the original encounters
         nenc_old = system%plplenc_list%nenc
         if (nenc_old > 0) then 
            allocate(plplenc_old, source=system%plplenc_list)
            call plplenc_old%copy(system%plplenc_list)
         end if

         ! Add in any new bodies
         if (nadd > 0) then
            ! Append the adds to the main pl object
            call pl%append(pl_adds, lsource_mask=[(.true., i=1, nadd)])

            allocate(ldump_mask(npl+nadd))  ! This mask is used only to append the original Fortran binary particle.dat file with new bodies. This is ignored for NetCDF output
            ldump_mask(1:npl) = .false.
            ldump_mask(npl+1:npl+nadd) = pl%status(npl+1:npl+nadd) == NEW_PARTICLE
            npl = pl%nbody
         else
            allocate(ldump_mask(npl))
            ldump_mask(:) = .false.
         end if

         ! Reset all of the status flags for this body
         pl%status(1:npl) = ACTIVE
         do i = 1, npl
            call pl%info(i)%set_value(status="ACTIVE")
         end do
         pl%ldiscard(1:npl) = .false.
         pl%lcollision(1:npl) = .false.
         pl%lmask(1:npl) = .true.

         select type(param)
         class is (symba_parameters)
            pl%lmtiny(1:npl) = pl%Gmass(1:npl) < param%GMTINY
            where(pl%lmtiny(1:npl))
               pl%info(1:npl)%particle_type = PL_TINY_TYPE_NAME 
            elsewhere
               pl%info(1:npl)%particle_type = PL_TYPE_NAME 
            end where
         end select

         call pl%write_info(param%nc, param)
         deallocate(ldump_mask)

         ! Reindex the new list of bodies 
         call pl%sort("mass", ascending=.false.)
         call pl%flatten(param)

         ! Reset the kinship trackers
         call pl%reset_kinship([(i, i=1, npl)])

         ! Re-build the zero-level encounter list, being sure to save the original level information for all bodies
         allocate(levelg_orig_pl, source=pl%levelg)
         allocate(levelm_orig_pl, source=pl%levelm)
         allocate(nplenc_orig_pl, source=pl%nplenc)
         lencounter = pl%encounter_check(param, system, param%dt, 0) 
         if (system%tp%nbody > 0) then
            select type(tp => system%tp)
            class is (symba_tp)
               allocate(ntpenc_orig_pl, source=pl%ntpenc)
               allocate(levelg_orig_tp, source=tp%levelg)
               allocate(levelm_orig_tp, source=tp%levelm)
               allocate(nplenc_orig_tp, source=tp%nplenc)
               lencounter = tp%encounter_check(param, system, param%dt, 0)
               call move_alloc(levelg_orig_tp, tp%levelg)
               call move_alloc(levelm_orig_tp, tp%levelm)
               call move_alloc(nplenc_orig_tp, tp%nplenc)
               call move_alloc(ntpenc_orig_pl, pl%ntpenc)
            end select
         end if
         call move_alloc(levelg_orig_pl, pl%levelg)
         call move_alloc(levelm_orig_pl, pl%levelm)
         call move_alloc(nplenc_orig_pl, pl%nplenc)

         ! Re-index the encounter list as the index values may have changed
         if (nenc_old > 0) then
            nencmin = min(system%plplenc_list%nenc, plplenc_old%nenc) 
            system%plplenc_list%nenc = nencmin
            do k = 1, nencmin
               idnew1 = system%plplenc_list%id1(k)
               idnew2 = system%plplenc_list%id2(k)
               idold1 = plplenc_old%id1(k)
               idold2 = plplenc_old%id2(k)
               if ((idnew1 == idold1) .and. (idnew2 == idold2)) then
                  ! This is an encounter we already know about, so save the old information
                  system%plplenc_list%lvdotr(k) = plplenc_old%lvdotr(k) 
                  system%plplenc_list%status(k) = plplenc_old%status(k) 
                  system%plplenc_list%x1(:,k) = plplenc_old%x1(:,k)
                  system%plplenc_list%x2(:,k) = plplenc_old%x2(:,k)
                  system%plplenc_list%v1(:,k) = plplenc_old%v1(:,k)
                  system%plplenc_list%v2(:,k) = plplenc_old%v2(:,k)
                  system%plplenc_list%tcollision(k) = plplenc_old%tcollision(k)
                  system%plplenc_list%level(k) = plplenc_old%level(k)
               else if (((idnew1 == idold2) .and. (idnew2 == idold1))) then
                  ! This is an encounter we already know about, but with the order reversed, so save the old information
                  system%plplenc_list%lvdotr(k) = plplenc_old%lvdotr(k) 
                  system%plplenc_list%status(k) = plplenc_old%status(k) 
                  system%plplenc_list%x1(:,k) = plplenc_old%x2(:,k)
                  system%plplenc_list%x2(:,k) = plplenc_old%x1(:,k)
                  system%plplenc_list%v1(:,k) = plplenc_old%v2(:,k)
                  system%plplenc_list%v2(:,k) = plplenc_old%v1(:,k)
                  system%plplenc_list%tcollision(k) = plplenc_old%tcollision(k)
                  system%plplenc_list%level(k) = plplenc_old%level(k)
               end if
               system%plplenc_list%index1(k) = findloc(pl%id(1:npl), system%plplenc_list%id1(k), dim=1)
               system%plplenc_list%index2(k) = findloc(pl%id(1:npl), system%plplenc_list%id2(k), dim=1)
            end do
            if (allocated(lmask)) deallocate(lmask)
            allocate(lmask(nencmin))
            nenc_old = nencmin
            if (any(system%plplenc_list%index1(1:nencmin) == 0) .or. any(system%plplenc_list%index2(1:nencmin) == 0)) then
               lmask(:) = system%plplenc_list%index1(1:nencmin) /= 0 .and. system%plplenc_list%index2(1:nencmin) /= 0
            else
               return
            end if
            nencmin = count(lmask(:))
            system%plplenc_list%nenc = nencmin
            if (nencmin > 0) then
               system%plplenc_list%index1(1:nencmin) = pack(system%plplenc_list%index1(1:nenc_old), lmask(1:nenc_old))
               system%plplenc_list%index2(1:nencmin) = pack(system%plplenc_list%index2(1:nenc_old), lmask(1:nenc_old))
               system%plplenc_list%id1(1:nencmin) = pack(system%plplenc_list%id1(1:nenc_old), lmask(1:nenc_old))
               system%plplenc_list%id2(1:nencmin) = pack(system%plplenc_list%id2(1:nenc_old), lmask(1:nenc_old))
               system%plplenc_list%lvdotr(1:nencmin) = pack(system%plplenc_list%lvdotr(1:nenc_old), lmask(1:nenc_old))
               system%plplenc_list%status(1:nencmin) = pack(system%plplenc_list%status(1:nenc_old), lmask(1:nenc_old))
               system%plplenc_list%tcollision(1:nencmin) = pack(system%plplenc_list%tcollision(1:nenc_old), lmask(1:nenc_old))
               system%plplenc_list%level(1:nencmin) = pack(system%plplenc_list%level(1:nenc_old), lmask(1:nenc_old))
               do i = 1, NDIM
                  system%plplenc_list%x1(i, 1:nencmin) = pack(system%plplenc_list%x1(i, 1:nenc_old), lmask(1:nenc_old))
                  system%plplenc_list%x2(i, 1:nencmin) = pack(system%plplenc_list%x2(i, 1:nenc_old), lmask(1:nenc_old))
                  system%plplenc_list%v1(i, 1:nencmin) = pack(system%plplenc_list%v1(i, 1:nenc_old), lmask(1:nenc_old))
                  system%plplenc_list%v2(i, 1:nencmin) = pack(system%plplenc_list%v2(i, 1:nenc_old), lmask(1:nenc_old))
               end do
            end if
         end if
      end associate

      return
   end subroutine symba_util_rearray_pl


   module subroutine symba_util_reset_kinship(self, idx)
      !! author: David A. Minton
      !! 
      !! Resets the kinship status of bodies.
      !!
      implicit none
      class(symba_pl),            intent(inout) :: self !! SyMBA massive body object
      integer(I4B), dimension(:), intent(in)    :: idx  !! Index array of bodies to reset
      ! Internals
      integer(I4B) :: i, j

      self%kin(idx(:))%parent = idx(:)
      self%kin(idx(:))%nchild = 0
      do j = 1, size(idx(:))
         i = idx(j)
         if (allocated(self%kin(i)%child)) deallocate(self%kin(i)%child)
      end do

      return
   end subroutine symba_util_reset_kinship
   

   module subroutine symba_util_resize_arr_kin(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of type character string. Array will only be resized if has previously been allocated. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      type(symba_kinship), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                                   intent(in)    :: nnew !! New size
      ! Internals
      type(symba_kinship), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size

      if (nnew < 0) return

      if (nnew == 0) then
         if (allocated(arr)) deallocate(arr)
         return
      end if
      
      if (allocated(arr)) then
         nold = size(arr)
      else
         nold = 0
      end if

      allocate(tmp(nnew))
      if (nnew > nold) then
         tmp(1:nold) = arr(1:nold)
      else
         tmp(1:nnew) = arr(1:nnew)
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine symba_util_resize_arr_kin


   module subroutine symba_util_resize_merger(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a SyMBA merger list against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(symba_merger), intent(inout) :: self  !! SyMBA massive body object
      integer(I4B),        intent(in)    :: nnew  !! New size neded

      call util_resize(self%ncomp, nnew)

      call symba_util_resize_pl(self, nnew)

      return
   end subroutine symba_util_resize_merger


   module subroutine symba_util_resize_pl(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a SyMBA massive body object against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(symba_pl), intent(inout) :: self  !! SyMBA massive body object
      integer(I4B),    intent(in)    :: nnew  !! New size neded

      call util_resize(self%lcollision, nnew)
      call util_resize(self%lencounter, nnew)
      call util_resize(self%lmtiny, nnew)
      call util_resize(self%nplenc, nnew)
      call util_resize(self%ntpenc, nnew)
      call util_resize(self%levelg, nnew)
      call util_resize(self%levelm, nnew)
      call util_resize(self%isperi, nnew)
      call util_resize(self%peri, nnew)
      call util_resize(self%atp, nnew)
      call util_resize(self%kin, nnew)

      call util_resize_pl(self, nnew)

      return
   end subroutine symba_util_resize_pl


   subroutine symba_util_save_storage(system, snapshot, t)
      !! author: David A. Minton
      !!
      !! Checks the current size of the encounter storage against the required size and extends it by a factor of 2 more than requested if it is too small.
      !! Note: The reason to extend it by a factor of 2 is for performance. When there are many enounters per step, resizing every time you want to add an 
      !! encounter takes significant computational effort. Resizing by a factor of 2 is a tradeoff between performance (fewer resize calls) and memory managment
      !! Memory usage grows by a factor of 2 each time it fills up, but no more. 
      implicit none
      ! Arguments
      type(symba_nbody_system),  intent(inout) :: system   !! SyMBA nbody system object 
      class(encounter_snapshot), intent(in)    :: snapshot !! Encounter snapshot object
      real(DP),                  intent(in)    :: t        !! The time of the snapshot
      ! Internals
      type(encounter_storage(nframes=:)), allocatable :: tmp
      integer(I4B) :: i, nnew, nold, nbig

      ! Advance the snapshot frame counter
      system%encounter_history%iframe = system%encounter_history%iframe + 1

      ! Check to make sure the current encounter_history object is big enough. If not, grow it by a factor of 2
      nnew = system%encounter_history%iframe
      nold = system%encounter_history%nframes

      if (nnew > nold) then
         nbig = nold
         do while (nbig < nnew)
            nbig = nbig * 2
         end do
         allocate(encounter_storage(nbig) :: tmp) 
         tmp%tvals(1:nold) = system%encounter_history%tvals(1:nold)
         tmp%tvals(nold+1:nbig) = huge(1.0_DP)
         tmp%tslot(1:nold) = system%encounter_history%tslot(1:nold)
         tmp%tslot(nold+1:nbig) = 0
         tmp%iframe = system%encounter_history%iframe
         call move_alloc(system%encounter_history%nce, tmp%nce)
         call move_alloc(system%encounter_history%ncc, tmp%ncc)

         do i = 1, nold
            if (allocated(system%encounter_history%frame(i)%item)) call move_alloc(system%encounter_history%frame(i)%item, tmp%frame(i)%item)
         end do
         deallocate(system%encounter_history)
         call move_alloc(tmp,system%encounter_history)
         nnew = nbig
      end if

      ! Find out which time slot this belongs in by searching for an existing slot
      ! with the same value of time or the first available one
      do i = 1, nnew
         if (t <= system%encounter_history%tvals(i)) then
            system%encounter_history%tvals(i) = t
            system%encounter_history%tslot(nnew) = i
            system%encounter_history%frame(nnew) = snapshot
            exit
         end if
      end do

      return
   end subroutine symba_util_save_storage


   module subroutine symba_util_resize_tp(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a test particle object against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(symba_tp), intent(inout) :: self  !! SyMBA test particle object
      integer(I4B),    intent(in)    :: nnew  !! New size neded

      call util_resize(self%nplenc, nnew)
      call util_resize(self%levelg, nnew)
      call util_resize(self%levelm, nnew)

      call util_resize_tp(self, nnew)

      return
   end subroutine symba_util_resize_tp

   
   module subroutine symba_util_set_renc(self, scale)
      !! author: David A. Minton
      !!
      !! Sets the critical radius for encounter given an input recursion depth
      !!
      implicit none
      ! Arguments
      class(symba_pl), intent(inout) :: self !! SyMBA massive body object
      integer(I4B),    intent(in)    :: scale !! Current recursion depth
      ! Internals
      integer(I4B) :: i
      real(DP)     :: rshell_irec

      associate(pl => self, npl => self%nbody)
         rshell_irec = 1._DP
         do i = 1, scale
            rshell_irec = rshell_irec * RSHELL
         end do
         pl%renc(1:npl) = pl%rhill(1:npl) * RHSCALE * rshell_irec
      end associate

      return
   end subroutine symba_util_set_renc


   module subroutine symba_util_sort_pl(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a SyMBA massive body object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(symba_pl), intent(inout) :: self      !! SyMBA massive body object
      character(*),    intent(in)    :: sortby    !! Sorting attribute
      logical,         intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B)                        :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(pl => self, npl => self%nbody)
         select case(sortby)
         case("nplenc")
            call util_sort(direction * pl%nplenc(1:npl), ind)
         case("ntpenc")
            call util_sort(direction * pl%ntpenc(1:npl), ind)
         case("levelg")
            call util_sort(direction * pl%levelg(1:npl), ind)
         case("levelm")
            call util_sort(direction * pl%levelm(1:npl), ind)
         case("peri")
            call util_sort(direction * pl%peri(1:npl), ind)
         case("atp")
            call util_sort(direction * pl%atp(1:npl), ind)
         case("lcollision", "lencounter", "lmtiny", "nplm", "nplplm", "kin", "info")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default ! Look for components in the parent class
            call util_sort_pl(pl, sortby, ascending)
            return
         end select

         call pl%rearrange(ind)

      end associate
      return
   end subroutine symba_util_sort_pl


   module subroutine symba_util_sort_tp(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a SyMBA test particle object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(symba_tp), intent(inout) :: self      !! SyMBA test particle object
      character(*),    intent(in)    :: sortby    !! Sorting attribute
      logical,         intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B)                        :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(tp => self, ntp => self%nbody)
         select case(sortby)
         case("nplenc")
            call util_sort(direction * tp%nplenc(1:ntp), ind)
         case("levelg")
            call util_sort(direction * tp%levelg(1:ntp), ind)
         case("levelm")
            call util_sort(direction * tp%levelm(1:ntp), ind)
         case default ! Look for components in the parent class
            call util_sort_tp(tp, sortby, ascending)
            return
         end select

         call tp%rearrange(ind)
      end associate

      return
   end subroutine symba_util_sort_tp



   module subroutine symba_util_sort_rearrange_arr_kin(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of particle kinship type in-place from an index list.
      implicit none
      ! Arguments
      type(symba_kinship),  dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B),         dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                                    intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      type(symba_kinship),  dimension(:), allocatable                :: tmp !! Temporary copy of array used during rearrange operation
      integer(I4B) :: i,j

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, source=arr)
      tmp(1:n) = arr(ind(1:n))

      do i = 1, n
         do j = 1, tmp(i)%nchild
            tmp(i)%child(j) = ind(tmp(i)%child(j))
         end do
      end do

      call move_alloc(tmp, arr)
      return
   end subroutine symba_util_sort_rearrange_arr_kin


   module subroutine symba_util_sort_rearrange_pl(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange SyMBA massive body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(symba_pl),               intent(inout) :: self !! SyMBA massive body object
      integer(I4B),    dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

      associate(pl => self, npl => self%nbody)
         call util_sort_rearrange(pl%lcollision, ind, npl)
         call util_sort_rearrange(pl%lencounter, ind, npl)
         call util_sort_rearrange(pl%lmtiny,     ind, npl)
         call util_sort_rearrange(pl%nplenc,     ind, npl)
         call util_sort_rearrange(pl%ntpenc,     ind, npl)
         call util_sort_rearrange(pl%levelg,     ind, npl)
         call util_sort_rearrange(pl%levelm,     ind, npl)
         call util_sort_rearrange(pl%isperi,     ind, npl)
         call util_sort_rearrange(pl%peri,       ind, npl)
         call util_sort_rearrange(pl%atp,        ind, npl)
         call util_sort_rearrange(pl%kin,        ind, npl)

         call util_sort_rearrange_pl(pl,ind)
      end associate

      return
   end subroutine symba_util_sort_rearrange_pl


   module subroutine symba_util_sort_rearrange_tp(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange SyMBA test particle object in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(symba_tp),               intent(inout) :: self !! SyMBA test particle object
      integer(I4B),    dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

      associate(tp => self, ntp => self%nbody)
         call util_sort_rearrange(tp%nplenc, ind, ntp)
         call util_sort_rearrange(tp%levelg, ind, ntp)
         call util_sort_rearrange(tp%levelm, ind, ntp)

         call util_sort_rearrange_tp(tp,ind)
      end associate
      
      return
   end subroutine symba_util_sort_rearrange_tp


   module subroutine symba_util_spill_arr_kin(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of particle kinships
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      type(symba_kinship), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      type(symba_kinship), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,             dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
      logical,                                        intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: nspill, nkeep, nlist
      type(symba_kinship), dimension(:), allocatable :: tmp

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
      if (.not.allocated(discards)) then
         allocate(discards(nspill))
      else if (size(discards) /= nspill) then
         deallocate(discards)
         allocate(discards(nspill))
      end if

      discards(:) = pack(keeps(1:nlist), lspill_list(1:nlist))
      if (ldestructive) then
         if (nkeep > 0) then
            allocate(tmp(nkeep))
            tmp(:) = pack(keeps(1:nlist), .not. lspill_list(1:nlist))
            call move_alloc(tmp, keeps)
         else
            deallocate(keeps)
         end if
      end if

      return
   end subroutine symba_util_spill_arr_kin


   module subroutine symba_util_spill_pl(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) SyMBA massive body particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(symba_pl),       intent(inout) :: self        !! SyMBA massive body object
      class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Spill all the common components
      associate(keeps => self)
         select type(discards)
         class is (symba_pl)
            call util_spill(keeps%lcollision, discards%lcollision, lspill_list, ldestructive)
            call util_spill(keeps%lencounter, discards%lencounter, lspill_list, ldestructive)
            call util_spill(keeps%lmtiny, discards%lmtiny, lspill_list, ldestructive)
            call util_spill(keeps%nplenc, discards%nplenc, lspill_list, ldestructive)
            call util_spill(keeps%ntpenc, discards%ntpenc, lspill_list, ldestructive)
            call util_spill(keeps%levelg, discards%levelg, lspill_list, ldestructive)
            call util_spill(keeps%levelm, discards%levelm, lspill_list, ldestructive)
            call util_spill(keeps%isperi, discards%isperi, lspill_list, ldestructive)
            call util_spill(keeps%peri, discards%peri, lspill_list, ldestructive)
            call util_spill(keeps%atp, discards%atp, lspill_list, ldestructive)
            call util_spill(keeps%kin, discards%kin, lspill_list, ldestructive)

            call util_spill_pl(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class symba_pl or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate
     
      return
   end subroutine symba_util_spill_pl


   module subroutine symba_util_spill_encounter_list(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) SyMBA encounter structure from active list to discard list
      !! Note: Because the symba_plplenc currently does not contain any additional variable components, this method can recieve it as an input as well.
      implicit none
      ! Arguments
      class(symba_encounter), intent(inout) :: self         !! SyMBA pl-tp encounter list 
      class(encounter_list),  intent(inout) :: discards     !! Discarded object 
      logical, dimension(:),  intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
  
      associate(keeps => self)
         select type(discards)
         class is (symba_encounter)
            call util_spill(keeps%level, discards%level, lspill_list, ldestructive)
            call encounter_util_spill_list(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class symba_encounter or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate
   
      return
   end subroutine symba_util_spill_encounter_list


   module subroutine symba_util_spill_tp(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) SyMBA test particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(symba_tp),       intent(inout) :: self         !! SyMBA test particle object
      class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Spill all the common components
      associate(keeps => self)
         select type(discards)
         class is (symba_tp)
            call util_spill(keeps%nplenc, discards%nplenc, lspill_list, ldestructive)
            call util_spill(keeps%levelg, discards%levelg, lspill_list, ldestructive)
            call util_spill(keeps%levelm, discards%levelm, lspill_list, ldestructive)

            call util_spill_tp(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class symba_tp or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate
     
      return
   end subroutine symba_util_spill_tp


   module subroutine symba_util_take_collision_snapshot(self, param, t, stage)
      !! author: David A. Minton
      !!
      !! Takes a minimal snapshot of the state of the system during an encounter so that the trajectories
      !! can be played back through the encounter
      implicit none
      ! Internals
      class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t     !! current time
      character(*),               intent(in)    :: stage !! Either before or after
      ! Arguments
      class(fraggle_collision_snapshot), allocatable:: snapshot
      integer(I4B) :: i,j

      select case(stage)
      case("before")
         ! Saves the states of the bodies involved in the collision before the collision is resolved
         associate (idx => self%colliders%idx, ncoll => self%colliders%ncoll)
            allocate(symba_pl :: self%colliders%pl)
            select type(pl => self%colliders%pl)
            class is (symba_pl)
               call pl%setup(ncoll, param)
               pl%id(:) = self%pl%id(idx(:))
               pl%Gmass(:) = self%pl%Gmass(idx(:))
               pl%radius(:) = self%pl%radius(idx(:))
               pl%rot(:,:) = self%pl%rot(:,idx(:))
               pl%Ip(:,:) = self%pl%Ip(:,idx(:))
               pl%rh(:,:) = self%pl%rh(:,idx(:))
               pl%vh(:,:) = self%pl%vh(:,idx(:))
               pl%info(:) = self%pl%info(idx(:))
            end select
         end associate
      case("after")
         allocate(fraggle_collision_snapshot :: snapshot)
         allocate(snapshot%colliders, source=self%colliders) 
         allocate(snapshot%fragments, source=self%fragments)
         !call symba_util_save_storage(self,snapshot,t)
      end select

      return
   end subroutine symba_util_take_collision_snapshot

   module subroutine symba_util_take_encounter_snapshot(self, param, t)
      !! author: David A. Minton
      !!
      !! Takes a minimal snapshot of the state of the system during an encounter so that the trajectories
      !! can be played back through the encounter
      implicit none
      ! Internals
      class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t     !! current time
      ! Arguments
      class(encounter_snapshot), allocatable :: snapshot
      integer(I4B) :: i, npl_snap, ntp_snap

      associate(npl => self%pl%nbody,  ntp => self%tp%nbody)

         if (self%plplenc_list%lcollision) then
            allocate(fraggle_collision_snapshot :: snapshot)
         else
            allocate(encounter_snapshot :: snapshot)
         end if
         snapshot%t = t
         snapshot%iloop = param%iloop

         if (npl + ntp == 0) return
         npl_snap = npl
         ntp_snap = ntp

         select type (pl => self%pl)
         class is (symba_pl)
            select type (tp => self%tp)
            class is (symba_tp)
               allocate(symba_pl :: snapshot%pl)
               allocate(symba_tp :: snapshot%tp)

               select type(pl_snap => snapshot%pl)
               class is (symba_pl)
                  select type(tp_snap => snapshot%tp)
                  class is (symba_tp)

                     if (npl > 0) then
                        pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE .and. pl%levelg(1:npl) == self%irec
                        npl_snap = count(pl%lmask(1:npl))
                     end if
                     if (ntp > 0) then
                        tp%lmask(1:ntp) = tp%status(1:ntp) /= INACTIVE .and. tp%levelg(1:ntp) == self%irec
                        ntp_snap = count(tp%lmask(1:ntp))
                     end if
                     pl_snap%nbody = npl_snap

                     ! Take snapshot of the currently encountering massive bodies
                     if (npl_snap > 0) then
                        allocate(pl_snap%id(npl_snap))
                        allocate(pl_snap%info(npl_snap))
                        allocate(pl_snap%Gmass(npl_snap))

                        allocate(pl_snap%levelg(npl_snap))
                        pl_snap%levelg(:) = pack(pl%levelg(1:npl), pl%lmask(1:npl))
                        pl_snap%id(:) = pack(pl%id(1:npl), pl%lmask(1:npl))
                        pl_snap%info(:) = pack(pl%info(1:npl), pl%lmask(1:npl))
                        pl_snap%Gmass(:) = pack(pl%Gmass(1:npl), pl%lmask(1:npl))
                        allocate(pl_snap%rh(NDIM,npl_snap))
                        allocate(pl_snap%vh(NDIM,npl_snap))
                        do i = 1, NDIM
                           pl_snap%rh(i,:) = pack(pl%rh(i,1:npl), pl%lmask(1:npl))
                           pl_snap%vh(i,:) = pack(pl%vb(i,1:npl), pl%lmask(1:npl))
                        end do
                        if (param%lclose) then
                           allocate(pl_snap%radius(npl_snap))
                           pl_snap%radius(:) = pack(pl%radius(1:npl), pl%lmask(1:npl))
                        end if

                        if (param%lrotation) then
                           allocate(pl_snap%Ip(NDIM,npl_snap))
                           allocate(pl_snap%rot(NDIM,npl_snap))
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
                        allocate(tp_snap%id(ntp_snap))
                        allocate(tp_snap%info(ntp_snap))
                        tp_snap%id(:) = pack(tp%id(1:ntp), tp%lmask(1:ntp))
                        tp_snap%info(:) = pack(tp%info(1:ntp), tp%lmask(1:ntp))
                        allocate(tp_snap%rh(NDIM,ntp_snap))
                        allocate(tp_snap%vh(NDIM,ntp_snap))
                        do i = 1, NDIM
                           tp_snap%rh(i,:) = pack(tp%rh(i,1:ntp), tp%lmask(1:ntp))
                           tp_snap%vh(i,:) = pack(tp%vh(i,1:ntp), tp%lmask(1:ntp))
                        end do
                     end if
                  end select
               end select

               ! Save the snapshot
               call symba_util_save_storage(self,snapshot,t)
            end select
         end select
      end associate

      return
   end subroutine symba_util_take_encounter_snapshot

end submodule s_symba_util
