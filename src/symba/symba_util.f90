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


   module subroutine symba_util_append_encounter(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one encounter list (pl-pl or pl-tp) body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(symba_encounter),    intent(inout) :: self         !! SyMBA encounter list object
      class(swiftest_encounter), intent(in)    :: source       !! Source object to append
      logical, dimension(:),     intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      associate(nold => self%nenc, nsrc => source%nenc)
         select type(source)
         class is (symba_encounter)
            call util_append(self%level, source%level, nold, nsrc, lsource_mask)
         end select
         call util_append_encounter(self, source, lsource_mask) 
      end associate

      return
   end subroutine symba_util_append_encounter


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


   module subroutine symba_util_copy_encounter(self, source)
      !! author: David A. Minton
      !!
      !! Copies elements from the source encounter list into self.
      implicit none
      ! Arguments
      class(symba_encounter),    intent(inout) :: self   !! Encounter list 
      class(swiftest_encounter), intent(in)    :: source !! Source object to copy into
  
      select type(source)
      class is (symba_encounter)
         associate(n => source%nenc)
            self%level(1:n) = source%level(1:n) 
         end associate
      end select

      call util_copy_encounter(self, source)
   
      return
   end subroutine symba_util_copy_encounter


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


   module subroutine symba_util_index_eucl_plpl(self, param)
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
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      ! Internals
      integer(I8B) :: i, j, counter, npl, nplm, nplpl, nplplm

      associate(pl => self)
         npl = int(self%nbody, kind=I8B)
         call pl%sort("mass", ascending=.false.)

         nplm = count(.not. pl%lmtiny(1:npl))
         pl%nplm = int(nplm, kind=I4B)

         pl%nplpl = (npl * (npl - 1) / 2) ! number of entries in a strict lower triangle, npl x npl, minus first column
         pl%nplplm = nplm * npl - nplm * (nplm + 1) / 2 ! number of entries in a strict lower triangle, npl x npl, minus first column including only mutually interacting bodies
         if (allocated(self%k_plpl)) deallocate(self%k_plpl) ! Reset the index array if it's been set previously
         allocate(self%k_plpl(2, pl%nplpl))
         do i = 1, npl
            counter = (i - 1_I8B) * npl - i * (i - 1_I8B) / 2_I8B + 1_I8B
            do j = i + 1_I8B, npl
               self%k_plpl(1, counter) = i
               self%k_plpl(2, counter) = j
               counter = counter + 1_I8B
            end do
         end do
      end associate

      return
   end subroutine symba_util_index_eucl_plpl


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
                     vdotr = dot_product(pl%xh(:,i), pl%vh(:,i))
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
                     vdotr = dot_product(pl%xh(:,i), pl%vh(:,i))
                     if (pl%isperi(i) == -1) then
                        if (vdotr >= 0.0_DP) then
                           pl%isperi(i) = 0
                           CALL orbel_xv2aeq(pl%mu(i), pl%xh(:,i), pl%vh(:,i), pl%atp(i), e, pl%peri(i))
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
                           CALL orbel_xv2aeq(system%Gmtot, pl%xb(:,i), pl%vb(:,i), pl%atp(i), e, pl%peri(i))
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
      integer(I4B) :: i, j, k, npl, nadd, nencmin, idnew1, idnew2, idold1, idold2
      logical, dimension(:), allocatable :: lmask, ldump_mask
      class(symba_plplenc), allocatable :: plplenc_old
      logical :: lencounter
      integer(I4B), dimension(:), allocatable :: levelg_orig_pl, levelm_orig_pl, levelg_orig_tp, levelm_orig_tp, nplenc_orig_pl, nplenc_orig_tp, ntpenc_orig_pl

      associate(pl => self, pl_adds => system%pl_adds)

         npl = pl%nbody
         nadd = pl_adds%nbody
         if (npl == 0) return
         ! Deallocate any temporary variables
         if (allocated(pl%xbeg)) deallocate(pl%xbeg)
         if (allocated(pl%xend)) deallocate(pl%xend)

         ! Remove the discards and destroy the list, as the system already tracks pl_discards elsewhere
         allocate(lmask(npl))
         lmask(1:npl) = pl%ldiscard(1:npl) .or. pl%status(1:npl) == INACTIVE
         allocate(tmp, mold=self)
         call pl%spill(tmp, lspill_list=lmask, ldestructive=.true.)
         npl = pl%nbody
         call tmp%setup(0,param)
         deallocate(tmp)
         deallocate(lmask)

         ! Store the original plplenc list so we don't remove any of the original encounters
         allocate(plplenc_old, source=system%plplenc_list)
         call plplenc_old%copy(system%plplenc_list)

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
         where(pl%status(1:npl) /= INACTIVE) 
            pl%status(1:npl) = ACTIVE
            pl%info(1:npl)%status = "ACTIVE"
            pl%ldiscard(1:npl) = .false.
            pl%lcollision(1:npl) = .false.
            pl%lmask(1:npl) = .true.
         elsewhere
            pl%ldiscard(1:npl) = .true.
            pl%lmask(1:npl) = .false.
         end where

         select type(param)
         class is (symba_parameters)
            pl%lmtiny(1:npl) = pl%Gmass(1:npl) < param%GMTINY
            where(pl%lmtiny(1:npl))
               pl%info(1:npl)%particle_type = PL_TINY_TYPE_NAME 
            elsewhere
               pl%info(1:npl)%particle_type = PL_TYPE_NAME 
            end where
         end select

         call pl%dump_particle_info(param, idx=pack([(i, i=1, npl)], ldump_mask))
         deallocate(ldump_mask)

         ! Reindex the new list of bodies 
         call pl%index(param)

         ! Reset the kinship trackers
         pl%kin(1:npl)%nchild = 0
         pl%kin(1:npl)%parent = [(i, i=1, npl)]
          do i = 1, npl
             if (allocated(pl%kin(i)%child)) deallocate(pl%kin(i)%child)
          end do

         ! Re-build the zero-level encounter list, being sure to save the original level information for all bodies
         allocate(levelg_orig_pl, source=pl%levelg)
         allocate(levelm_orig_pl, source=pl%levelm)
         allocate(nplenc_orig_pl, source=pl%nplenc)
         lencounter = pl%encounter_check(system, param%dt, 0) 
         if (system%tp%nbody > 0) then
            select type(tp => system%tp)
            class is (symba_tp)
               allocate(ntpenc_orig_pl, source=pl%ntpenc)
               allocate(levelg_orig_tp, source=tp%levelg)
               allocate(levelm_orig_tp, source=tp%levelm)
               allocate(nplenc_orig_tp, source=tp%nplenc)
               lencounter = tp%encounter_check(system, param%dt, 0)
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
         nencmin = min(system%plplenc_list%nenc, plplenc_old%nenc) 
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
               system%plplenc_list%t(k) = plplenc_old%t(k)
               system%plplenc_list%level(k) = plplenc_old%level(k)
            else if (((idnew1 == idold2) .and. (idnew2 == idold1))) then
               ! This is an encounter we already know about, but with the order reversed, so save the old information
               system%plplenc_list%lvdotr(k) = plplenc_old%lvdotr(k) 
               system%plplenc_list%status(k) = plplenc_old%status(k) 
               system%plplenc_list%x1(:,k) = plplenc_old%x2(:,k)
               system%plplenc_list%x2(:,k) = plplenc_old%x1(:,k)
               system%plplenc_list%v1(:,k) = plplenc_old%v2(:,k)
               system%plplenc_list%v2(:,k) = plplenc_old%v1(:,k)
               system%plplenc_list%t(k) = plplenc_old%t(k)
               system%plplenc_list%level(k) = plplenc_old%level(k)
            end if
         end do
      end associate

      return
   end subroutine symba_util_rearray_pl


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

      if (.not. allocated(arr) .or. nnew < 0) return

      nold = size(arr)
      if (nnew == nold) return

      if (nnew == 0) then
         deallocate(arr)
         return
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
         allocate(ind(npl))
         select case(sortby)
         case("nplenc")
            call util_sort(direction * pl%nplenc(1:npl), ind(1:npl))
         case("ntpenc")
            call util_sort(direction * pl%ntpenc(1:npl), ind(1:npl))
         case("levelg")
            call util_sort(direction * pl%levelg(1:npl), ind(1:npl))
         case("levelm")
            call util_sort(direction * pl%levelm(1:npl), ind(1:npl))
         case("peri")
            call util_sort(direction * pl%peri(1:npl), ind(1:npl))
         case("atp")
            call util_sort(direction * pl%atp(1:npl), ind(1:npl))
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
         allocate(ind(ntp))
         select case(sortby)
         case("nplenc")
            call util_sort(direction * tp%nplenc(1:ntp), ind(1:ntp))
         case("levelg")
            call util_sort(direction * tp%levelg(1:ntp), ind(1:ntp))
         case("levelm")
            call util_sort(direction * tp%levelm(1:ntp), ind(1:ntp))
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
      ! Internals
      integer(I4B) :: i, j

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
            keeps(:) = pack(keeps(1:nlist), .not. lspill_list(1:nlist))
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


   module subroutine symba_util_spill_encounter(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) SyMBA encounter structure from active list to discard list
      !! Note: Because the symba_plplenc currently does not contain any additional variable components, this method can recieve it as an input as well.
      implicit none
      ! Arguments
      class(symba_encounter),      intent(inout) :: self         !! SyMBA pl-tp encounter list 
      class(swiftest_encounter), intent(inout) :: discards     !! Discarded object 
      logical, dimension(:),     intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                   intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
  
      associate(keeps => self)
         select type(discards)
         class is (symba_encounter)
            call util_spill(keeps%level, discards%level, lspill_list, ldestructive)
            call util_spill_encounter(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class symba_encounter or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate
   
      return
   end subroutine symba_util_spill_encounter


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

end submodule s_symba_util
