submodule (symba_classes) s_symba_encounter_check
   use swiftest
contains

   module function symba_encounter_check_pl(self, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between massive bodies.
      !!
      implicit none
      ! Arguments
      class(symba_pl),           intent(inout)  :: self           !! SyMBA test particle object  
      class(symba_nbody_system), intent(inout)  :: system         !! SyMBA nbody system object
      real(DP),                  intent(in)     :: dt             !! step size
      integer(I4B),              intent(in)     :: irec           !! Current recursion level
      ! Result
      logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      ! Internals
      integer(I8B)                              :: k
      integer(I4B)                              :: nenc
      real(DP),     dimension(NDIM)             :: xr, vr
      logical,      dimension(:),   allocatable :: lencounter, loc_lvdotr
  
      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody, nplpl => self%nplpl)
         allocate(lencounter(nplpl), loc_lvdotr(nplpl))
         lencounter(:) = .false.
   
         do k = 1, nplpl
            associate(i => pl%k_plpl(1, k), j => pl%k_plpl(2, k))
               xr(:) = pl%xh(:, j) - pl%xh(:, i)
               vr(:) = pl%vh(:, j) - pl%vh(:, i)
               call symba_encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%rhill(i), pl%rhill(j), dt, irec, lencounter(k), loc_lvdotr(k))
            end associate
         end do

         nenc = count(lencounter(:))
         lany_encounter = nenc > 0
         if (lany_encounter) then 
            associate(plplenc_list => system%plplenc_list)
               call plplenc_list%resize(nenc)
               plplenc_list%lvdotr(1:nenc) = pack(loc_lvdotr(:), lencounter(:))
               plplenc_list%index1(1:nenc) = pack(pl%k_plpl(1,:), lencounter(:))
               plplenc_list%index2(1:nenc) = pack(pl%k_plpl(2,:), lencounter(:))
               plplenc_list%id1(1:nenc) = pl%id(plplenc_list%index1(1:nenc))
               plplenc_list%id2(1:nenc) = pl%id(plplenc_list%index2(1:nenc))
               do k = 1, nenc
                  plplenc_list%status(k) = ACTIVE
                  plplenc_list%level(k) = irec
                  pl%lencounter(plplenc_list%index1(k)) = .true.
                  pl%lencounter(plplenc_list%index2(k)) = .true.
               end do
            end associate
         end if
      end associate
      return
   end function symba_encounter_check_pl


   module function symba_encounter_check(self, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between test particles and massive bodies in the pltpenc list.
      !! Note: This method works for the polymorphic symba_pltpenc and symba_plplenc types.
      !!
      !! Adapted from portions of David E. Kaufmann's Swifter routine: symba_step_recur.f90
      implicit none
      ! Arguments
      class(symba_encounter),      intent(inout) :: self       !! SyMBA pl-pl encounter list object
      class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt         !! step size
      integer(I4B),              intent(in)    :: irec       !! Current recursion level 
      logical                                  :: lany_encounter !! Returns true if there is at least one close encounter  
      ! Internals
      integer(I4B)              :: k
      real(DP), dimension(NDIM) :: xr, vr
      logical                   :: lencounter, isplpl
      real(DP)                  :: rlim2, rji2
      logical, dimension(:), allocatable :: lencmask

      lany_encounter = .false.
      if (self%nenc == 0) return

      select type(self)
      class is (symba_plplenc)
         isplpl = .true.
      class is (symba_pltpenc)
         isplpl = .false.
      end select

      select type(pl => system%pl)
      class is (symba_pl)
         select type(tp => system%tp)
         class is (symba_tp)
            allocate(lencmask(self%nenc))
            lencmask(:) = (self%status(1:self%nenc) == ACTIVE) .and. (self%level(1:self%nenc) == irec - 1)
            if (.not.any(lencmask(:))) return
            associate(ind1 => self%index1, ind2 => self%index2) 
               do concurrent(k = 1:self%nenc, lencmask(k))
                  if (isplpl) then
                     xr(:) = pl%xh(:,ind2(k)) - pl%xh(:,ind1(k))
                     vr(:) = pl%vb(:,ind2(k)) - pl%vb(:,ind1(k))
                     call symba_encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%rhill(ind1(k)), pl%rhill(ind2(k)), dt, irec, lencounter, self%lvdotr(k))
                  else
                     xr(:) = tp%xh(:,ind2(k)) - pl%xh(:,ind1(k))
                     vr(:) = tp%vb(:,ind2(k)) - pl%vb(:,ind1(k))
                     call symba_encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%rhill(ind1(k)), 0.0_DP, dt, irec, lencounter, self%lvdotr(k))
                  end if
                  if (lencounter) then
                     if (isplpl) then
                        rlim2 = (pl%radius(ind1(k)) + pl%radius(ind2(k)))**2
                     else
                        rlim2 = (pl%radius(ind1(k)))**2
                     end if
                     rji2 = dot_product(xr(:), xr(:))! Check to see if these are physically overlapping bodies first, which we should ignore
                     if (rji2 > rlim2) then
                        lany_encounter = .true.
                        pl%levelg(ind1(k)) = irec
                        pl%levelm(ind1(k)) = MAX(irec, pl%levelm(ind1(k)))
                        if (isplpl) then
                           pl%levelg(ind2(k)) = irec
                           pl%levelm(ind2(k)) = MAX(irec, pl%levelm(ind2(k)))
                        else
                           tp%levelg(ind2(k)) = irec
                           tp%levelm(ind2(k)) = MAX(irec, tp%levelm(ind2(k)))
                        end if
                        self%level(k) = irec
                     end if
                  end if   
               end do
            end associate
         end select
      end select

      return
   end function symba_encounter_check


   module function symba_encounter_check_tp(self, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between test particles and massive bodies.
      !!
      implicit none
      ! Arguments
      class(symba_tp),           intent(inout) :: self       !! SyMBA test particle object  
      class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt         !! step size
      integer(I4B),              intent(in)    :: irec           !! Current recursion level
      ! Result
      logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      ! Internals
      real(DP)                                  :: r2crit, vdotr, r2, v2, tmin, r2min, term2
      integer(I4B)                              :: i, j, k,nenc
      real(DP),     dimension(NDIM)             :: xr, vr
      logical,      dimension(:,:), allocatable :: lencounter, loc_lvdotr
  
      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody, pl => system%pl, npl => system%pl%nbody)
         allocate(lencounter(ntp, npl), loc_lvdotr(ntp, npl))
         lencounter(:,:) = .false.
   
         do j = 1, npl
            do i = 1, ntp
               xr(:) = tp%xh(:, i) - pl%xh(:, j)
               vr(:) = tp%vh(:, i) - pl%vh(:, j)
               call symba_encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%rhill(j), 0.0_DP, dt, irec, lencounter(i,j), loc_lvdotr(i,j))
            end do
         end do

         nenc = count(lencounter(:,:))
         lany_encounter = nenc > 0
         if (lany_encounter) then 
            associate(pltpenc_list => system%pltpenc_list)
               call pltpenc_list%resize(nenc)
               pltpenc_list%status(1:nenc) = ACTIVE
               pltpenc_list%level(1:nenc) = irec
               pltpenc_list%lvdotr(1:nenc) = pack(loc_lvdotr(:,:), lencounter(:,:))
               pltpenc_list%index1(1:nenc) = pack(spread([(i, i = 1, npl)], dim=1, ncopies=ntp), lencounter(:,:)) 
               pltpenc_list%index2(1:nenc) = pack(spread([(i, i = 1, ntp)], dim=2, ncopies=npl), lencounter(:,:))
               pltpenc_list%id1(1:nenc) = pl%id(pltpenc_list%index1(1:nenc))
               pltpenc_list%id2(1:nenc) = tp%id(pltpenc_list%index2(1:nenc))
               select type(pl)
               class is (symba_pl)
                  pl%lencounter(:) = .false.
                  do k = 1, nenc
                     pl%lencounter(pltpenc_list%index1(k)) = .true.
                  end do
               end select
            end associate
         end if
      end associate

      return
   end function symba_encounter_check_tp


   module pure elemental subroutine symba_encounter_check_one(xr, yr, zr, vxr, vyr, vzr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for an encounter.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_chk.f90
      !! Adapted from Hal Levison's Swift routine symba5_chk.f
      implicit none
      ! Arguments
      real(DP),     intent(in)  :: xr, yr, zr, vxr, vyr, vzr
      real(DP),     intent(in)  :: rhill1, rhill2, dt
      integer(I4B), intent(in)  :: irec
      logical,      intent(out) :: lencounter, lvdotr
      ! Internals
      real(DP)     :: r2, v2, rcrit, r2crit, vdotr

      rcrit = (rhill1 + rhill2)*RHSCALE*(RSHELL**(irec))
      r2crit = rcrit**2
      r2 = xr**2 + yr**2 + zr**2
      v2 = vxr**2 + vyr**2 + vzr**2
      vdotr = xr * vxr + yr * vyr + zr * vzr
      lencounter = rmvs_chk_ind(r2, v2, vdotr, dt, r2crit)
      lvdotr = (vdotr < 0.0_DP)

      return
   end subroutine symba_encounter_check_one

end submodule s_symba_encounter_check