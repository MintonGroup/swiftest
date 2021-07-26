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
      integer(I4B)                              :: nenc_old
      integer(I8B)                              :: k
      real(DP),     dimension(NDIM)             :: xr, vr
      logical,      dimension(:),   allocatable :: lencounter, loc_lvdotr
   
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

         lany_encounter = any(lencounter(:))
         if (lany_encounter) then 
            associate(plplenc_list => system%plplenc_list, nenc => system%plplenc_list%nenc)
               nenc_old = nenc
               call plplenc_list%resize(nenc_old + count(lencounter(:)))
               plplenc_list%status(nenc_old+1:nenc) = ACTIVE
               plplenc_list%level(nenc_old+1:nenc) = irec
               plplenc_list%lvdotr(nenc_old+1:nenc) = pack(loc_lvdotr(:), lencounter(:))
               plplenc_list%index1(nenc_old+1:nenc) = pack(pl%k_plpl(1,:), lencounter(:))
               plplenc_list%index2(nenc_old+1:nenc) = pack(pl%k_plpl(2,:), lencounter(:))
               pl%lencounter(plplenc_list%index1(nenc_old+1:nenc)) = .true.
               pl%lencounter(plplenc_list%index2(nenc_old+1:nenc)) = .true.
            end associate
         end if
      end associate
      return
   end function symba_encounter_check_pl

   module function symba_encounter_check_pltpenc(self, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between test particles and massive bodies in the pltpenc list.
      !! Note: This method works for the polymorphic symba_pltpenc and symba_plplenc types.
      !!
      !! Adapted from portions of David E. Kaufmann's Swifter routine: symba_step_recur.f90
      implicit none
      ! Arguments
      class(symba_pltpenc),      intent(inout) :: self       !! SyMBA pl-pl encounter list object
      class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt         !! step size
      integer(I4B),              intent(in)    :: irec       !! Current recursion level 
      logical                                  :: lany_encounter !! Returns true if there is at least one close encounter  
      ! Internals
      integer(I4B)              :: i
      real(DP), dimension(NDIM) :: xr, vr
      logical                   :: lencounter, isplpl
      real(DP)                  :: rlim2, rji2

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
            do i = 1, self%nenc
               associate(index_i => self%index1(i), index_j => self%index2(i))
                  if (isplpl) then
                     xr(:) = pl%xh(:,index_j) - pl%xh(:,index_i)
                     vr(:) = pl%vb(:,index_j) - pl%vb(:,index_i)
                     call symba_encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%rhill(index_i), pl%rhill(index_j), dt, irec, lencounter, self%lvdotr(i))
                  else
                     xr(:) = tp%xh(:,index_j) - pl%xh(:,index_i)
                     vr(:) = tp%vb(:,index_j) - pl%vb(:,index_i)
                     call symba_encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%rhill(index_i), 0.0_DP, dt, irec, lencounter, self%lvdotr(i))
                  end if
                  if (lencounter) then
                     if (isplpl) then
                        rlim2 = (pl%radius(index_i) + pl%radius(index_j))**2
                     else
                        rlim2 = (pl%radius(index_i))**2
                     end if
                     rji2 = dot_product(xr(:), xr(:))! Check to see if these are physically overlapping bodies first, which we should ignore
                     if (rji2 > rlim2) then
                        lany_encounter = .true.
                        pl%levelg(index_i) = irec
                        pl%levelm(index_i) = MAX(irec, pl%levelm(index_i))
                        if (isplpl) then
                           pl%levelg(index_j) = irec
                           pl%levelm(index_j) = MAX(irec, pl%levelm(index_j))
                        else
                           tp%levelg(index_j) = irec
                           tp%levelm(index_j) = MAX(irec, tp%levelm(index_j))
                        end if
                        self%level(i) = irec
                     end if
                  end if   
               end associate
            end do
      end select
   end select
   end function symba_encounter_check_pltpenc

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
      integer(I4B)                              :: i, j, nenc_old
      real(DP),     dimension(NDIM)             :: xr, vr
      logical,      dimension(:,:), allocatable :: lencounter, loc_lvdotr
   
      associate(tp => self, ntp => self%nbody, pl => system%pl, npl => system%pl%nbody)
         allocate(lencounter(npl, ntp), loc_lvdotr(npl, ntp))
         lencounter(:,:) = .false.
   
         do j = 1, ntp
            do i = 1, npl
               xr(:) = tp%xh(:, j) - pl%xh(:, i)
               vr(:) = tp%vh(:, j) - pl%vh(:, i)
               call symba_encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%rhill(i), 0.0_DP, dt, irec, lencounter(i,j), loc_lvdotr(i,j))
            end do
         end do

         lany_encounter = any(lencounter(:,:))
         if (lany_encounter) then 
            associate(pltpenc_list => system%pltpenc_list, nenc => system%pltpenc_list%nenc)
               nenc_old = nenc
               call pltpenc_list%resize(nenc_old + count(lencounter(:,:)))
               pltpenc_list%status(nenc_old+1:nenc) = ACTIVE
               pltpenc_list%level(nenc_old+1:nenc) = irec
               pltpenc_list%lvdotr(nenc_old+1:nenc) = pack(loc_lvdotr(:,:), lencounter(:,:))
               pltpenc_list%index1(nenc_old+1:nenc) = pack(spread([(i, i = 1, npl)], dim=2, ncopies=ntp), lencounter(:,:)) 
               pltpenc_list%index2(nenc_old+1:nenc) = pack(spread([(j, j = 1, ntp)], dim=1, ncopies=npl), lencounter(:,:))
               select type(pl)
               class is (symba_pl)
                  pl%lencounter(pltpenc_list%index1(nenc_old+1:nenc)) = .true.
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
      integer(I4B) :: iflag
      real(DP)     :: r2, v2, rcrit, r2crit, vdotr

      lencounter = .false.
      rcrit = (rhill1 + rhill2)*RHSCALE*(RSHELL**(irec))
      r2crit = rcrit**2
      r2 = xr**2 + yr**2 + zr**2
      v2 = vxr**2 + vyr**2 + vzr**2
      vdotr = xr * vxr + yr * vyr + zr * vzr
      iflag = rmvs_chk_ind(r2, v2, vdotr, dt, r2crit)
      if (iflag /= 0) lencounter = .true.
      lvdotr = (vdotr < 0.0_DP)

      return
   end subroutine symba_encounter_check_one


end submodule s_symba_encounter_check