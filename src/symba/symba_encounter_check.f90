submodule (symba_classes) s_symba_encounter_check
   use swiftest
contains

   subroutine symba_encounter_check_all_flat(nplplm, k_plpl, x, v, rhill,  dt, irec, lencounter, loc_lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. Split off from the main subroutine for performance.
      !! This is the flat (single loop) version.
      implicit none
      integer(I8B), intent(in) :: nplplm
      integer(I4B), dimension(:,:), intent(in) :: k_plpl
      real(DP), dimension(:,:), intent(in) :: x, v
      real(DP), dimension(:), intent(in) :: rhill
      real(DP), intent(in) :: dt
      integer(I4B), intent(in) :: irec
      logical, dimension(:), intent(out) :: lencounter, loc_lvdotr
      ! Internals
      integer(I8B) :: k
      integer(I4B) :: i, j
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, rhill1, rhill2
      
      !$omp parallel do simd default(private) schedule(static)&
      !$omp shared(nplplm, k_plpl, x, v, rhill,  dt, irec, lencounter, loc_lvdotr) &
      !$omp lastprivate(xr, yr, zr, vxr, vyr, vzr, rhill1, rhill2)
      do k = 1_I8B, nplplm
         i = k_plpl(1, k)
         j = k_plpl(2, k)
         xr = x(1, j) - x(1, i)
         yr = x(2, j) - x(2, i)
         zr = x(3, j) - x(3, i)
         vxr = v(1, j) - v(1, i)
         vyr = v(2, j) - v(2, i)
         vzr = v(3, j) - v(3, i)
         rhill1 = rhill(i)
         rhill2 = rhill(j)
         call symba_encounter_check_one(xr, yr, zr, vxr, vyr, vzr, rhill1, rhill2, dt, irec, lencounter(k), loc_lvdotr(k))
      end do
      !$omp end parallel do simd

      return
   end subroutine symba_encounter_check_all_flat


   subroutine symba_encounter_check_all_triangular(npl, nplm, x, v, rhill,  dt, irec, loc_lvdotr, k_plplenc, nenc)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. Split off from the main subroutine for performance
      !! This is the upper triangular (double loop) version.
      implicit none
      integer(I4B), intent(in) :: npl, nplm
      real(DP), dimension(:,:), intent(in) :: x, v
      real(DP), dimension(:), intent(in) :: rhill
      real(DP), intent(in) :: dt
      integer(I4B), intent(in) :: irec
      logical, dimension(:), allocatable, intent(out) :: loc_lvdotr
      integer(I4B), dimension(:,:), allocatable, intent(out) :: k_plplenc
      integer(I4B), intent(out) :: nenc
      ! Internals
      integer(I4B) :: i, j, nenci, j0, j1
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, rhill1, rhill2
      logical, dimension(npl) :: lencounteri, loc_lvdotri
      integer(I4B), dimension(npl) :: ind_arr
      type lenctype
         logical, dimension(:), allocatable :: loc_lvdotr
         integer(I4B), dimension(:), allocatable :: index2
         integer(I4B) :: nenc
      end type
      type(lenctype), dimension(nplm) :: lenc
     
      ind_arr(:) = [(i, i = 1, npl)]
      !$omp parallel do simd default(private) schedule(static)&
      !$omp shared(npl, nplm, x, v, rhill,  dt, irec, lenc) &
      !$omp lastprivate(xr, yr, zr, vxr, vyr, vzr, rhill1, rhill2, lencounteri, loc_lvdotri)
      do i = 1, nplm
         do concurrent(j = i+1:npl)
            xr = x(1, j) - x(1, i)
            yr = x(2, j) - x(2, i)
            zr = x(3, j) - x(3, i)
            vxr = v(1, j) - v(1, i)
            vyr = v(2, j) - v(2, i)
            vzr = v(3, j) - v(3, i)
            rhill1 = rhill(i)
            rhill2 = rhill(j)
            call symba_encounter_check_one(xr, yr, zr, vxr, vyr, vzr, rhill1, rhill2, dt, irec, lencounteri(j), loc_lvdotri(j))
         end do
         lenc(i)%nenc = count(lencounteri(i+1:npl))
         if (lenc(i)%nenc > 0) then
            allocate(lenc(i)%loc_lvdotr(nenci), lenc(i)%index2(nenci))
            lenc(i)%loc_lvdotr(:) = pack(loc_lvdotri(i+1:npl), lencounteri(i+1:npl)) 
            lenc(i)%index2(:) = pack(ind_arr(i+1:npl), lencounteri(i+1:npl)) 
         end if
      end do
      !$omp end parallel do simd

      associate(nenc_arr => lenc(:)%nenc)
         nenc = sum(nenc_arr(1:nplm))
      end associate
      if (nenc > 0) then
         allocate(loc_lvdotr(nenc))
         allocate(k_plplenc(2,nenc))
         j0 = 1
         do i = 1, nplm
            if (lenc(i)%nenc > 0) then
               j1 = j0 + lenc(i)%nenc - 1
               loc_lvdotr(j0:j1) = lenc(i)%loc_lvdotr(:)
               k_plplenc(1,j0:j1) = i
               k_plplenc(2,j0:j1) = lenc(i)%index2(:)
               j0 = j1 + 1
            end if
         end do
      end if

      return
   end subroutine symba_encounter_check_all_triangular


   module function symba_encounter_check_pl(self, param, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between massive bodies.
      !!
      implicit none
      ! Arguments
      class(symba_pl),            intent(inout)  :: self   !! SyMBA test particle object  
      class(swiftest_parameters), intent(in)     :: param  !! Current swiftest run configuration parameters
      class(symba_nbody_system),  intent(inout)  :: system !! SyMBA nbody system object
      real(DP),                   intent(in)     :: dt     !! step size
      integer(I4B),               intent(in)     :: irec   !! Current recursion level
      ! Result
      logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      ! Internals
      integer(I8B) :: k, nplplm, kenc
      integer(I4B) :: i, j, nenc, npl
      logical, dimension(:), allocatable :: lencounter, loc_lvdotr, lvdotr
      integer(I4B), dimension(:), allocatable :: index1, index2
  
      if (self%nbody == 0) return

      associate(pl => self)
         nplplm = pl%nplplm
         npl = pl%nbody
         allocate(lencounter(nplplm))
         allocate(loc_lvdotr(nplplm))
  
         call symba_encounter_check_all_flat(nplplm, pl%k_plpl, pl%xh, pl%vh, pl%rhill, dt, irec, lencounter, loc_lvdotr)

         nenc = count(lencounter(:))

         lany_encounter = nenc > 0
         if (lany_encounter) then 
            associate(plplenc_list => system%plplenc_list)
               call plplenc_list%resize(nenc)
               allocate(lvdotr(nenc))					
               allocate(index1(nenc))					
               allocate(index2(nenc))					
               lvdotr(:) = pack(loc_lvdotr(:), lencounter(:))					
               index1(:) = pack(pl%k_plpl(1,1:nplplm), lencounter(:)) 					
               index2(:) = pack(pl%k_plpl(2,1:nplplm), lencounter(:)) 					
               deallocate(lencounter, loc_lvdotr)					
               call move_alloc(lvdotr, plplenc_list%lvdotr)					
               call move_alloc(index1, plplenc_list%index1) 					
               call move_alloc(index2, plplenc_list%index2) 					
               do k = 1, nenc
                  i = plplenc_list%index1(k)
                  j = plplenc_list%index2(k)
                  call util_index_eucl_ij_to_k(npl, i, j, kenc)
                  plplenc_list%kidx(k) = kenc
                  plplenc_list%id1(k) = pl%id(plplenc_list%index1(k))					
                  plplenc_list%id2(k) = pl%id(plplenc_list%index2(k))					
                  plplenc_list%status(k) = ACTIVE
                  plplenc_list%level(k) = irec
                  pl%lencounter(i) = .true.
                  pl%lencounter(j) = .true.
                  pl%levelg(i) = irec
                  pl%levelm(i) = irec
                  pl%levelg(j) = irec
                  pl%levelm(j) = irec
                  pl%nplenc(i) = pl%nplenc(i) + 1
                  pl%nplenc(j) = pl%nplenc(j) + 1
               end do
            end associate
         end if
      end associate
      return
   end function symba_encounter_check_pl


   module function symba_encounter_check(self, param, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between test particles and massive bodies in the pltpenc list.
      !! Note: This method works for the polymorphic symba_pltpenc and symba_plplenc types.
      !!
      !! Adapted from portions of David E. Kaufmann's Swifter routine: symba_step_recur.f90
      implicit none
      ! Arguments
      class(symba_encounter),     intent(inout) :: self           !! SyMBA pl-pl encounter list object
      class(swiftest_parameters), intent(in)    :: param          !! Current swiftest run configuration parameters
      class(symba_nbody_system),  intent(inout) :: system         !! SyMBA nbody system object
      real(DP),                   intent(in)    :: dt             !! step size
      integer(I4B),               intent(in)    :: irec           !! Current recursion level 
      logical                                   :: lany_encounter !! Returns true if there is at least one close encounter  
      ! Internals
      integer(I4B)              :: i, j, k, lidx, nenc_enc
      real(DP), dimension(NDIM) :: xr, vr
      logical                   :: lencounter, isplpl
      real(DP)                  :: rlim2, rji2
      logical, dimension(:), allocatable :: lencmask
      integer(I4B), dimension(:), allocatable :: encidx

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
            nenc_enc = count(lencmask(:))
            if (nenc_enc == 0) return

            allocate(encidx(nenc_enc))
            encidx(:) = pack([(k, k = 1, self%nenc)], lencmask(:))

            do lidx = 1, nenc_enc
               k = encidx(lidx)
               i = self%index1(k)
               j = self%index2(k)
               if (isplpl) then
                  xr(:) = pl%xh(:,j) - pl%xh(:,i)
                  vr(:) = pl%vb(:,j) - pl%vb(:,i)
                  call symba_encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%rhill(i), pl%rhill(j), dt, irec, lencounter, self%lvdotr(k))
               else
                  xr(:) = tp%xh(:,j) - pl%xh(:,i)
                  vr(:) = tp%vb(:,j) - pl%vb(:,i)
                  call symba_encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%rhill(i), 0.0_DP, dt, irec, lencounter, self%lvdotr(k))
               end if
               if (lencounter) then
                  if (isplpl) then
                     rlim2 = (pl%radius(i) + pl%radius(j))**2
                  else
                     rlim2 = (pl%radius(i))**2
                  end if
                  rji2 = dot_product(xr(:), xr(:))! Check to see if these are physically overlapping bodies first, which we should ignore
                  if (rji2 > rlim2) then
                     lany_encounter = .true.
                     pl%levelg(i) = irec
                     pl%levelm(i) = MAX(irec, pl%levelm(i))
                     if (isplpl) then
                        pl%levelg(j) = irec
                        pl%levelm(j) = MAX(irec, pl%levelm(j))
                     else
                        tp%levelg(j) = irec
                        tp%levelm(j) = MAX(irec, tp%levelm(j))
                     end if
                     self%level(k) = irec
                  end if
               end if   
            end do
         end select
      end select

      return
   end function symba_encounter_check


   module function symba_encounter_check_tp(self, param, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between test particles and massive bodies.
      !!
      implicit none
      ! Arguments
      class(symba_tp),            intent(inout) :: self   !! SyMBA test particle object  
      class(swiftest_parameters), intent(in)    :: param  !! Current swiftest run configuration parameters
      class(symba_nbody_system),  intent(inout) :: system !! SyMBA nbody system object
      real(DP),                   intent(in)    :: dt     !! step size
      integer(I4B),               intent(in)    :: irec   !! Current recursion level
      ! Result
      logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      ! Internals
      real(DP)                                  :: r2crit, vdotr, r2, v2, tmin, r2min, term2
      integer(I4B)                              :: i, j, k,nenc, plind, tpind
      real(DP),     dimension(NDIM)             :: xr, vr
      logical,      dimension(:,:), allocatable :: lencounter, loc_lvdotr
  
      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody, pl => system%pl, npl => system%pl%nbody)
         allocate(lencounter(ntp, npl), loc_lvdotr(ntp, npl))
         lencounter(:,:) = .false.
   
         do j = 1, npl
            do i = 1, ntp
               xr(1) = tp%xh(1, i) - pl%xh(1, j)
               xr(2) = tp%xh(2, i) - pl%xh(2, j)
               xr(3) = tp%xh(3, i) - pl%xh(3, j)
               vr(1) = tp%vh(1, i) - pl%vh(1, j)
               vr(2) = tp%vh(2, i) - pl%vh(2, j)
               vr(3) = tp%vh(3, i) - pl%vh(3, j)
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
               pltpenc_list%lvdotr(1:nenc) = pack(loc_lvdotr(1:ntp, 1:npl), lencounter(1:ntp, 1:npl))
               pltpenc_list%index1(1:nenc) = pack(spread([(i, i = 1, npl)], dim=1, ncopies=ntp), lencounter(1:ntp, 1:npl)) 
               pltpenc_list%index2(1:nenc) = pack(spread([(i, i = 1, ntp)], dim=2, ncopies=npl), lencounter(1:ntp, 1:npl))
               pltpenc_list%id1(1:nenc) = pl%id(pltpenc_list%index1(1:nenc))
               pltpenc_list%id2(1:nenc) = tp%id(pltpenc_list%index2(1:nenc))
               select type(pl)
               class is (symba_pl)
                  pl%lencounter(1:npl) = .false.
                  do k = 1, nenc
                     plind = pltpenc_list%index1(k)
                     tpind = pltpenc_list%index2(k)
                     pl%lencounter(plind) = .true.
                     pl%levelg(plind) = irec
                     pl%levelm(plind) = irec
                     tp%levelg(tpind) = irec
                     tp%levelm(tpind) = irec
                     pl%ntpenc(plind) = pl%ntpenc(plind) + 1
                     tp%nplenc(tpind) = tp%nplenc(tpind) + 1
                  end do
               end select
            end associate
         end if
      end associate

      return
   end function symba_encounter_check_tp


   pure subroutine symba_encounter_check_one(xr, yr, zr, vxr, vyr, vzr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
      !$omp declare simd(symba_encounter_check_one)
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
      real(DP)     :: r2crit, rshell_irec
      integer(I4B) :: i

      rshell_irec = 1._DP
      do i = 1, irec
         rshell_irec = rshell_irec * RSHELL
      end do
      r2crit = (rhill1 + rhill2) * RHSCALE * rshell_irec
      r2crit = r2crit**2
      call rmvs_chk_ind(xr, yr, zr, vxr, vyr, vzr, dt, r2crit, lencounter, lvdotr)

      return
   end subroutine symba_encounter_check_one

end submodule s_symba_encounter_check
