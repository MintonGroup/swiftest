submodule(swiftest_classes) s_kick
   use swiftest
contains

   module subroutine kick_getacch_int_pl(self, param)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of massive bodies
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int.f90
      implicit none
      ! Arguments
      class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
      class(swiftest_parameters), intent(inout) :: param !! Current swiftest run configuration parameters
      ! Internals
      type(interaction_timer), save :: itimer
      logical, save :: lfirst = .true.

      if (param%ladaptive_interactions) then
         if (self%nplpl > 0) then
            if (lfirst) then
               write(itimer%loopname, *) "kick_getacch_int_pl"
               write(itimer%looptype, *) "INTERACTION"
               call itimer%time_this_loop(param, self%nplpl, self)
               lfirst = .false.
            else
               if (itimer%check(param, self%nplpl)) call itimer%time_this_loop(param, self%nplpl, self)
            end if
         else
            param%lflatten_interactions = .false.
         end if
      end if

      if (param%lflatten_interactions) then
         call kick_getacch_int_all_flat_pl(self%nbody, self%nplpl, self%k_plpl, self%xh, self%Gmass, self%radius, self%ah)
      else
         call kick_getacch_int_all_triangular_pl(self%nbody, self%nbody, self%xh, self%Gmass, self%radius, self%ah)
      end if

      if (param%ladaptive_interactions .and. self%nplpl > 0) then 
         if (itimer%is_on) call itimer%adapt(param, self%nplpl, self)
      end if

      return
   end subroutine kick_getacch_int_pl


   module subroutine kick_getacch_int_tp(self, param, GMpl, xhp, npl)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of test particles by massive bodies
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int_tp.f90
      implicit none
      ! Arguments
      class(swiftest_tp),         intent(inout) :: self  !! Swiftest test particle object
      class(swiftest_parameters), intent(inout) :: param !! Current swiftest run configuration parameters
      real(DP), dimension(:),     intent(in)    :: GMpl  !! Massive body masses
      real(DP), dimension(:,:),   intent(in)    :: xhp   !! Massive body position vectors
      integer(I4B),               intent(in)    :: npl   !! Number of active massive bodies

      if ((self%nbody == 0) .or. (npl == 0)) return

      call kick_getacch_int_all_tp(self%nbody, npl, self%xh, xhp, GMpl, self%lmask, self%ah)
      
      return
   end subroutine kick_getacch_int_tp


   module subroutine kick_getacch_int_all_flat_pl(npl, nplpl, k_plpl, x, Gmass, radius, acc)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations for massive bodies, with parallelization.
      !! This is the flattened (single loop) version that uses the k_plpl interaction pair index array
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int.f9
      implicit none
      integer(I4B),                 intent(in)    :: npl    !! Number of massive bodies
      integer(I8B),                 intent(in)    :: nplpl  !! Number of massive body interactions to compute
      integer(I4B), dimension(:,:), intent(in)    :: k_plpl !! Array of interaction pair indices (flattened upper triangular matrix)
      real(DP),     dimension(:,:), intent(in)    :: x      !! Position vector array
      real(DP),     dimension(:),   intent(in)    :: Gmass  !! Array of massive body G*mass
      real(DP),     dimension(:),   intent(in)    :: radius !! Array of massive body radii
      real(DP),     dimension(:,:), intent(inout) :: acc    !! Acceleration vector array 
      ! Internals
      integer(I8B)                      :: k
      real(DP), dimension(NDIM,npl) :: ahi, ahj
      integer(I4B) :: i, j
      real(DP)     :: rji2, rlim2
      real(DP)     :: xr, yr, zr

      ahi(:,:) = 0.0_DP
      ahj(:,:) = 0.0_DP

      !$omp parallel do default(private) schedule(static)&
      !$omp shared(nplpl, k_plpl, x, Gmass, radius) &
      !$omp lastprivate(rji2, rlim2, xr, yr, zr) &
      !$omp reduction(+:ahi) &
      !$omp reduction(-:ahj) 
      do k = 1_I8B, nplpl
         i = k_plpl(1, k)
         j = k_plpl(2, k)
         xr = x(1, j) - x(1, i) 
         yr = x(2, j) - x(2, i) 
         zr = x(3, j) - x(3, i) 
         rji2 = xr**2 + yr**2 + zr**2
         rlim2 = (radius(i) + radius(j))**2
         if (rji2 > rlim2) call kick_getacch_int_one_pl(rji2, xr, yr, zr, Gmass(i), Gmass(j), &
                                 ahi(1,i), ahi(2,i), ahi(3,i), ahj(1,j), ahj(2,j), ahj(3,j))
      end do
      !$omp end parallel do 
     
      do concurrent(i = 1:npl)
         acc(:,i) = acc(:,i) + ahi(:,i) + ahj(:,i)
      end do

      return
   end subroutine kick_getacch_int_all_flat_pl


   module subroutine kick_getacch_int_all_triangular_pl(npl, nplm, x, Gmass, radius, acc)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations for massive bodies, with parallelization.
      !! This is the upper triangular matrix (double loop) version.
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int.f9
      implicit none
      integer(I4B),                 intent(in)    :: npl    !! Total number of massive bodies
      integer(I4B),                 intent(in)    :: nplm   !! Number of fully interacting massive bodies 
      real(DP),     dimension(:,:), intent(in)    :: x      !! Position vector array
      real(DP),     dimension(:),   intent(in)    :: Gmass  !! Array of massive body G*mass
      real(DP),     dimension(:),   intent(in)    :: radius !! Array of massive body radii
      real(DP),     dimension(:,:), intent(inout) :: acc    !! Acceleration vector array 
      ! Internals
      real(DP), dimension(NDIM,npl) :: ahi, ahj
      integer(I4B) :: i, j
      real(DP)     :: rji2, rlim2
      real(DP)     :: xr, yr, zr

      ahi(:,:) = 0.0_DP
      ahj(:,:) = 0.0_DP

      !$omp parallel do default(private) schedule(static)&
      !$omp shared(npl, nplm, x, Gmass, radius) &
      !$omp lastprivate(rji2, rlim2, xr, yr, zr) &
      !$omp reduction(+:ahi) &
      !$omp reduction(-:ahj) 
      do i = 1, nplm
         do concurrent(j = i+1:npl)
            xr = x(1, j) - x(1, i) 
            yr = x(2, j) - x(2, i) 
            zr = x(3, j) - x(3, i) 
            rji2 = xr**2 + yr**2 + zr**2
            rlim2 = (radius(i) + radius(j))**2
            if (rji2 > rlim2) call kick_getacch_int_one_pl(rji2, xr, yr, zr, Gmass(i), Gmass(j), &
                                    ahi(1,i), ahi(2,i), ahi(3,i), ahj(1,j), ahj(2,j), ahj(3,j))
         end do
      end do
      !$omp end parallel do

      do concurrent(i = 1:npl)
         acc(:,i) = acc(:,i) + ahi(:,i) + ahj(:,i)
      end do

      return
   end subroutine kick_getacch_int_all_triangular_pl


   module subroutine kick_getacch_int_all_tp(ntp, npl, xtp, xpl, GMpl, lmask, acc)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of test particles by massive bodies with parallelisim
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int_tp.f99
      implicit none
      integer(I4B),                 intent(in)    :: ntp    !! Number of test particles
      integer(I4B),                 intent(in)    :: npl    !! Number of massive bodies
      real(DP),     dimension(:,:), intent(in)    :: xtp    !! Test particle position vector array
      real(DP),     dimension(:,:), intent(in)    :: xpl    !! Massive body particle position vector array
      real(DP),     dimension(:),   intent(in)    :: GMpl   !! Array of massive body G*mass
      logical,      dimension(:),   intent(in)    :: lmask  !! Logical mask indicating which test particles should be computed
      real(DP),     dimension(:,:), intent(inout) :: acc    !! Acceleration vector array 
      ! Internals
      real(DP)     :: rji2
      real(DP)     :: xr, yr, zr
      integer(I4B) :: i, j

      !$omp parallel do default(private) schedule(static)&
      !$omp shared(npl, ntp, lmask, xtp, xpl, GMpl) &
      !$omp reduction(-:acc)
      do i = 1, ntp
         if (lmask(i)) then
            do j = 1, npl
               xr = xtp(1, i) - xpl(1, j)
               yr = xtp(2, i) - xpl(2, j)
               zr = xtp(3, i) - xpl(3, j)
               rji2 = xr**2 + yr**2 + zr**2
               call kick_getacch_int_one_tp(rji2, xr, yr, zr, GMpl(j), acc(1,i), acc(2,i), acc(3,i))
            end do
         end if
      end do
      !$omp end parallel do
      
      return
   end subroutine kick_getacch_int_all_tp


   module pure subroutine kick_getacch_int_one_pl(rji2, xr, yr, zr, Gmi, Gmj, axi, ayi, azi, axj, ayj, azj)
      !!$omp declare simd(kick_getacch_int_one_pl)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations for a single pair of massive bodies
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int.f9
      implicit none
      real(DP), intent(in)  :: rji2            !! Square of distance between the two bodies
      real(DP), intent(in)  :: xr, yr, zr      !! Distances between the two bodies in x, y, and z directions
      real(DP), intent(in)  :: Gmi             !! G*mass of body i
      real(DP), intent(in)  :: Gmj             !! G*mass of body j
      real(DP), intent(inout) :: axi, ayi, azi !! Acceleration vector components of body i
      real(DP), intent(inout) :: axj, ayj, azj !! Acceleration vector components of body j
      ! Internals
      real(DP) :: faci, facj, irij3

      irij3 = 1.0_DP / (rji2 * sqrt(rji2))
      faci = Gmi * irij3
      facj = Gmj * irij3
      axi = axi + facj * xr
      ayi = ayi + facj * yr
      azi = azi + facj * zr
      axj = axj - faci * xr
      ayj = ayj - faci * yr
      azj = azj - faci * zr

      return
   end subroutine kick_getacch_int_one_pl


   module pure subroutine kick_getacch_int_one_tp(rji2, xr, yr, zr, GMpl, ax, ay, az)
      !!$omp declare simd(kick_getacch_int_one_tp)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of a single test particle massive body pair.
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int_tp.f90
      implicit none
      real(DP), intent(in)  :: rji2         !! Square of distance between the test particle and massive body
      real(DP), intent(in)  :: xr, yr, zr   !! Distances between the two bodies in x, y, and z directions
      real(DP), intent(in)  :: Gmpl         !! G*mass of massive body
      real(DP), intent(inout) :: ax, ay, az !! Acceleration vector components of test particle
      ! Internals
      real(DP) :: fac

      fac = GMpl * sqrt(1.0_DP / (rji2*rji2*rji2))
      ax = ax - fac * xr
      ay = ay - fac * yr
      az = az - fac * zr

      return
   end subroutine kick_getacch_int_one_tp

end submodule s_kick
