submodule(swiftest_classes) s_kick
   use swiftest
contains

   module subroutine kick_getacch_int_pl(self)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of massive bodies
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int.f90
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self

      call kick_getacch_int_all_pl(self%nbody, self%nplpl, self%k_plpl, self%xh, self%Gmass, self%radius, self%ah)

      return
   end subroutine kick_getacch_int_pl


   module pure subroutine kick_getacch_int_tp(self, GMpl, xhp, npl)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of test particles by massive bodies
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int_tp.f90
      implicit none
      ! Arguments
      class(swiftest_tp),       intent(inout) :: self !! Swiftest test particle
      real(DP), dimension(:),   intent(in)    :: GMpl !! Massive body masses
      real(DP), dimension(:,:), intent(in)    :: xhp  !! Massive body position vectors
      integer(I4B),             intent(in)    :: npl  !! Number of active massive bodies
      ! Internals
      integer(I4B)              :: i, j

      if ((self%nbody == 0) .or. (npl == 0)) return

      associate(tp => self, ntp => self%nbody)
         do i = 1, ntp
            if (tp%lmask(i)) then
               do j = 1, npl
                  block
                     real(DP) :: rji2
                     real(DP) :: xr, yr, zr
                     xr = tp%xh(1, i) - xhp(1, j)
                     yr = tp%xh(2, i) - xhp(1, j)
                     zr = tp%xh(3, i) - xhp(1, j)
                     rji2 = xr**2 + yr**2 + zr**2
                     call kick_getacch_int_one_tp(rji2, xr, yr, zr, GMpl(i), tp%ah(1,i), tp%ah(2,i), tp%ah(3,i))
                  end block
               end do
            end if
         end do
      end associate
      
      return
   end subroutine kick_getacch_int_tp


   module subroutine kick_getacch_int_all_pl(npl, nplpl, k_plpl, x, Gmass, radius, acc)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations for massive bodies, with parallelization
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
      real(DP)                          :: rji2, rlim2
      real(DP)                          :: xr, yr, zr
      integer(I4B) :: i, j
      real(DP), dimension(NDIM,npl) :: ahi, ahj

      ahi(:,:) = 0.0_DP
      ahj(:,:) = 0.0_DP
      !$omp parallel do default(private)&
      !$omp shared(nplpl, k_plpl, x, Gmass, radius)  &
      !$omp reduction(+:ahi) &
      !$omp reduction(-:ahj) 
      do k = 1_I8B, nplpl
         i = k_plpl(1,k)
         j = k_plpl(2,k)
         xr = x(1, j) - x(1, i)
         yr = x(2, j) - x(2, i)
         zr = x(3, j) - x(3, i)
         rji2 = xr**2 + yr**2 + zr**2
         rlim2 = (radius(i) + radius(j))**2
         if (rji2 > rlim2) call kick_getacch_int_one_pl(rji2, xr, yr, zr, Gmass(i), Gmass(j), ahi(1,i), ahi(2,i), ahi(3,i), ahj(1,j), ahj(2,j), ahj(3,j))
      end do
      !$omp end parallel do
      !$omp parallel workshare
      acc(:,1:npl) = acc(:,1:npl) + ahi(:,1:npl) + ahj(:,1:npl)
      !$omp end parallel workshare
      return
   end subroutine kick_getacch_int_all_pl


   module pure subroutine kick_getacch_int_one_pl(rji2, xr, yr, zr, Gmi, Gmj, axi, ayi, azi, axj, ayj, azj)
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

      fac = GMpl / (rji2 * sqrt(rji2))
      ax = ax - fac * xr
      ay = ay - fac * yr
      az = az - fac * zr
      return
   end subroutine kick_getacch_int_one_tp

end submodule s_kick
