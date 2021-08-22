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
      ! Internals
      integer(I8B)                      :: k, nplpl
      real(DP)                          :: rji2, rlim2
      real(DP)                          :: dx, dy, dz
      integer(I4B) :: i, j
      real(DP), dimension(:,:), pointer :: ah, xh
      real(DP), dimension(NDIM,self%nbody) :: ahi, ahj
      integer(I4B), dimension(:,:), pointer :: k_plpl
      logical, dimension(:), pointer :: lmask
      real(DP), dimension(:), pointer :: Gmass

      associate(ah => self%ah, xh => self%xh, k_plpl => self%k_plpl, lmask => self%lmask, Gmass => self%Gmass)
         nplpl = self%nplpl
         ahi(:,:) = 0.0_DP
         ahj(:,:) = 0.0_DP
         !$omp parallel do default(shared)&
         !$omp private(k, i, j, dx, dy, dz, rji2)  &
         !$omp reduction(+:ahi) &
         !$omp reduction(-:ahj) 
         do k = 1_I8B, nplpl
            i = k_plpl(1,k)
            j = k_plpl(2,k)
            dx = xh(1, j) - xh(1, i)
            dy = xh(2, j) - xh(2, i)
            dz = xh(3, j) - xh(3, i)
            rji2 = dx**2 + dy**2 + dz**2
            if (lmask(i) .and. lmask(j)) call kick_getacch_int_one_pl(rji2, dx, dy, dz, Gmass(i), Gmass(j), ahi(1,i), ahi(2,i), ahi(3,i), ahj(1,j), ahj(2,j), ahj(3,j))
         end do
         !$omp end parallel do
         ah(:,:) = ah(:,:) + ahi(:,:) + ahj(:,:)
      end associate

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
      !real(DP), dimension(:,:), allocatable :: aht

      if ((self%nbody == 0) .or. (npl == 0)) return

      associate(tp => self, ntp => self%nbody)
         do i = 1, ntp
            if (tp%lmask(i)) then
               do j = 1, npl
                  block
                     real(DP) :: rji2
                     real(DP) :: dx, dy, dz
                     dx = tp%xh(1, i) - xhp(1, j)
                     dy = tp%xh(2, i) - xhp(1, j)
                     dz = tp%xh(3, i) - xhp(1, j)
                     rji2 = dx**2 + dy**2 + dz**2
                     call kick_getacch_int_one_tp(rji2, dx, dy, dz, GMpl(i), tp%ah(1,i), tp%ah(2,i), tp%ah(3,i))
                  end block
               end do
            end if
         end do
         !call move_alloc(aht, tp%ah)
      end associate
      
      return
   end subroutine kick_getacch_int_tp

   module pure elemental subroutine kick_getacch_int_one_pl(rji2, dx, dy, dz, Gmi, Gmj, axi, ayi, azi, axj, ayj, azj)
      implicit none
      real(DP), intent(in)  :: rji2
      real(DP), intent(in)  :: dx, dy, dz
      real(DP), intent(in)  :: Gmi
      real(DP), intent(in)  :: Gmj
      real(DP), intent(inout) :: axi, ayi, azi
      real(DP), intent(inout) :: axj, ayj, azj
      ! Internals
      real(DP) :: faci, facj, irij3

      irij3 = 1.0_DP / (rji2 * sqrt(rji2))
      faci = Gmi * irij3
      facj = Gmj * irij3
      axi = axi + facj * dx
      ayi = ayi + facj * dy
      azi = azi + facj * dz
      axj = axj - faci * dx
      ayj = ayj - faci * dy
      azj = azj - faci * dz
      return
   end subroutine kick_getacch_int_one_pl


   !module pure elemental subroutine kick_getacch_int_one_tp(rji2, dx, dy, dz, GMpl, ax, ay, az)
   module pure subroutine kick_getacch_int_one_tp(rji2, dx, dy, dz, GMpl, ax, ay, az)
      implicit none
      real(DP), intent(in)  :: rji2
      real(DP), intent(in)  :: dx, dy, dz
      real(DP), intent(in)  :: GMpl
      real(DP), intent(inout) :: ax, ay, az
      ! Internals
      real(DP) :: fac

      fac = GMpl / (rji2 * sqrt(rji2))
      ax = ax - fac * dx
      ay = ay - fac * dy
      az = az - fac * dz
      return
   end subroutine kick_getacch_int_one_tp

end submodule s_kick
