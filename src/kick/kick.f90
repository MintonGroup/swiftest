submodule(swiftest_classes) s_kick
   use swiftest
contains

   module pure subroutine kick_getacch_int_pl(self)
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
      integer(I4B)                      :: k
      real(DP)                          :: rji2, irij3, faci, facj
      real(DP), dimension(NDIM)         :: dx

      associate(pl => self, npl => self%nbody, nplpl => self%nplpl)
         do k = 1, nplpl
            associate(i => pl%k_plpl(1, k), j => pl%k_plpl(2, k))
               if (pl%lmask(i) .and. pl%lmask(j)) then
                  dx(:) = pl%xh(:, j) - pl%xh(:, i)
                  rji2  = dot_product(dx(:), dx(:))
                  irij3 = 1.0_DP / (rji2 * sqrt(rji2))
                  faci = pl%Gmass(i) * irij3
                  facj = pl%Gmass(j) * irij3
                  pl%ah(:, i) = pl%ah(:, i) + facj * dx(:)
                  pl%ah(:, j) = pl%ah(:, j) - faci * dx(:)
               end if
            end associate
         end do
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
      real(DP)                  :: rji2, irij3, fac, r2
      real(DP), dimension(NDIM) :: dx

      associate(tp => self, ntp => self%nbody)
         do concurrent(i = 1:ntp, tp%lmask(i))
            do j = 1, npl
               dx(:) = tp%xh(:,i) - xhp(:, j)
               r2 = dot_product(dx(:), dx(:))
               fac = GMpl(j) / (r2 * sqrt(r2))
               tp%ah(:, i) = tp%ah(:, i) - fac * dx(:)
            end do
         end do
      end associate
      
      return
   end subroutine kick_getacch_int_tp

end submodule s_kick
