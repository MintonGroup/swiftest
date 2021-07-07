submodule(swiftest_classes) s_kick
   use swiftest
contains
   module subroutine kickvh_body(self, dt)
      !! author: David A. Minton
      !!
      !! Kick heliocentric velocities of bodies
      !!
      !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f and kickvh_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kickvh.f90 and whm_kickvh_tp.f90
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self !! Swiftest generic body object
      real(DP),                     intent(in)    :: dt   !! Stepsize
      ! Internals
      integer(I4B) :: i

      associate(n => self%nbody, vh => self%vh, ah => self%ah, status => self%status)
         if (n == 0) return
         do i = 1, n
            if (status(i) == ACTIVE) vh(:, i) = vh(:, i) + ah(:, i) * dt
         end do
      end associate

      return
   end subroutine kickvh_body

end submodule s_kick
