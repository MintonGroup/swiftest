submodule(swiftest_classes) s_kick
   use swiftest
contains
   module subroutine kick_vh_body(self, dt)
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
   end subroutine kick_vh_body

   module subroutine kick_vb_body(self, dt)
      !! author: David A. Minton
      !!
      !! Kick barycentric velocities of bodies
      !!
      !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f and kickvh_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_kickvb.f90 and helio_kickvb_tp.f90
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self !! Swiftest generic body object
      real(DP),                     intent(in)    :: dt   !! Stepsize
      ! Internals
      integer(I4B) :: i

      associate(n => self%nbody, vb => self%vb, ah => self%ah, status => self%status)
         if (n ==0) return
         do concurrent(i = 1:n, status(i) == ACTIVE) 
            vb(:, i) = vb(:, i) + ah(:, i) * dt
         end do
      end associate
   
      return
   
   end subroutine kick_vb_body
end submodule s_kick
