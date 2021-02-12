submodule(swiftest_classes) kick_implementations
contains
   module procedure kick_vh_body 
   !! author: David A. Minton
   !!
   !! Kick heliocentric velocities of bodies
   !!
   !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f and kickvh_tp.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_kickvh.f90 and whm_kickvh_tp.f90
   use swiftest
   implicit none
   integer(I4B) :: i

   associate(n => self%nbody, vh => self%vh, ah => self%ah, status => self%status)
      if (n == 0) return
      do concurrent(i = 1:n, status(i) == ACTIVE) 
      !do i = 1, n
         vh(:, i) = vh(:, i) + ah(:, i) * dt
      end do
   end associate

   return

   end procedure kick_vh_body

   module procedure kick_vb_body 
      !! author: David A. Minton
      !!
      !! Kick barycentric velocities of bodies
      !!
      !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f and kickvh_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_kickvb.f90 and helio_kickvb_tp.f90
      use swiftest
      implicit none
      integer(I4B) :: i

      associate(n => self%nbody, vb => self%vb, ah => self%ah, status => self%status)
         if (n ==0) return
         !do concurrent(i = 1:n, status(i) == ACTIVE) 
         do i = 1, n
            vb(:, i) = vb(:, i) + ah(:, i) * dt
         end do
      end associate
   
      return
   
      end procedure kick_vb_body
end submodule kick_implementations
