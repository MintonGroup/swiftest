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

   associate(n => self%nbody, vh => self%vh, ah => self%ah, status => self%status)
      if (n == 0) return
      where(status(1:n) == ACTIVE)
         vh(1, :) = vh(1, :) + ah(1, :) * dt
         vh(2, :) = vh(2, :) + ah(2, :) * dt
         vh(3, :) = vh(3, :) + ah(3, :) * dt
      end where
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
   
      associate(n => self%nbody, vb => self%vb, ah => self%ah, status => self%status)
         if (n ==0) return
         where(status(1:n) == ACTIVE)
            vb(1, :) = vb(1, :) + ah(1, :) * dt
            vb(2, :) = vb(2, :) + ah(2, :) * dt
            vb(3, :) = vb(3, :) + ah(3, :) * dt
         end where
      end associate
   
      return
   
      end procedure kick_vb_body
end submodule kick_implementations
