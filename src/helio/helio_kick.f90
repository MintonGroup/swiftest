submodule(helio_classes) s_helio_kick
   use swiftest
contains
   module subroutine helio_kickvb_pl(self, dt)
      !! author: David A. Minton
      !!
      !! Kick barycentric velocities of bodies
      !!
      !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f 
      !! Adapted from David E. Kaufmann's Swifter routine helio_kickvb.f90
      implicit none
      ! Arguments
      class(helio_pl), intent(inout) :: self !! Swiftest generic body object
      real(DP),        intent(in)    :: dt   !! Stepsize
      ! Internals
      integer(I4B) :: i

      associate(pl => self, npl => self%nbody)
         if (npl ==0) return
         do concurrent(i = 1:npl, pl%status(i) == ACTIVE) 
            pl%vb(:, i) = pl%vb(:, i) + pl%ah(:, i) * dt
         end do
      end associate
   
      return
   
   end subroutine helio_kickvb_pl

   module subroutine helio_kickvb_tp(self, dt)
      !! author: David A. Minton
      !!
      !! Kick barycentric velocities of bodies
      !!
      !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_kickvb_tp.f90
      implicit none
      ! Arguments
      class(helio_tp), intent(inout) :: self !! Swiftest generic body object
      real(DP),        intent(in)    :: dt   !! Stepsize
      ! Internals
      integer(I4B) :: i

      associate(tp => self, ntp => self%nbody)
         if (ntp ==0) return
         do concurrent(i = 1:ntp, tp%status(i) == ACTIVE) 
            tp%vb(:, i) = tp%vb(:, i) + tp%ah(:, i) * dt
         end do
      end associate
   
      return
   
   end subroutine helio_kickvb_tp
end submodule s_helio_kick